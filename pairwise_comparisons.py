#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 14:46:03 2021

@author: Dr Charles Foster, Virology Research Lab

To do:
- Fix up or remove the amino acid translation option:
    * Currently weird/wrong when deletions occur
    * Not really useful for this program?

- Add option to _only_ get SNP matrix as output

- Add option to work with _all_ seqs in an alignment?
    * Would eliminate the need for a comparisons.tsv file

- Update README to require networkx installation, or formulate a different way
  of converting to matrix & writing to csv
  
- Update README in general
    
"""
import pandas as pd
from itertools import groupby, count
import networkx as nx
from Bio import AlignIO
from Bio.Seq import Seq
import argparse
import re
import multiprocessing as mp
import sys
try:
    import tqdm
    tqdm_installed = True
except:
    tqdm_installed = False
version = '1.1.0'
#%%
def printc(thing, level):
    '''
    Print in colour :)
    '''
    cols = {'green':'\033[1;32m', 'blue':'\033[96m', 'yellow':'\033[93m'}
    col = cols[level]
    print(f"{col}{thing}\033[0m")
    return()

#%%
def check_input(comparisons, MSA):
    '''
    Quick check to make sure that the comparison names match names in the alignment
    '''
    flat_comp = [x.split(sep='\t') for x in comparisons]
    flat_comp = [item for sublist in flat_comp for item in sublist]
    flat_comp = [x.split()[0] for x in flat_comp]
    seqs = [x.id for x in MSA]
    if not set(flat_comp).issubset(seqs):
        sys.exit('''
Error: Some of your comparison IDs are not present as headers in your alignment file. Check the spelling and try again
''')

#%%
def intervals(data):
    out = []
    counter = count()
    for key, group in groupby(data, key = lambda x: x-next(counter)):
        block = list(group)
        out.append([block[0], block[-1], len(block)])
    return out

#%%
def find_consecutive_gaps(reduced_differences):
    indels = reduced_differences[reduced_differences['vartype']=='INDEL']['summary']
    s1_indels = [int(re.sub(r"\D", "",x)) for x in indels if x.startswith('-')]
    s2_indels = [int(re.sub(r"\D", "",x)) for x in indels if x.endswith('-')]
    s1_intervals = intervals(s1_indels)
    s2_intervals = intervals(s2_indels)
    return {"seq1":s1_intervals, "seq2":s2_intervals}

#%%
def convert_to_matrix(results_table):
    '''
    Convert results table into a SNP matrix
    '''
    reduced_results = results_table[["Seq1","Seq2","Num_Differences"]]
    G = nx.from_pandas_edgelist(
    reduced_results,
    source='Seq1',
    target='Seq2',
    edge_attr='Num_Differences')

    adjacency_df = pd.DataFrame(
        nx.adjacency_matrix(G, weight='Num_Differences').todense(),
        index=G.nodes,
        columns=G.nodes
    )
    return adjacency_df

#%%
def find_insertions(MSA, annotations):
    refseq = []
    for record in MSA:
        if record.id == annotations.iloc[0][0]:
            refseq.append(record)
    insertions = [i for i, ltr in enumerate(refseq[0].seq) if ltr == '-' and i < 29902] # need to change to args.reflength
    return insertions

#%%
def fix_bed_insertions(annotations, insertions):
    for i in annotations.index:
        for site in insertions:
            if site < annotations.iloc[i][1] and site < (annotations.iloc[i][2]-1):
                annotations.at[i,1] = annotations.at[i,1] + 1
                annotations.at[i,2] = annotations.at[i,2] + 1
            elif site > annotations.iloc[i][1]  and site < (annotations.iloc[i][2]-1):
                annotations.at[i,2] = annotations.at[i,2] + 1
 
#%%
def indices(string,bad_sites=['N','-']):
    '''
    Find the index of the first and last non-N/non-gap bases in a DNA sequence.
    '''
    base_list = [x for x in string]
    non_N = [i for i,x in enumerate(base_list) if x not in bad_sites]
    start_end = [non_N[0],non_N[-1]]
    return start_end
#%%
def trim_seqs(seq1,seq2):
    idx1 = indices(seq1)
    idx2 = indices(seq2)
    if idx1[0] >= idx2[0]:
        start = idx1[0]
    else:
        start = idx2[0]
    if idx1[1] >= idx2[1]:
        end = idx2[1]+1
    else:
        end = idx1[1]+1
    new_seq1 = seq1[start:end]
    new_seq2 = seq2[start:end]
    return [new_seq1,new_seq2]

#%%
def prep_aa_seqs(fix_seq1, fix_seq2, new_site):
    '''
    Remove internal gaps before attempting translation
    '''
    deletions = []
    keep = []
    for idx,(i,j) in enumerate(zip(fix_seq1,fix_seq2)):
        if i == '-' or j == '-':
            deletions.append(idx)
        if i != '-' or j != '-':
            keep.append(idx)
    aa_seq1 = Seq(''.join([fix_seq1[x] for x in keep]))
    aa_seq2 = Seq(''.join([fix_seq2[x] for x in keep]))
    for delsite in deletions:
        if delsite < new_site:
            new_site = new_site-1
    res = [aa_seq1, aa_seq2, new_site]
    return res
 
#%%
def get_aa_differences(seq1,seq2,annotations,difference_sites):
    aa_results = []
    for i in range(len(annotations.index)):
        gene_range = range(int(annotations.iloc[i][1]),int(annotations.iloc[i][2]))
        reduced_df = annotations.iloc[i]
        start = int(reduced_df[1])
        end = int(reduced_df[2])
        gene_name = reduced_df.iloc[3]
        for site in difference_sites:
            if site in gene_range:
                new_site = site
                if reduced_df[4] == True:
                    slip_site = int(reduced_df[5])
                    fix_seq1a = seq1[start:slip_site]
                    fix_seq1b = seq1[slip_site-1:end]
                    fix_seq1 = fix_seq1a+fix_seq1b
                    fix_seq2a = seq2[start:slip_site]
                    fix_seq2b = seq2[slip_site-1:end]
                    fix_seq2 = fix_seq2a+fix_seq2b
                else:
                    fix_seq1 = seq1[start:end]
                    fix_seq2 = seq2[start:end]

                if reduced_df[4] == True and int(reduced_df[5]) < new_site:
                    new_site = new_site +1
                
                new_site = new_site-start
                new_seqs = prep_aa_seqs(fix_seq1, fix_seq2, new_site)
                aa_seq1 = new_seqs[0]
                aa_seq2 = new_seqs[1]
                new_site = new_seqs[2]
                 
                if len(aa_seq1) % 3 !=0:
                    aa_results.append(str(site+1)+':'+gene_name+":FRAMESHIFT")
                else:
                    mod = list(range(0,(len(aa_seq1)+1))).index(new_site+1) % 3
                    if mod == 1:
                        codon_pos = 1
                        slice_range = range(new_site,new_site+4)
                    elif mod == 2:
                        codon_pos = 2
                        slice_range = range(new_site-1,new_site+3)
                    else:
                        codon_pos = 3
                        slice_range = range(new_site-2,new_site+2)
                    codon_num = int(slice_range[-1]/3)
                    codon1 = aa_seq1[slice_range[0]:slice_range[-1]]
                    codon2 = aa_seq2[slice_range[0]:slice_range[-1]]
                    codon1_AA = codon1.translate()
                    codon2_AA = codon2.translate()
                    if codon1_AA == codon2_AA:
                        result = str(site+1)+':'+gene_name+':'+str(codon1_AA)+str(codon_num)+str(codon2_AA)+":s"
                    elif codon1_AA != codon2_AA:
                        result = str(site+1)+':'+gene_name+':'+str(codon1_AA)+str(codon_num)+str(codon2_AA)+":n"
                    aa_results.append(result)
    return aa_results

#%%
def calc_differences(name1,seq1,name2,seq2,removals, annotations, gap_identity):
    ''''
    Calculate total pairwise differences. These are later pruned based on trimming conditions.
    '''
    full_length = len(seq1) #aligned, so both seqs same length
    sites = []
    summary = []
    status = []
    vartype = []
    for idx,(i,j) in enumerate(zip(seq1,seq2)):
        sites.append(idx)
        summary.append(i+str(idx+1)+j)
        if i==j:
            status.append('MATCH')
            if i=="-" and j=="-":
                vartype.append('SHARED_DELELTION')
            else:
                vartype.append('NA')
        else:
            if i == "-" or j == "-":
                vartype.append('INDEL')
            else:
                vartype.append('SNP')
            status.append('DIFFERENCE')
    summary_df = pd.concat([pd.Series(sites,name='sites'),pd.Series(summary,name='summary'),pd.Series(status,name='status'),pd.Series(vartype,name='vartype')], axis=1)
    full_differences = summary_df[summary_df.status.isin(['DIFFERENCE'])]
    if len(full_differences.index) > 0:
        reduced_differences = full_differences[~full_differences.sites.isin(removals)]
        num_reduced_differences = len(reduced_differences.index)
        if gap_identity == 'compressed': 
            compressed_differences = reduced_differences
            consecutive_indels = find_consecutive_gaps(reduced_differences)
            for key in consecutive_indels:
                if len(consecutive_indels[key]) != 0:
                    for i in range(len(consecutive_indels[key])):
                        num_reduced_differences = (num_reduced_differences - consecutive_indels[key][i][2]) + 1
                        beginning_site = consecutive_indels[key][i][0] 
                        end_site = consecutive_indels[key][i][1]
                        destroy = list(range(beginning_site,end_site))
                        compressed_differences = compressed_differences[~compressed_differences.isin(destroy).iloc[:,0]]
                        replacement_summary = compressed_differences.at[beginning_site-1, 'summary'].replace('-','del'+str(consecutive_indels[key][i][2]))
                        compressed_differences.at[beginning_site-1, 'summary'] = replacement_summary
            final_differences = '; '.join([x for x in compressed_differences['summary']])
        else:
            final_differences = '; '.join([x for x in reduced_differences['summary']])
    else:
        num_reduced_differences = 0
        final_differences = ''
    full_matches = summary_df[summary_df.status.isin(['MATCH'])]
    reduced_matches = full_matches[~full_matches.sites.isin(removals)]
    num_reduced_matches = len(reduced_matches.index)
    reduced_total_length = num_reduced_matches+num_reduced_differences
    PID = "{:.2f}".format((num_reduced_matches/reduced_total_length)*100)
    if reduced_total_length > 0 and PID == '100.00':
        PID = (num_reduced_matches/reduced_total_length)*100
    if args.translate:
        if num_reduced_differences != 0:
            difference_sites = [x for x in reduced_differences.loc[reduced_differences['status'] == 'DIFFERENCE']['sites'].to_list()]
            aa_differences = get_aa_differences(seq1,seq2,annotations,difference_sites)
            my_order = [int(re.sub(pattern=':.*',repl='',string=x)) for x in aa_differences]
            my_order = [i[0] for i in sorted(enumerate(my_order), key=lambda x:x[1])]
            aa_differences = '; '.join([aa_differences[x] for x in my_order])
        else:
            aa_differences = ''
        calc_dict = {'Seq1':[name1],
                     'Seq2':[name2],
                     'Full_Alignment_Length':[full_length],
                     'Trimmed_Alignment_Length':[reduced_total_length],
                     'Num_Differences':[num_reduced_differences],
                     'PID':[PID],
                     'Differences':[final_differences],
                     'AA_Differences':[aa_differences]}
    else:
        calc_dict = {'Seq1':[name1],
                    'Seq2':[name2],
                    'Full_Alignment_Length':[full_length],
                    'Trimmed_Alignment_Length':[reduced_total_length],
                    'Num_Differences':[num_reduced_differences],
                    'PID':[PID],
                    'Differences':[final_differences]}      
    return(calc_dict)
#%%
def find_removals(seq1,seq2,ambiguity_codes,ignore_gaps,ignore_ambiguous,trim_ends,ignore_indels):
    '''
    Define the sites to be removed for PID calculations
    '''
    gaps = ['N','-']
    if ignore_gaps == True and ignore_ambiguous == True:
        bad_sites = ambiguity_codes+gaps
    elif ignore_gaps == True and ignore_ambiguous == False:
        bad_sites = gaps
    elif ignore_gaps == False and ignore_ambiguous == True:
        bad_sites = ambiguity_codes
    trim_sites = []
    if trim_ends == True:
        idx1 = indices(seq1)
        idx2 = indices(seq2)
        if idx1[0] >= idx2[0]:
            start = idx1[0]
        else:
            start = idx2[0]
        if idx1[1] >= idx2[1]:
            end = idx2[1]+1
        else:
            end = idx1[1]+1
        trim_sites = list(range(0,start+1))+list(range(end,len(seq1)+1))
    removals = []
    if ignore_gaps == True or ignore_ambiguous == True:
        a = [i for i,x in enumerate([y for y in seq1]) if x in bad_sites]
        b = [i for i,x in enumerate([y for y in seq2]) if x in bad_sites]
        removals = list(set().union(a,b))
    removals = list(set().union(removals,trim_sites))
    indel_sites = []
    for idx,(i,j) in enumerate(zip(seq1,seq2)):
        if (i not in ['-','N'] and j=='-') or (j not in ['-','N'] and i=='-'):
            if idx not in trim_sites:
                indel_sites.append(idx)
    if ignore_indels == False:
        removals = [x for x in removals if x not in indel_sites]
    removals.sort()
    return(removals)

#%%
def compare_seqs_parallel_old(pair,MSA,accession_list,gap_identity,ignore_indels,ignore_gaps,ignore_ambiguous,trim_ends,seqtype, annotations= ''):
#    print('pair: ',pair)
#    print('msa: ',type(MSA))
#    print('gap identity: ',gap_identity)
#    print('ignore indels: ',ignore_indels)
#    print('ignore gaps: ',ignore_gaps)
#    print('ignore ambig: ',ignore_ambiguous)
#    print('trim ends: ',trim_ends)
#    print('seqtype: ',seqtype)
#    print('annotations: ', annotations)
    ppair = pair.split(sep='\t')
    ppair = [x.split()[0] for x in ppair]
    sequences = []
    for record in MSA:
        if record.id in ppair:
            sequences.append(record)
    if seqtype != 'nucleotide':
        ambiguity_codes = []
    else:
        ambiguity_codes = ['M','R','W','S','Y','K','V','H','D','B']
    seq1 = sequences[0].seq
    name1 = sequences[0].id
    seq2 = sequences[1].seq
    name2 = sequences[1].id
    removals = find_removals(seq1,seq2,ambiguity_codes,ignore_gaps,ignore_ambiguous,trim_ends,ignore_indels)
    result = calc_differences(name1,seq1,name2,seq2,removals, annotations, gap_identity)
    return pd.DataFrame.from_dict(result, orient='columns')

#%%
def compare_seqs_parallel(pair,MSA,accession_list,gap_identity,ignore_indels,ignore_gaps,ignore_ambiguous,trim_ends,seqtype, annotations= ''):
    ppair = pair.split(sep='\t')
    ppair = [x.split()[0] for x in ppair]
    i1 = accession_list.index(ppair[0])
    i2 = accession_list.index(ppair[1])
    if seqtype != 'nucleotide':
        ambiguity_codes = []
    else:
        ambiguity_codes = ['M','R','W','S','Y','K','V','H','D','B']
    seq1 = MSA[i1].seq
    name1 = MSA[i1].id
    seq2 = MSA[i2].seq
    name2 = MSA[i2].id
    removals = find_removals(seq1,seq2,ambiguity_codes,ignore_gaps,ignore_ambiguous,trim_ends,ignore_indels)
    result = calc_differences(name1,seq1,name2,seq2,removals, annotations, gap_identity)
    return pd.DataFrame.from_dict(result, orient='columns')

#%%
def parse_commands():
    '''
    Parse command-line arguments.
    '''
    global args
    #"Parse the input arguments, use -h for help"
    epilog = '''
    Compares pairs of sequences described within a tab-delimited text file to
    determine the pairwise identity between each pair, and report the differences.
    Optional flags determine how differences are calculated: indels, gaps, and
    ambiguous sites can be ignored, and runs of 'N' or '-' from the ends of sequences
    can also be ignored in calculations. Each of these flags can be used in combination.
    Internal gaps are defined as 'N' if indels are included in calculations (default), but are re-defined
    as '-' or 'N' if indels are ignored (--ignore_indels).
    '''
    parser=argparse.ArgumentParser(description="Sequence Pairwise Comparison Script (version "+version+')', formatter_class=argparse.ArgumentDefaultsHelpFormatter, epilog=epilog)
    #Read inputs
    parser.add_argument('-a','--annotations', required=False, default=False, help='Modified bed file with annotations for CDS regions (see example file)')
    parser.add_argument('-c','--comparisons', required=True, help='Tab-delimited text file (no header) with seqs to compare: col1 = seq1, col2 = seq2', metavar='<.tsv file>')
    parser.add_argument('-d','--data_type', required=False, default='nucleotide', help="Type of data in alignment", metavar="<nucleotide / protein>")
    parser.add_argument('-f','--fasta', required=True, help="Fasta file containing aligned sequences")
    parser.add_argument('-ii','--ignore_indels', required=False, default=False, action='store_true', help="Ignore indels when calculating pairwise identity")
    parser.add_argument('-ig','--ignore_gaps', required=False, default=False, action='store_true', help="Ignore gaps when calculating pairwise identity")
    parser.add_argument('-ia','--ignore_ambiguous', required=False, default=False, action='store_true', help="Ignore ambiguous sites when calculating pairwise identity")
    parser.add_argument('-te','--trim_ends', required=False,  action='store_true', default=False, help='Trim runs of gaps ("N" or "-") at the beginning and end of alignments')
    parser.add_argument('-tr','--translate', required=False,  action='store_true', default=False, help='EXPERIMENTAL: Translate CDS regions from annotation file to calculate amino acid consequences')
    parser.add_argument('--gap_identity', required=False, default='compressed', help='Method to score gaps: "compressed" (default) or "blast". See: https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity')
    parser.add_argument('-o','--outfile_stem', required=False, default='results', help='Stem name for the results outfiles. Default = "results" --> "results.all.csv"; "results.matrix.csv"')
    args=parser.parse_args()
    return()

#%%
def main():
    printc("\nSequence Pairwise Comparison Script (version "+version+')', 'green')
    parse_commands()
    if args.ignore_indels:
        ignore_indels = True
    else:
        ignore_indels = False
    if args.ignore_gaps:
        ignore_gaps = True
    else:
        ignore_gaps = False
    if args.ignore_ambiguous:
        ignore_ambiguous = True
    else:
        ignore_ambiguous = False
    if args.trim_ends:
        trim_ends = True
    else:
        trim_ends = False
    if args.gap_identity:
        if args.gap_identity not in ['compressed','blast']:
            print('\nError: gap identity must be one of "compressed" or "blast"\n')
            sys.exit(-1)      
        gap_identity = args.gap_identity
    else:
        gap_identity = "compressed"

    seqtype = args.data_type.lower()
    MSA = AlignIO.read(args.fasta,'fasta')
    
    accession_list = [record.id for record in MSA]
    
    if args.translate and not args.annotations:
        print('\nError: cannot calculate the amino acid consequences of differences without an annotation file\n')
        sys.exit(-1)
    
    annotations = ''
    if args.translate:
        printc("\nWarning: amino acid translation results are experimental.\nI don't recommend you use this.\n", 'yellow')
        annotations = pd.read_csv(args.annotations, header=None, sep="\t")
        insertions = find_insertions(MSA, annotations)
        fix_bed_insertions(annotations, insertions)

    with open(args.comparisons,'r') as file:
        comparisons = file.read().splitlines()
    
    check_input(comparisons, MSA)
    printc('Calculating pairwise similarities','blue')
    print('Number of sequences: {0}'.format(len(MSA)))
    print('Number of comparisons: {0}'.format(len(comparisons)))
    if tqdm_installed == True:
        with mp.Pool(mp.cpu_count()) as pool:
            results = pool.starmap(compare_seqs_parallel,
                                   tqdm.tqdm([(x, MSA, accession_list, gap_identity,ignore_indels,ignore_gaps,ignore_ambiguous,trim_ends,seqtype, annotations) for x in comparisons], total=len(comparisons)))#,
                                   #chunksize=1)
    else:
        print("If you install 'tqdm', you can get a nice progress bar while the program runs")
        print("Try running 'python3 -m pip install tqdm' before you run this program again")
        with mp.Pool(mp.cpu_count()) as pool:
            results = pool.starmap(compare_seqs_parallel,
                                   [(x,MSA,accession_list,gap_identity,ignore_indels,ignore_gaps,ignore_ambiguous,trim_ends,seqtype,annotations) for x in comparisons])#,
                                   #chunksize=1)
    final = pd.concat(results)
    final.to_csv(args.outfile_stem+".all.csv", index = False)
    results_matrix = convert_to_matrix(final)
    results_matrix.to_csv(args.outfile_stem+".matrix.csv")
    PID = pd.Series([float(x) for x in final['PID']],name='PID')
    mean_pid = PID.mean()
    median_pid = PID.median()
    sd_pid = PID.std()
    print("Mean PID: {0}".format(mean_pid))
    print("Median PID: {0}".format(median_pid))
    print("Stdev PID: {0}".format(sd_pid))
    printc('Analysis complete: check result in {0} and {1}\n'.format(args.outfile_stem+".all.csv", args.outfile_stem+".matrix.csv"),'blue')
#%%
if __name__ == '__main__':
    main()
