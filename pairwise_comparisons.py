#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 14:46:03 2021

@author: Dr Charles Foster, Virology Research Lab
"""
import pandas as pd
from Bio import AlignIO
import argparse
import multiprocessing as mp
try:
    import tqdm
    tqdm_installed = True
except:
    tqdm_installed = False
version = '1.0.0'
#%%
def printc(thing, level):
    '''
    Print in colour :)
    '''
    cols = {'green':'\033[1;32m', 'blue':'\033[96m'}
    col = cols[level]
    print(f"{col}{thing}\033[0m")
    return()

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
def calc_differences(name1,seq1,name2,seq2,removals):
    ''''
    Calculate total pairwise differences. These are later pruned based on trimming conditions.
    '''
    full_length = len(seq1) #aligned, so both seqs same length
    sites = []
    summary = []
    status = []
    for idx,(i,j) in enumerate(zip(seq1,seq2)):
        sites.append(idx)
        summary.append(i+str(idx+1)+j)
        if i==j:
            status.append('MATCH')
        else: 
            status.append('DIFFERENCE')
    summary_df = pd.concat([pd.Series(sites,name='sites'),pd.Series(summary,name='summary'),pd.Series(status,name='status')], axis=1)
    full_differences = summary_df[summary_df.status.isin(['DIFFERENCE'])]
    if len(full_differences.index) > 0:
        reduced_differences = full_differences[~full_differences.sites.isin(removals)]
        num_reduced_differences = len(reduced_differences.index)
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
def compare_seqs_parallel(pair,MSA,ignore_indels=False,ignore_gaps=False,ignore_ambiguous=False,trim_ends=True,seqtype='nucleotide'):
    ppair = pair.split(sep='\t')
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
    result = calc_differences(name1,seq1,name2,seq2,removals)
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
    parser.add_argument('-c','--comparisons', required=True, help='Tab-delimited text file (no header) with seqs to compare: col1 = seq1, col2 = seq2', metavar='<.tsv file>')
    parser.add_argument('-d','--data_type', required=False, default='nucleotide', help="Type of data in alignment", metavar="<nucleotide / protein>")
    parser.add_argument('-f','--fasta', required=True, help="Fasta file containing aligned sequences")
    parser.add_argument('-ii','--ignore_indels', required=False, default=False, action='store_true', help="Ignore indels when calculating pairwise identity")
    parser.add_argument('-ig','--ignore_gaps', required=False, default=False, action='store_true', help="Ignore gaps when calculating pairwise identity")
    parser.add_argument('-ia','--ignore_ambiguous', required=False, default=False, action='store_true', help="Ignore ambiguous sites when calculating pairwise identity")
    parser.add_argument('-te','--trim_ends', required=False,  action='store_true', default=False, help='Trim runs of gaps ("N" or "-") at the beginning and end of alignments')
    parser.add_argument('-o','--outfile', required=False, default='results.csv', help='Name of the outfile to store results')
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
    seqtype = args.data_type.lower()
    MSA = AlignIO.read(args.fasta,'fasta')
    with open(args.comparisons,'r') as file:
        comparisons = file.read().splitlines()
    printc('Calculating pairwise similarities','blue')
    print('Number of sequences: {0}'.format(len(MSA)))
    print('Number of comparisons: {0}'.format(len(comparisons)))
    if tqdm_installed == True:
        with mp.Pool(mp.cpu_count()) as pool:
            results = pool.starmap(compare_seqs_parallel, 
                                   tqdm.tqdm([(x, MSA,ignore_indels,ignore_gaps,ignore_ambiguous,trim_ends,seqtype) for x in comparisons], total=len(comparisons)),
                                   chunksize=1)
    else:
        print("If you install 'tqdm', you can get a nice progress bar while the program runs")
        print("Try running 'python3 -m pip install tqdm' before you run this program again")
        with mp.Pool(mp.cpu_count()) as pool:
            results = pool.starmap(compare_seqs_parallel, 
                                   [(x, MSA,ignore_indels,ignore_gaps,ignore_ambiguous,trim_ends,seqtype) for x in comparisons],
                                   chunksize=1)
    final = pd.concat(results)
    mean_pid = final['PID'].mean()
    median_pid = final['PID'].median()
    sd_pid = final['PID'].std()
    print("Mean PID: {0}".format(mean_pid))
    print("Median PID: {0}".format(median_pid))
    print("Stdev PID: {0}".format(sd_pid))
    final.to_csv(args.outfile, index = False)
    printc('Analysis complete: check result in {0}\n'.format(args.outfile),'blue')
#%%
if __name__ == '__main__':
    main()
