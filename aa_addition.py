#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 11:34:08 2021

@author: cfos
"""
annotation_bed = "/home/cfos/Programs/pairwise_comparisons/SARS-CoV_2_annotations.bed"
annotations = pd.read_csv(annotation_bed, header=None, sep="\t")
fasta = '/home/cfos/Documents/Collaboration/Elsa/all_samples.fasta'

MSA = AlignIO.read(fasta,'fasta')


ppair = ["24_S44","EC7_S19"]
sequences = []
for record in MSA:
    if record.id in ppair:
        sequences.append(record)
seq1 = sequences[0].seq
name1 = sequences[0].id
seq2 = sequences[1].seq
name2 = sequences[1].id

ignore_gaps = False
ignore_ambiguous = False
trim_end = False
ignore_indels = False

difference_sites = [x for x in reduced_differences.loc[reduced_differences['status'] == 'DIFFERENCE']['sites'].to_list()]

import pandas as pd

def parse_annotation_bed(annotation_bed):
    annotations = pd.read_csv(annotation_bed, header=None, sep="\t")
    anno_list = []
    for i in range(len(annotations.index)):
        gene_dict = {'name':annotations.iloc[i][3],
                     'start':int(annotations.iloc[i][1]),
                     'end':int(annotations.iloc[i][2])-1}
        anno_list.append(gene_dict)
    return anno_list

aa_results = []
for site in difference_sites:
    for i in range(len(annotations.index)):
        gene_range = range(int(annotations.iloc[i][1]),int(annotations.iloc[i][2]))
        if site in gene_range:
            reduced_df = annotations.iloc[i]
            start = reduced_df[1]
            end = int(reduced_df[2])
            gene_name = reduced_df.iloc[3]
            if reduced_df[4] == True:
                slip_site = int(reduced_df[5])
                aa_seq1a = seq1[start:slip_site+1]
                aa_seq1b = seq1[slip_site:end]
                aa_seq1 = aa_seq1a+aa_seq1b
                aa_seq2a = seq2[start:slip_site+1]
                aa_seq2b = seq2[slip_site:end]
                aa_seq2 = aa_seq2a+aa_seq2b
            else:
                aa_seq1 = seq1[start:end]
                aa_seq2 = seq2[start:end]
            aa_seq1 = aa_seq1.translate()
            aa_seq2 = aa_seq2.translate()
            for idx,(i,j) in enumerate(zip(aa_seq1,aa_seq2)):
                if i!=j:
                    result = gene_name+':'+str(i)+str(idx)+str(j)
                    aa_results.append(result)
            
aa_results = []
for i in range(len(annotations.index)):
    gene_range = range(int(annotations.iloc[i][1]),int(annotations.iloc[i][2]))
    reduced_df = annotations.iloc[i]
    start = reduced_df[1]
    end = int(reduced_df[2])
    gene_name = reduced_df.iloc[3]
    for site in difference_sites:
        if site in gene_range:
            mod = list(gene_range).index(site) % 3
            if mod == 1:
                codon_pos = 1
                slice_range = range(site,site+3)
            elif mod == 2:
                codon_pos = 2
                slice_range = range(site-1,site+2)           
            else:
                codon_pos = 3
                slice_range = range(site-2,site+1)           
            if reduced_df[4] == True:
                slip_site = int(reduced_df[5])
                aa_seq1a = seq1[start:slip_site+1]
                aa_seq1b = seq1[slip_site:end]
                aa_seq1 = aa_seq1a+aa_seq1b
                aa_seq2a = seq2[start:slip_site+1]
                aa_seq2b = seq2[slip_site:end]
                aa_seq2 = aa_seq2a+aa_seq2b
            else:
                aa_seq1 = seq1[start:end]
                aa_seq2 = seq2[start:end]            
            
            aa_seq1 = ''.join([x.replace('-','N') for x in list(aa_seq1)])
            aa_seq1 = SeqRecord(Seq(aa_seq1))
            aa_seq1 = aa_seq1.translate()
            aa_seq2 = ''.join([x.replace('-','N') for x in list(aa_seq2)])
            aa_seq2 = SeqRecord(Seq(aa_seq2))
            aa_seq2 = aa_seq2.translate()

            for idx,(i,j) in enumerate(zip(aa_seq1,aa_seq2)):
                if i!=j:
                    result = gene_name+':'+str(i)+str(idx)+str(j)
                    aa_results.append(result)


######
            ##
aa_results = []
for i in range(len(annotations.index)):
    gene_range = range(int(annotations.iloc[i][1]),int(annotations.iloc[i][2]))
    reduced_df = annotations.iloc[i]
    start = reduced_df[1]
    end = int(reduced_df[2])
    gene_name = reduced_df.iloc[3]
    for site in difference_sites:
        if site in gene_range:
            if reduced_df[4] == True and int(reduced_df[5]) < site:
                site = site +1
                start = start-1
            mod = list(gene_range).index(site+1) % 3
            if mod == 1:
                codon_pos = 1
                slice_range = range(site,site+4)
            elif mod == 2:
                codon_pos = 2
                slice_range = range(site-1,site+3)
            else:
                codon_pos = 3
                slice_range = range(site-2,site+2)
            codon_num = int((slice_range[-2] - (start-1))/3)
            if reduced_df[4] == True:
                slip_site = int(reduced_df[5])
                aa_seq1a = seq1[0:slip_site+1]
                aa_seq1b = seq1[slip_site:]
                aa_seq1 = aa_seq1a+aa_seq1b
                aa_seq2a = seq2[0:slip_site+1]
                aa_seq2b = seq2[slip_site:]
                aa_seq2 = aa_seq2a+aa_seq2b
            aa_seq1 = ''.join([x.replace('-','N') for x in list(aa_seq1)])
            aa_seq1 = SeqRecord(Seq(aa_seq1))
            aa_seq2 = ''.join([x.replace('-','N') for x in list(aa_seq2)])
            aa_seq2 = SeqRecord(Seq(aa_seq2))
            codon1 = aa_seq1[slice_range[0]:slice_range[-1]]
            codon2 = aa_seq2[slice_range[0]:slice_range[-1]]
            codon1_AA = codon1.translate()
            codon2_AA = codon2.translate()
            if codon1_AA.seq == codon2_AA.seq:
                result = gene_name+':'+str(codon1_AA.seq)+str(codon_num)+str(codon2_AA.seq)+":s"
            elif codon1_AA.seq != codon2_AA.seq:
                result = gene_name+':'+str(codon1_AA.seq)+str(codon_num)+str(codon2_AA.seq)+":n"
            aa_results.append(result)


#####
refseq = []
for record in MSA:
    if record.id == annotations.iloc[0][0]:
        print('hooray')
        refseq.append(record)

refseq = [x for x in refseq[0]]
refseq.index('-')

insertions = [i for i, ltr in enumerate(refseq[0].seq) if ltr == '-' and i < 29902] # need to change to args.reflength
annotations = pd.read_csv(annotation_bed, header=None, sep="\t")

def fix_bed_insertions(annotations, insertions):
    for i in annotations.index:
        print('Checking {0} which spans {1} to {2}'.format(annotations.iloc[i][3],annotations.iloc[i][1],annotations.iloc[i][2]))
        for site in insertions:
            print('Checking site: {}'.format(site))
            if site < annotations.iloc[i][1] and site < (annotations.iloc[i][2]-1):
                print('Changing coords of {}'.format(annotations.iloc[i][3]))
                annotations.at[i,1] = annotations.at[i,1] + 1
                annotations.at[i,2] = annotations.at[i,2] + 1
            elif site > annotations.iloc[i][1]  and site < (annotations.iloc[i][2]-1):
                print('Changing coords of {}'.format(annotations.iloc[i][3]))
                annotations.at[i,2] = annotations.at[i,2] + 1
            else:
                print('No need to update {}'.format(annotations.iloc[i][3]))

def check_is_consecutive(l):
    maximum = max(l)
    if sum(l) == maximum * (maximum+1) /2 : 
        return True
    return False

def ranges(nums):
    nums = sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    return list(zip(edges, edges))
