# Pairwise Comparisons
A program written in python3 to allow pairwise comparisons of DNA/protein sequences.

## Installation
Firstly, clone this repository:

```
git clone https://github.com/charlesfoster/pairwise_comparisons.git

cd pairwise_comparisons
```

Make the python script executable:

```
chmod +x ./pairwise_comparisons.py
```

Optional: install the 'tqdm' python3 package to get a nice progress bar while the program runs:

```
python3 -m pip install tqdm
```

If you choose not to install 'tqdm' the program will still run just as well, but you'll be missing out on a visual experience.

## What does the program do?
Let's say you have sequences of molecular data (DNA or protein sequences) and you want to compare them to determine the pairwise identity between the sequences, and report what the differences are. That's when this program comes in handy. The program can work with any number of sequences within a given alignment, as long as you tell it which pairs of sequences you want to compare. All desired comparisons are run in parallel: with my 10-core processor (Intel(R) Core(TM) i9-10900X CPU @ 3.70GHz), I can run 253 pairwise comparisons in just over 2 seconds (111.78 comparisons per second).

By default, the program considers all sites in the alignment. However, there are a number of options that allow you to successively ignore certain sites in the alignment. You can choose to (a) trim runs of 'N' or '-' from the end of the alignment (common in consensus genomes, and/or aligned sequences of different lengths), (b) ignore indel sites, (a) ignore gaps in the alignment, and/or (d) ignore sites with IUPAC ambiguity codes.

Each of these options can be used in combination. Internal gaps are defined as 'N' if indels are included in
calculations (default), but are re-defined as '-' or 'N' if indels are ignored.

## What do you need to run the program?
* <alignment>.fasta: a multiple sequence alignment in fasta format (the sequences must be aligned beforehand)
* <comparisons>.tsv: a tab-delimited text file names of sequences to be compared (the names must exactly match the sequence headers in the fasta file)

## How do you run the program?
General usage:

```
Sequence Pairwise Comparison Script (version 1.0.0)
usage: pairwise_comparisons.py [-h] [-a ANNOTATIONS] -c <.tsv file>
                               [-d <nucleotide / protein>] -f FASTA [-ii]
                               [-ig] [-ia] [-te] [-tr] [-o OUTFILE]

Sequence Pairwise Comparison Script (version 1.0.0)

optional arguments:
  -h, --help            show this help message and exit
  -a ANNOTATIONS, --annotations ANNOTATIONS
                        Modified bed file with annotations for CDS regions
                        (see example file) (default: False)
  -c <.tsv file>, --comparisons <.tsv file>
                        Tab-delimited text file (no header) with seqs to
                        compare: col1 = seq1, col2 = seq2 (default: None)
  -d <nucleotide / protein>, --data_type <nucleotide / protein>
                        Type of data in alignment (default: nucleotide)
  -f FASTA, --fasta FASTA
                        Fasta file containing aligned sequences (default:
                        None)
  -ii, --ignore_indels  Ignore indels when calculating pairwise identity
                        (default: False)
  -ig, --ignore_gaps    Ignore gaps when calculating pairwise identity
                        (default: False)
  -ia, --ignore_ambiguous
                        Ignore ambiguous sites when calculating pairwise
                        identity (default: False)
  -te, --trim_ends      Trim runs of gaps ("N" or "-") at the beginning and
                        end of alignments (default: False)
  -tr, --translate      Translate CDS regions from annotation file to
                        calculate amino acid consequences (default: False)
  -o OUTFILE, --outfile OUTFILE
                        Name of the outfile to store results (default:
                        results.csv)

Compares pairs of sequences described within a tab-delimited text file to
determine the pairwise identity between each pair, and report the differences.
Optional flags determine how differences are calculated: indels, gaps, and
ambiguous sites can be ignored, and runs of 'N' or '-' from the ends of
sequences can also be ignored in calculations. Each of these flags can be used
in combination. Internal gaps are defined as 'N' if indels are included in
calculations (default), but are re-defined as '-' or 'N' if indels are ignored
(--ignore_indels).
```

## What do you get?
The output file (in .csv format) contains statistics about the sequences that were compared. The columns, in order, contain:

* Name of sequence1
* Name of sequence2
* Full alignment length
* Alignment length after trimming based on command line options
* Number of differences between the two sequences
* Pairwise identity percentage (to two d.p., unless rounding would falsely bring the identity to 100%, in which case you'll get the results to many d.p.)
* The differences between the sequences (in the form of <Base_in_Seq1>Position<Base_in_Seq2>, e.g. A111T)

At the end of the analysis, the mean, median, and stdev pairwise identity (across all comparisons) will be printed to stdout in your terminal.

## Future improvements
* I've added an option to translate sequences and get the amino acid consequences. It works, but could be improved. I'll work on it when I have time.

## Citation
If you use this software and find it useful, I'd appreciate some kind of attribution, e.g.:

Foster, C.S.P, Pairwise Comparisons, (2021), GitHub repository, https://github.com/charlesfoster/pairwise_comparisons

I might add a zenodo citation. If I build on the program at some stage I'll also consider publishing it somewhere, and I'll add a citation here.
