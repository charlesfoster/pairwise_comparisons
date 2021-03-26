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
pairwise_comparisons.py [-h] -c <.txt file> [-d <nucleotide / protein>] -f FASTA [-ii] [-ig] [-ia] [-te] [-o OUTFILE]
```

optional arguments:
  -h, --help            show help message and exit
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
  -o OUTFILE, --outfile OUTFILE
                        Name of the outfile to store results (default:
                        results.csv)

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
* At some stage I will update the program to have an option to calculate the effects of nucleotide differences at the amino acid level (i.e., whether differences are synonymous or non-synonymous). This approach will be via providing gene coordinates, or via implementing a blast annotation step. Haven't decided yet.

## Citation
If you use this software and find it useful, I'd appreciate some kind of attribution, e.g.:

Foster, C.S.P, Pairwise Comparisons, (2021), GitHub repository, https://github.com/charlesfoster/pairwise_comparisons

I might add a zenodo citation. If I build on the program at some stage I'll also consider publishing it somewhere, and I'll add a citation here.
