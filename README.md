# BWT Sequence Aligner-CSE185-SP24-

## Overview
The BWT Sequence Aligner is a tool that implements a sequence aligner using the Burrows-Wheeler Transform (BWT) and FM-index for efficient pattern matching and sequence alignment. This approach is similar to `bwa mem`, a tool we used in lab. The tool reads sequencing reads from a FASTQ file, aligns them to a reference genome, and outputs the alignment results in SAM (Sequence Alignment/Map) format.

## Features
- BWT transformation
- FM index
- Local alignment
- Command-line utility for easy use
- Outputs results to a specified SAM file

## Table of Contents

- [Contributors](#contributors)
- [So...what is an Alignment](#sowhat-is-an-alignment)
  - [Application of Alignments](#application-of-alignments)
  - [What is BWT](#what-is-bwt)
- [How do I run it?](#how-do-i-run-it)
  - [Installation](#installation)
  - [Usage](#usage)
  - [File Output Format](#file-output-format)
- [Why should I run it? (Benchmarking Analysis)](#why-should-i-run-it-benchmarking-analysis)

## Contributors

- Gary Lin (A16915179)
- Lorenzo Olmo Marchal (A17013640)
- Nabeeha Rashid (A16851960)

This project is intended to be our final project in CSE 185: Advanced Bioinformatics Laboratory for the Spring 2024 quarter, to be presented to our peers and Professor Melissa Gyrmrek.

## File Descriptions
Here is a brief description of the files:

#### Code Files
- The `ReferenceIndexing.py` file implements the indexing of a reference genome sequence using the Burrows-Wheeler Transform (BWT) and the FM-index. This indexing process efficiently preprocesses the reference genome sequence to enable fast pattern matching and retrieval of substring occurrences.
- The `seedAndExtend.py` file implements a local alignment approach. It first identifies potential seed matches between the query sequence and the reference genome using exact matching with a k-mer approach. Then, it extends these seeds locally to find the best alignment between the sequences. This local alignment strategy allows for more accurate alignments, particularly in regions where sequences may differ due to insertions, deletions, or substitutions.

#### Data
- 

## So...what is an Alignment?
An alignment is a process of arranging DNA, RNA, or protein sequences to identify regions of similarity. These similarities can be due to functional, structural, or evolutionary relationships. Alignments help in comparing sequences to find conserved regions, which are important for understanding biological functions and evolutionary histories.

### Application of Alignments
Alignments are used in many fields. Here are some examples of many:
- **Comparative Genomics:** Identifying conserved sequences across different species
- **Phylogenetics:** Constructing evolutionary trees
- **Evolutionary Studies:** Investigating the evolutionary relationships between organisms
- **Medical Research:** Identifying mutations associated with diseases

### What is a BWT


# How do I Run it?

## Installation
Clone the repository and navigate to the project directory:
```sh
git clone https://github.com/ADDOURLINKHERe
cd ADD
```

ADD REQUIREMENTS AND VERSIONS HEREEEEEEEEE:
```sh
pip install numpy
pip install matplotlib
```

## Usage
To use the tool, run the script with the required arguments--an example is shown below:
```sh
python example.py sequence_name.fasta -m 1 -s -1 -d -2 -a -o output.txt -p
```

### Arguments
- `seq_files`: Path to the FASTA file containing sequences to be aligned.
- `-m`, `--match`: Match reward (positive integer).
- `-s`, `--mismatch`: Mismatch penalty (negative integer).
- `-d`, `--indel`: Indel penalty (negative integer).
- `-a`, `--alignment`: Include aligned sequences in the output.
- `-o`, `--output`: Output file path (default: `./Alignment.txt`).
- `-p`, `--plot`: Plot the distribution of alignment lengths.

### Test Dataset Command
Here is the command to run our test dataset. Reminder: If you do not include `-a`, then the output file will only show alignment scores and lengths, but will not show alignments!

An example FASTA file from Lab 1, the public 'Genomes' folder (`hg19.fa`), is included in the `data` directory:
```sh
python run_alignmentfilenameputhere.py data/hg19.fa -m 1 -s -1 -d -2 -a -o output.txt -p
```

## File Output Format
The output file will contain the aligned sequences with their scores and lengths. Example format:
```plaintext
>seq 1|score: 10|length: 50
AGCTGA...
A-C-GA...

>seq 2|score: 8|length: 45
CGATGC...
C-ATG-...
```

## Why should I run it? (Benchmarking Analysis)
TODO

