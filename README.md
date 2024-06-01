# Alignment-CSE185-SP24-

## Overview
This project implements a method to perform local alignment of DNA sequences, following the BWT (Burrows-Wheeler Transform) and exact matching filter approach similar to `bwa mem`. The tool reads sequences from a FASTA file, performs BWT, applies exact matching, and then aligns the sequences, outputting the results, including alignment scores and lengths.

## Features
- BWT transformation and exact matching filter
- Local alignment
- Command-line utility for easy use
- Outputs alignment results to a specified file
- Optionally plots the distribution of alignment lengths

## Table of Contents

- [Contributors](#contributors)
- [So...what is an Alignment](#sowhat-is-an-alignment)
  - [Application of Alignments](#application-of-alignments)
  - [What is BWT and exact matching filter?](#what-is-bwt-and-exact-matching-filter)
- [How do I run it?](#how-do-i-run-it)
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Usage](#usage)
  - [File Output Format](#file-output-format)
- [Why should I run it? (Benchmarking Analysis)](#why-should-i-run-it-benchmarking-analysis)

## Contributors

- Gary Lin (A16915179)
- Lorenzo Olmo Marchal (A17013640)
- Nabeeha Rashid (A16851960)

This project is intended to be our final project in CSE 185: Advanced Bioinformatics Laboratory for the Spring 2024 quarter, to be presented to our peers and Professor Melissa Gyrmrek.

## So...what is an Alignment?
An alignment is a process of arranging DNA, RNA, or protein sequences to identify regions of similarity. These similarities can be due to functional, structural, or evolutionary relationships. Alignments help in comparing sequences to find conserved regions, which are important for understanding biological functions and evolutionary histories.

### Application of Alignments
Alignments are used in many fields. Here are some examples of many:
- **Comparative Genomics:** Identifying conserved sequences across different species
- **Phylogenetics:** Constructing evolutionary trees
- **Evolutionary Studies:** Investigating the evolutionary relationships between organisms
- **Medical Research:** Identifying mutations associated with diseases

### What is BWT and exact matching filter?
TODO

# How do I Run it?

## Requirements
TODO

## Installation
Clone the repository and navigate to the project directory:
```sh
git clone https://github.com/ADDOURLINKHERe
cd ADD
```

Ensure you have Python 3.x and the necessary packages installed:
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

