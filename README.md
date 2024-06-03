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
- [File Descriptions](#file-descriptions)
- [So...what is an Alignment](#sowhat-is-an-alignment)
  - [Application of Alignments](#application-of-alignments)
  - [What is BWT and FM Index](#what-is-bwt-and-fm-index)
- [How do I run it?](#how-do-i-run-it)
  - [Installation](#installation)
  - [Usage](#usage)
  - [File Output Format](#file-output-format)

## Contributors

- Gary Lin (A16915179)
- Lorenzo Olmo Marchal (A17013640)
- Nabeeha Rashid (A16851960)

This project is intended to be our final project in CSE 185: Advanced Bioinformatics Laboratory for the Spring 2024 quarter, to be presented to our peers and Professor Melissa Gyrmrek.

## File Descriptions
Here is a brief description of the files:

#### Code Files
- The `ReferenceIndexing.py` file implements the indexing of a reference genome sequence using the Burrows-Wheeler Transform (BWT) and the FM-index. 
- The `seedAndExtend.py` file implements a local alignment approach. It first identifies potential seed matches between the query sequence and the reference genome using exact matching with a k-mer approach. Then, it extends these seeds locally to find the best alignment between the sequences. 

#### Data
- The `ERR10021327.fastq` file is a publically available query file containing the reads that are to be aligned. This reads file is from [ENA project PRJEB37886](https://www.ebi.ac.uk/ena/browser/view/PRJEB37886). 
- The `sequence.fasta` file is a publically available reference genome file from [NCBI accession ERR10021327](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2).

## So...what is an Alignment?
An alignment is a process of arranging DNA, RNA, or protein sequences to identify regions of similarity. These similarities can be due to functional, structural, or evolutionary relationships. Alignments help in comparing sequences to find conserved regions, which are important for understanding biological functions and evolutionary histories.

### Application of Alignments
Alignments are used in many fields. Here are some examples of many:
- **Comparative Genomics:** Identifying conserved sequences across different species
- **Phylogenetics:** Constructing evolutionary trees
- **Evolutionary Studies:** Investigating the evolutionary relationships between organisms
- **Medical Research:** Identifying mutations associated with diseases

### What is a BWT and FM Index

1. **Burrows-Wheeler Transform (BWT):**
   - It groups identical or similar characters closer to each other; useful for pattern matching.

2. **FM-index:**
   - The FM-index is a compressed data structure that can be used to find the number of occurences of a pattern within the compressed text, as well as identify the positions of the occurences.

Overall, these are known to improve sequence alignment efficiency/speed.

# How do I Run it?

## Installation
Clone the repository and navigate to the project directory:
```sh
git clone https://github.com/ADDOURLINKHERe
cd ADD
```

ADD REQUIREMENTS AND VERSIONS HEREEEEEEEEE. Also what python version:
```sh
pip install numpy
```

## Usage

To run the seed and extend BWT aligner tool, use the following command:

```
python seedAndExtend.py [options]
```

**Options:**

- `-k, --k`: Length of k-mer for seeding. Default is 10.
- `-mt, --mismatch_threshold`: Threshold for mismatches allowed during alignment. Default is 3.
- `-d, --gap_penalty`: Penalty for introducing a gap during alignment. Default is -2.
- `-s, --mismatch_penalty`: Penalty for a mismatch during alignment. Default is -1.
- `-m, --match_score`: Score for a match during alignment. Default is 2.
- `-r, --reference_path`: Path to the reference genome file. **Required**.
- `-i, --read_path`: Path to the reads file. **Required**.
- `-o, --output_file`: Path to the output SAM file. **Required**.

## Test Dataset Command:

To run a test example with default parameters, copy the following command:

```
python seedAndExtend.py -r reference.fasta -i reads.fastq -o output.sam
```

### Public Dataset 
The publically available FASTq file we are using is `ERR10021327.fastq`, which was described earlier, while the reference we are using is titled `sequence.fasta`. You can find these files in the `data` directory. Simply copy and modify the names of the files into the command above (from Test Dataset Command). If you have time, give it a try :)

## File Output Format
The SAM output file contains the alignment results of the sequencing reads against the reference genome. Each line represents a single alignment record, and the fields are tab-separated similar to how it appears using `bwa`. Example format:
```plaintext
@HD	VN:1.0	SO:unsorted
@SQ	SN:NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome	LN:0
ERR10021327.3116638	0	NC_045512.2	11625	60	81M	*	0	0	GCTATTTTTGTACTTGTTACTTTGGCCTCTTTTGTTTACTCAACCGCTACTTTAGACTGACTCGTGGTGTTTAGGATTTCTGAGTTTCTACACTGGCGTTTTGATATAGGAGTTCACAGGGAGTACTCCCGCCCGCGAATAGCATAGCTGCGTTCTGCGTGGGGGGTGATTTGTTGGGGGGTGGTGGGGGTGGGGGGTGGGCTGGGGGCGCGGGGGCGGGG	F:::FF:F,:FFFFFFFFFF:FFF:FFFF:FFFF,:FFF:,,,FFFFFF,:FF,:F:FF,:,F:,FF:F::F,,F,:,,:F,,FFF,FFFFFF,FF:F,,:,F,F:F,,F,,,:,:F,FFF,,,,F,,:F,,:F,,:,,,,FF:,:F,FF,,,,,,,,,,,,,,,F,,,,,,F,,,FFF,F:FF,:F,,,,FFF,:,,F::,,FF,FF,,,,F:,F,F,,:	NM:i:4	MD:Z:63T9T4A2T	AS:i:152	XS:i:0
```

Yay! Let's break down some (not all!) parts of this example a little:

1. **Header Section (@HD, @SQ):**
   - `VN:1.0`: SAM format version.
   - `SO:unsorted`: Sorting order of alignments--here, it is unsorted.

2. **Sequence Dictionary (@SQ):**
   - `SN:NC_045512.2`: Reference sequence name 
   - `LN:0`: Length of the reference sequence 

3. **Alignment Record:**
   - `ERR10021327.3116638`: Query template name 
   - `NC_045512.2`: Reference sequence name that the read aligns to
   - `11625`: The index the reference was first found
   - `60`: Mapping quality 
   - `81M`: CIGAR (match/mismatch/indel) string representing the alignment. Here, `81M` means the read is aligned as a match that is 81 bases long.
   - `GCTATTTTTGTACTTGTTACTTTGGCCTCTTTTGTTTACTCAACCGCTACTTTAGACTGACTCGTGGTGTTTAGGATTTCTGAGTTTCTACACTGGCGTTTTGATATAGGAGTTCACAGGGAGTACTCCCGCCCGCGAATAGCATAGCTGCGTTCTGCGTGGGGGGTGATTTGTTGGGGGGTGGTGGGGGTGGGGGGTGGGCTGGGGGCGCGGGGGCGGGG`: Read sequence
   - `NM:i:4`: Number of mismatches in the alignment.
   - `AS:i:152`: Alignment score.

To sum it up, this alignment record represents the alignment of a sequencing read (`ERR10021327.3116638`) to the reference sequence (`NC_045512.2`) starting at position 11625 with a mapping quality of 60. The alignment consists of a match (81 bases) with 4 mismatches at specific positions.
