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

- [BWT Sequence Aligner](#bwt-sequence-aligner-cse185-sp24-)
  - [Overview](#overview)
  - [Features](#features)
  - [Table of Contents](#table-of-contents)
  - [Contributors](#contributors)
  - [File Descriptions](#file-descriptions)
    - [Source Code](#source-code)
    - [Example Data](#example-data)
    - [Example Results](#example-results)
  - [So...what is an alignment](#sowhat-is-an-alignment)
    - [Application of Alignments](#application-of-alignments)
    - [What is a BWT and an FM Index?](#what-is-a-bwt-and-an-fm-index)
  - [How do I run it?](#how-do-i-run-it)
    - [Installation](#installation)
    - [Usage](#usage)
    - [File Output Format](#file-output-format)

[](#bwt-sequence-aligner-cse185-sp24-)

## Contributors

- Gary Lin (A16915179)
- Lorenzo Olmo Marchal (A17013640)
- Nabeeha Rashid (A16851960)

This project is intended to be our final project in CSE 185: Advanced Bioinformatics Laboratory for the Spring 2024 quarter, to be presented to our peers and Professor Melissa Gyrmrek.

[🔼 Back to top](#bwt-sequence-aligner-cse185-sp24-)

## File Descriptions

File Tree:

```txt
│   README.md
|   requirements.txt
├───benchmarking
├───data
│       chimpref.fa
│       ERR10021327.fastq
│       sequence.fasta
├───results
│       outputcovid.sam
└───src
        reference_indexing.py
        seed_and_extend.py
```

Here is a brief description of our files:

[`README.md`](/README.md): This file you're reading! 👋 We hope you're enjoying it.

[`requirements.txt`](/requirements.txt): The modules one is *required* to install in order to run this tool. More details about installation can be found [here](#installation).

### [Source Code](/src/)

- [`reference_indexing.py`](/src/reference_indexing.py): implements the indexing of a reference genome sequence using the Burrows-Wheeler Transform (BWT) and the FM-index.
- [`seed_and_extend.py`](/src/seed_and_extend.py): implements local alignment, first identifying potential seed matches between the query sequence and the reference genome using exact k-mer matching, then extending these seeds locally to find the most optimal alignment.

### [Example Data](/data/)

- [`chimpref.fa`](/data/chimpref.fa):
- [`ERR10021327.fastq`](/data/ERR10021327.fastq): a publically available query file containing the reads to be aligned, from the [ENA project PRJEB37886](https://www.ebi.ac.uk/ena/browser/view/PRJEB37886).
- [`sequence.fasta`](/data/sequence.fasta): a publically available reference genome file from [NCBI accession ERR10021327](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2).

### [Example Results](/results/)

- [`outputcovid.sam`](/results/outputcovid.sam): the result of running `seed_and_extend.py` on our example data files. Specifically, we used the following command (from the root directory):

```shell
python ./src/seed_and_extend.py -r ./data/sequence.fasta -i ./data/ERR10021327.fastq -o ./results/outputcovid.sam
```

### [Benchmarking Analysis](/benchmarking/)

Read more about how our tool compares to the standard tools out there: [Why Should I Run It?](#why-should-i-run-it)

[🔼 Back to top](#bwt-sequence-aligner-cse185-sp24-)

## So...what is an alignment?

An alignment is a process of arranging DNA, RNA, or protein sequences to identify regions of similarity. These similarities can be due to functional, structural, or evolutionary relationships. Alignments help in comparing sequences to find conserved regions, which are important for understanding biological functions and evolutionary histories.

On a simple level, suppose you wanted to align the sequences `ATCG` and `TTCG`.

```text
ATCG                       _ATCG
 |||  one mismatch           |||  two indels (insertions/deletions)
TTCG                       T_TCG
```

Both are valid alignments, but one may be more biologically accurate. Typically alignment generates "scores" given some metric, which may reward or penalize certain mutations more than others to approach biological accuracy. The exact metric one may use varies depending on the exact situation, which is why our tool allows one to [specify](#options) how much to reward or penalize the mutation type of their choice.

### Application of Alignments

Alignments are used in many fields in bioinformatics and system biology:

- **Comparative Genomics:** Identifying conserved sequences across different species
- **Phylogenetics:** Constructing evolutionary trees
- **Evolutionary Studies:** Investigating the evolutionary relationships between organisms
- **Medical Research:** Identifying mutations associated with diseases

### What is a BWT and an FM Index?

- **Burrows-Wheeler Transform (BWT):** a data structure which groups identical or similar characters closer to each other; useful for pattern matching and compression.
- **FM-index:**  a compressed data structure that can be used to find the number of occurences of a pattern within the compressed text, as well as identify the positions of the occurences.

Overall, these data structures are known to improve sequence alignment efficiency/speed.

[🔼 Back to top](#bwt-sequence-aligner-cse185-sp24-)

## How Do I Run It?

### Installation

The lovely UI/UX and front-end engineers at GitHub have bestowed upon us a glorious green button that says `<> Code`, where one may choose a plethora of options to build this tool. To avoid installing locally, one may open a GitHub Codespace. Otherwise, one may choose to clone the repository and navigate to the project directory:

```sh
git clone https://github.com/lolmomarchal/Alignment-CSE185-SP24-.git
cd ./Alignment-CSE185-SP24-.git
```

ADD REQUIREMENTS AND VERSIONS HEREEEEEEEEE. Also what python version:

```sh
pip install -r requirements.txt
```

### Usage

To run the seed and extend BWT aligner tool, use the following command in the `src` directory:

```sh
python seed_and_extend.py [options]
```

#### Options

- `-h, --help`: Shows a help message.
- `-k, --k`: Length of k-mer for seeding. Default is 10.
- `-mt, --mismatch_threshold`: Threshold for mismatches allowed during alignment. Default is 3.
- `-d, --gap_penalty`: Penalty for introducing a gap during alignment. Default is -2.
- `-s, --mismatch_penalty`: Penalty for a mismatch during alignment. Default is -1.
- `-m, --match_score`: Score for a match during alignment. Default is 2.
- `-r, --reference_path`: Path to the reference genome file. **Required**.
- `-i, --read_path`: Path to the reads file. **Required**.
- `-o, --output_file`: Path to the output SAM file. **Required**.

#### Test Dataset Command

To run a test example with default parameters, copy and run the following command in the `src` directory replacing `<REFERENCE>` and `<READS>` with your reference sequence and your reads file, respectively:

```shell
python seed_and_extend.py -r <REFERENCE>.fasta -i <READS>.fastq -o output.sam
```

### Public Dataset

The publically available FASTQ file we are using is `ERR10021327.fastq`, which was described earlier, while the reference we are using is titled `sequence.fasta`. You can find these files in the `data` directory. Simply copy and modify the names of the files into the command above (from Test Dataset Command). If you have time, give it a try! 😊

### File Output Format

The SAM output file contains the alignment results of the sequencing reads against the reference genome. Each line represents a single alignment record, and the fields are tab-separated similar to how it appears using `bwa`. Example format:

```plaintext
@HD VN:1.0 SO:unsorted
@SQ SN:NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome LN:0
ERR10021327.3116638 0 NC_045512.2 11625 60 81M * 0 0 GCTATTTTTGTACTTGTTACTTTGGCCTCTTTTGTTTACTCAACCGCTACTTTAGACTGACTCGTGGTGTTTAGGATTTCTGAGTTTCTACACTGGCGTTTTGATATAGGAGTTCACAGGGAGTACTCCCGCCCGCGAATAGCATAGCTGCGTTCTGCGTGGGGGGTGATTTGTTGGGGGGTGGTGGGGGTGGGGGGTGGGCTGGGGGCGCGGGGGCGGGG F:::FF:F,:FFFFFFFFFF:FFF:FFFF:FFFF,:FFF:,,,FFFFFF,:FF,:F:FF,:,F:,FF:F::F,,F,:,,:F,,FFF,FFFFFF,FF:F,,:,F,F:F,,F,,,:,:F,FFF,,,,F,,:F,,:F,,:,,,,FF:,:F,FF,,,,,,,,,,,,,,,F,,,,,,F,,,FFF,F:FF,:F,,,,FFF,:,,F::,,FF,FF,,,,F:,F,F,,: NM:i:4 MD:Z:63T9T4A2T AS:i:152 XS:i:0
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
   - `NM:i:4`: Number of mismatches in the alignment
   - `AS:i:152`: Alignment score

To sum it up, this alignment record represents the alignment of a sequencing read (`ERR10021327.3116638`) to the reference sequence (`NC_045512.2`) starting at position 11625 with a mapping quality of 60. The alignment consists of a match (81 bases) with 4 mismatches at specific positions.

Read more about the SAM file format [here](https://samtools.github.io/hts-specs/SAMv1.pdf).

[🔼 Back to top](#bwt-sequence-aligner-cse185-sp24-)

## Why Should I Run It?

[🔼 Back to top](#bwt-sequence-aligner-cse185-sp24-)
