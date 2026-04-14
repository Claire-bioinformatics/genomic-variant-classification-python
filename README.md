# genomic-variant-filter-python
# Genomic Variant Classification Tool (Python)

## Overview
This project implements a Python command-line tool for processing genomic variant data from VCF, GFF, and FASTA files. The script filters variants based on quality scores and classifies high-quality variants as non-coding, synonymous, or non-synonymous by integrating genomic annotations and sequence information.

## Project Objective
The aim of this project was to build a reproducible bioinformatics workflow that:

- reads genomic variants from a VCF file  
- uses genome annotation from a GFF file  
- uses a reference genome FASTA file  
- filters variants based on QUAL score  
- determines whether variants occur in coding regions  
- predicts whether coding variants alter the encoded amino acid sequence  

## Inputs
The script takes the following input files via a command-line interface:

- VCF file (`--vcffile`)  
- GFF file (`--gfffile`)  
- genome FASTA file (`--fastafile`)  

All three inputs are required.

## Analysis Workflow
For each variant in the VCF file, the script:

1. checks whether the QUAL score is greater than 20  
2. counts variants with QUAL ≤ 20 and records this in a log file  
3. for variants with QUAL > 20:
   - determines whether the variant lies in a coding region  
   - classifies coding variants as:
     - **synonymous** (no amino acid change)  
     - **non-synonymous** (amino acid change)  
   - labels variants outside coding regions as **non-coding**  

## Outputs

### 1. Log File
The log file contains:
- confirmation of input filenames  
- count of variants with QUAL ≤ 20  
- output file locations  
- error messages if the script fails  

### 2. Bar Plot
A bar plot showing the proportion of high-quality variants (QUAL > 20) classified as:
- Non-coding  
- Synonymous  
- Non-synonymous  

### 3. Tab-Separated Results Table
For each variant with QUAL > 20, the output table includes:

- Chrom  
- Pos  
- Ref  
- Alt  
- Type (Non-coding / Synonymous / Non-synonymous)  
- Transcript  
- Protein Location  
- Ref AA  
- Alt AA  

## Example Classification
- **Non-coding**: variant is not located in a coding region  
- **Synonymous**: variant is in a coding region but does not change the amino acid  
- **Non-synonymous**: variant is in a coding region and changes the amino acid  

## How to Run

Run the script using:

```bash
python3 script_name.py --vcffile input.vcf --gfffile annotation.gff --fastafile genome.fasta
