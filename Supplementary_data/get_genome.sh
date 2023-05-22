#!/bin/bash

# Download Genome frome ensembl
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa hg38.fa
#samtools faidx hg38.fa

# Donwload annotation from ensembl
wget ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
gzip -d Homo_sapiens.GRCh38.84.gtf.gz
mv Homo_sapiens.GRCh38.84.gtf hg38.gtf


bwa index hg38.fa