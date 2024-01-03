#! /usr/bin/bash


rm -rf ./downloads
mkdir -p ./downloads/bin


echo "Downloading and installing blat..."
wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat -O ./downloads/bin/blat
chmod +x ./downloads/bin/blat

echo "Downloading and installing blast..."
wget -q https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary -O ./downloads/bin/bedtools
chmod +x ./downloads/bin/bedtools


echo "Downloading annotation files..."
pip install bedparse
wget -q https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz -O ./downloads/Homo_sapiens.GRCh38.110.gtf.gz
gunzip ./downloads/Homo_sapiens.GRCh38.110.gtf.gz
bedparse gtf2bed ./downloads/Homo_sapiens.GRCh38.110.gtf > ./downloads/Ensembl_Homo_sapiens.GRCh38.110.bed
wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O ./downloads/hg38.fa.gz
gunzip ./downloads/hg38.fa.gz
wget -p https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz -O ./downloads/Homo_sapiens.GRCh38.cds.all.fa.gz
gunzip ./downloads/Homo_sapiens.GRCh38.cds.all.fa.gz

awk '/^>/ {printf("\n%s\t", $1); next; } { printf("%s", $0);} END {printf("\n");}' Homo_sapiens.GRCh38.cds.all.fa | sed '1d' > ./downloads/Homo_sapiens.GRCh38.cds.all.tsv