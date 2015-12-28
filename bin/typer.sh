#!/bin/bash
S3=$1
OUT=$2

set -ex

#samtools view -u $S3 chr6:29723339-33089696 | samtools view -u -L data/hg38.hla_region.refGene.bed - | bamtools convert -format fastq | sga preprocess -q20 -f5 -m 100 - > $OUT.fq
samtools view -u $S3 chr6:26085440-33091696 | samtools view -u -L data/hla.bed - | bamtools convert -format fastq | sga preprocess -q20 -f5 -m 100 - > $OUT.fq
#diamond blastx -C 10000 --top 10 -s 1 -c 1 -v -d data/hla -q $OUT.fq -a $OUT.daa 
#diamond blastx -C 10000 --seg no --min-score 30 --top 3 -s 1 -c 1 -v -d data/hla -q $OUT.fq -a $OUT.daa
diamond blastx -C 20000 --index-mode 2 --seg no --min-score 10 --top 5 -c 1 -v -d data/hla -q $OUT.fq -a $OUT.daa
./bin/by.msa.pl $OUT.fq data/hla.tsv $OUT.daa 2>$OUT.log >$OUT.tsv
./bin/typing.r $OUT.tsv $OUT.hla
cat $OUT.hla
