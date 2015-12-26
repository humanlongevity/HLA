#!/bin/bash
S3=$1
OUT=$2

samtools view -u $S3 chr6:29723339-33089696 | samtools view -u -L data/hg38.hla_region.refGene.bed - | bamtools convert -format fastq | sga preprocess -q20 -f5 -m 100 - > $OUT.fq
#diamond blastx -C 10000 --top 10 -s 1 -c 1 -v -d data/hla -q $OUT.fq -a $OUT.daa 
diamond blastx -C 10000 --seg no --min-score 30 --top 3 -s 1 -c 1 -v -d data/hla -q $OUT.fq -a $OUT.daa
./bin/by.msa.pl $OUT.fq data/hla.tsv $OUT.daa 2>$OUT.log >$OUT.tsv
./bin/typing.r $OUT.tsv $OUT.hla
cat $OUT.hla
