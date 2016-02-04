#!/bin/bash
set -e

S3=$1
ID=$2
OUT=hla-$ID

[[ $# -ne 2 ]] && {
    echo "usage: $(basename "$0") [S3://path.bam] [sample_id]";
    exit 1;
}
BIN="`dirname \"$0\"`"

mkdir -p $OUT

echo "Extracting reads from S3"
samtools view -u $S3 chr6:29886751-33090696 | samtools view -L $BIN/../data/hla.bed - | $BIN/preprocess.pl > $OUT/$ID.fq
echo "Aligning reads to IMGT database"
diamond blastx -C 20000 --index-mode 2 --seg no --min-score 10 --top 20 -c 1 -v -d $BIN/../data/hla -q $OUT/$ID.fq -a $OUT/$ID.daa
echo "Preparing data for the typing algorithm"
$BIN/by.msa.pl $OUT/$ID.fq $OUT/$ID.daa $OUT/$ID.tsv
echo "Typing"
$BIN/typing.r $OUT/$ID.tsv $OUT/$ID.hla
echo "Reporting"
$BIN/report.py -in $OUT/$ID.hla -out $OUT/$ID.json -subject $ID -sample $ID
