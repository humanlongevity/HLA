#!/bin/bash
set -ex

IN=$1
OUT=$2
NCORE=`grep -c '^processor' /proc/cpuinfo`
BIN="`dirname \"$0\"`"

samtools view -u $IN \
	chr6:29844528-33100696 chr6_GL000250v2_alt chr6_GL000251v2_alt chr6_GL000252v2_alt chr6_GL000253v2_alt chr6_GL000254v2_alt chr6_GL000255v2_alt chr6_GL000256v2_alt chr6_GL383533v1_alt chr6_KB021644v2_alt chr6_KI270758v1_alt chr6_KI270797v1_alt chr6_KI270798v1_alt chr6_KI270799v1_alt chr6_KI270800v1_alt chr6_KI270801v1_alt chr6_KI270802v1_alt \
	| sambamba sort -p -n -t $NCORE -o /dev/stdout /dev/stdin \
	| bamToFastq -i /dev/stdin -fq ${OUT}.1.fq -fq2 ${OUT}.2.fq
bwa mem -t $NCORE $BIN/../data/chr6/hg38.chr6.fna ${OUT}.1.fq ${OUT}.2.fq | samtools view -b - | sambamba sort -m 60GB -t $NCORE -o ${OUT}.full.bam /dev/stdin
sambamba index ${OUT}.full.bam
samtools view -o $OUT -b ${OUT}.full.bam chr6:29844528-33100696
sambamba index $OUT
rm ${OUT}.1.fq
rm ${OUT}.2.fq
rm ${OUT}.full.bam
rm ${OUT}.full.bam.bai

[[ $# -eq 3 ]] && {
	aws s3 cp ${OUT} $3 --sse
	aws s3 cp ${OUT}.bai $3 --sse
}
