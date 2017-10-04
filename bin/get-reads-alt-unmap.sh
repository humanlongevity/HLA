#!/bin/bash
set -eu -o pipefail

[[ -z $@ ]] && echo \
"
Filter reads from a .bam file for HLA typing.\n
Usage: ./get-reads-alt-unmap.sh <input.bam> <output.bam>"

IN="$1"
OUT="$2"

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    NCORE=$(grep -c '^processor' /proc/cpuinfo)
    # POSIX compliant system memory query
    MEM=$(awk 'BEGIN{for (i=1; i<ARGC;i++)
       printf "%.0fGB\n", ARGV[i]}' $(echo "$(grep MemTotal /proc/meminfo | awk '{print $2}') * 0.95 / 1000 / 1000" | bc))
elif [[ "$OSTYPE" == darwin* ]]; then
    NCORE=$(sysctl -n hw.ncpu)
    MEM=8GB
fi

echo "NCORE=$NCORE MEM=$MEM"

BIN="`dirname \"$0\"`"

cleanup(){

  rm "${OUT}.1.fq" || true
  rm "${OUT}.2.fq" || true
  rm "${OUT}.full.bam" || true
  rm "${OUT}.full.bam.bai" || true

}

trap cleanup EXIT

sambamba view \
  -f "bam" -h -p -l 0 -t $NCORE \
  -F "unmapped or mate_is_unmapped or (ref_name == 'chr6' and (position > 29844528 and position < 33100696)) or ref_name =~ /^HLA|chr6.*alt/" \
  $IN | sambamba sort -p -n -t $NCORE -o - /dev/stdin |
  bamToFastq -i /dev/stdin -fq "${OUT}.1.fq" -fq2 "${OUT}.2.fq"

bwa mem -t $NCORE "$BIN/../data/chr6/hg38.chr6.fna" "${OUT}.1.fq" "${OUT}.2.fq" | samtools view -b - | sambamba sort -m $MEM -t $NCORE -o "${OUT}.full.bam" /dev/stdin

sambamba index "${OUT}.full.bam"

samtools view -o "$OUT" -b "${OUT}.full.bam" chr6:29844528-33100696
sambamba index "$OUT"

# S3 upload
[[ $# -eq 3 ]] && {
	aws s3 cp ${OUT} $3 --sse
	aws s3 cp ${OUT}.bai $3 --sse
}
