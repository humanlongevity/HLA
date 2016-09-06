#!/usr/bin/env Rscript
# Author: Xie Chao

args <- commandArgs(T)
if(length(args) != 3) {
	cat("usage: full.r input.tsv.dna input.hla output.tsv\n")
	q()
}
align.path <- args[1]
hla.path <- args[2]
out.path <- args[3]
library(data.table)

dna <- fread(align.path)
setnames(dna, c('q', 'type4', 'type6'))
type <- fread(hla.path)
solution <- type[rank == 0, solution]
dna[, type := sub('-.+', '', type4)]
dna[, exon := sub('.+-(E\\d+)-.+', '\\1', type4)]
dna[, full := sub('-.+', '', type6)]
exons <- unique(dna[, .(full, exon)])[, .(total = .N), by = full]
setkey(exons, full)
dna <- dna[type %in% solution]
setkey(dna, type, exon)
dna <- dna[, .N, by = .(type, exon, full)]
dna <- dna[, .SD[rank(-N, tie = 'min') == 1], by = .(type, exon)]
dna <- dna[, .(good = length(exon)), by = .(type, full)]
setkey(dna, full)
dna <- exons[dna]
dna[, rank := rank(-good/total, tie = 'min'), keyby = type]
dna <- dna[, full := paste(full[rank==1], collapse = ';'), by = type]
dna <- dna[, .SD[1], by = type]
write.table(dna, row = F, sep = '\t', quo = F, file = out.path)
