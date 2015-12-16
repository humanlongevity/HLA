#!/usr/bin/env Rscript
args <- commandArgs(T)
library(parallel)
options(mc.cores = detectCores())
library(data.table)

dt <- fread(args[1])
setnames(dt, c('q', 'qpos', 't', 'tlen', 'ts', 'te', 'type', 'msa', 'exon', 'specific', 'left', 'right', 'start', 'end'))
print(nrow(dt))

# for HLA alleles with frame shift variants, we require reads span over the frame shift site
frame.shift <- fread('data/hla.shift')
setnames(frame.shift, c('t', 'exon', 'shift'))
frame.shift[, type := sub('-E.+', '', t)]
frame.shift[, t := NULL]
setkey(frame.shift, type, exon)
frame.shift <- unique(frame.shift)
setkey(dt, type, exon)
dt <- frame.shift[dt]
print(nrow(dt))
spanned <- dt[ts < shift-1 & te > shift+1]
print(nrow(spanned))
print(sort(table(spanned$type)))

dt <- dt[type %in% spanned$type | !(type %in% frame.shift$type)]
print(nrow(dt))

# filter pair end matches
dt[, qp := sub('/\\d$', '', q)]
nr <- dt[, .(pair1 = length(unique(q))), keyby = qp]
setkey(dt, qp)
dt <- nr[dt]
nr <- dt[, .(pair2 = length(unique(q))), keyby = .(qp, t)]
setkey(dt, qp, t)
dt <- nr[dt]
dt <- dt[pair1 == pair2]

# filter non-specific matching
#dt <- dt[specific==1 & left==0 & right==0]
dt <- dt[left==0 & right==0]
# TODO, filter core exons

#library(IRanges)
#setkey(dt, t)
#cov <- dt[, .(
#	n = .N, 
#	cov = sum(width(reduce(IRanges(pos, width = len)))),
#), keyby = t]

mat <- dcast(dt, q ~ type, value.var = 'specific', fun.aggregate = max, fill = 0)
qs <- mat$q
mat[, q := NULL]
allele.names <- colnames(mat)
allele.genes <- unique(sub('\\*.+', '', allele.names))
alleles <- 1:ncol(mat)
reads <- 1:nrow(mat)
na <- length(alleles)
nr <- length(reads)

library(lpSolve)
f.obj <- c(rep(0, na), rep(1, nr))

# constraints for hit incidence matrix
f.con <- do.call(rbind, mclapply(1:nr, function(i){
    con <- c(mat[i,], rep(0, nr))
    con[na + i] <- -1
    con
}))
f.dir <- rep('>=', nr)
f.rhs <- rep(0, nr)

# constraints for 1 or 2 alleles per gene
f.con2 <- do.call(rbind, mclapply(allele.genes, function(gene){
    con <- c(rep(0, na), rep(0, nr))
    con[grep(sprintf("^%s", gene), allele.names)] <- 1
    rbind(con, con)
}))
f.dir2 <- rep(c('>=', '<='), length(allele.genes))
f.rhs2 <- rep(c( 1,    2  ), length(allele.genes))

# final constraints
f.con <- rbind(f.con, f.con2)
f.dir <- c(f.dir, f.dir2)
f.rhs <- c(f.rhs, f.rhs2)

lps <- lp('max', f.obj, f.con, f.dir, f.rhs, all.bin = T)
solution <- lps$solution[alleles]
names(solution) <- allele.names
solution <- solution[order(-solution)]
solution <- solution[solution > 0]
print(solution)

save.image(file = 'temp.rda')
