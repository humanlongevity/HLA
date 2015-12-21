#!/usr/bin/env Rscript
args <- commandArgs(T)
library(parallel)
options(mc.cores = detectCores())
library(data.table)

dt <- fread(args[1])
setnames(dt, c('q', 'qpos', 't', 'tlen', 'ts', 'te', 'type', 'msa', 'exon', 'specific', 'left', 'right', 'start', 'end'))

# for HLA alleles with frame shift variants, we require reads span over the frame shift site
frame.shift <- fread('data/hla.shift')
setnames(frame.shift, c('t', 'exon', 'shift'))
frame.shift[, type := sub('-E.+', '', t)]
frame.shift[, t := NULL]
setkey(frame.shift, type, exon)
frame.shift <- unique(frame.shift)
setkey(dt, type, exon)
dt <- frame.shift[dt]
spanned <- dt[ts < shift-1 & te > shift+1]

dt <- dt[type %in% spanned$type | !(type %in% frame.shift$type)]

# only keep Class I (A, B, C) and DRB1, DQB1, and DPB1.
dt <- dt[msa %in% c('ClassI', 'DRB1', 'DQB1', 'DPB1')]

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
mat <- as.matrix(mat)

# filter out types with too few reads
counts <- colSums(mat)
summary(counts)
cand <- counts > quantile(counts, 0.25) 
mat <- mat[, cand]
## filter out reads with no alleles mapped to
counts <- rowSums(mat)
summary(counts)
cand <- counts > 0
mat <- mat[cand, ]
qs <- qs[cand]

allele.names <- colnames(mat)
allele.genes <- unique(sub('\\*.+', '', allele.names))
n.genes <- length(allele.genes)
alleles <- 1:ncol(mat)
reads <- 1:nrow(mat)
na <- length(alleles)
nr <- length(reads)
gamma <- 0.01
beta <- 0.009

library(lpSolve)
f.obj <-  c(rep(-gamma, na), rep(1,  nr), rep(-beta, nr), 0  )
f.type <- c(rep('b', na),    rep('b',nr), rep('i',   nr), 'i')

all.zero <- c(rep(0, na), rep(0, 2 * nr), 0)
heter <- length(all.zero)

# constraints for 1 or 2 alleles per gene
#f.con.bound <- do.call(rbind, mclapply(allele.genes, function(gene){
#    con <- all.zero
#    con[grep(sprintf("^%s", gene), allele.names)] <- 1
#    rbind(con, con)
#}))
#f.dir.bound <- rep(c('>=', '<='), n.genes)
#f.rhs.bound <- rep(c( 1,    2  ), n.genes)
f.con.bound <- t(matrix(all.zero, nrow = heter, ncol = n.genes))
for(g in seq_along(allele.genes)){
	this.gene <- grep(sprintf("^%s", allele.genes[g]), allele.names)
	f.con.bound[g, this.gene] <- 1
}
f.dir.bound <- rep(c('>=', '<='), each = n.genes)
f.rhs.bound <- rep(c( 1,    2  ), each = n.genes)

zero.m <- t(matrix(all.zero, nrow = heter, ncol = nr))
yindex <- matrix(c(1:nr, na + 1:nr), ncol = 2)
gindex <- matrix(c(1:nr, na + nr + 1:nr), ncol = 2)
# constraints for hit incidence matrix
#system.time(
#f.con.hit <- do.call(rbind, mclapply(1:nr, function(i){
#    con <- all.zero
#    con[1:na] <- mat[i,]
#    con[na + i] <- -1
#    con
#}))
#)
system.time(f.con.hit <- zero.m)
system.time(f.con.hit[yindex] <- -1)
system.time(f.con.hit[, 1:na] <- mat)
f.dir.hit <- rep('>=', nr)
f.rhs.hit <- rep(0, nr)


# num of heterozygous genes
f.con.heter <- all.zero
f.con.heter[1:na] <- 1
f.con.heter[heter] <- -1
f.dir.heter <- '=='
f.rhs.heter <- n.genes

# regularization for heterozygous genes
#system.time(f.con.reg <- do.call(rbind, mclapply(1:nr, function(i){
#    con <- rbind(all.zero, all.zero, all.zero)
#    con[ , na + nr + i] <- 1
#    con[1, na + i] <- -n.genes
#    con[3, na + i] <- -n.genes
#    con[2, heter] <- -1
#    con[3, heter] <- -1
#    con
#})))
#f.dir.reg <- rep(c('<=', '<=', '>='), nr)
#f.rhs.reg <- rep(c(0, 0, -n.genes), nr)

zero.m[gindex] <- 1
f.con.reg1 <- zero.m
f.con.reg2 <- zero.m
f.con.reg3 <- zero.m
f.con.reg2[, heter] <- -1
f.con.reg3[, heter] <- -1
f.con.reg1[yindex] <- -n.genes
f.con.reg3[yindex] <- -n.genes
f.dir.reg <- rep(c('<=', '<=', '>='), each = nr)
f.rhs.reg <- rep(c(0, 0, -n.genes), each = nr)

# final constraints
f.con <- rbind(f.con.hit, f.con.bound, f.con.bound, f.con.heter, f.con.reg1, f.con.reg2, f.con.reg3)
f.dir <- c(f.dir.hit, f.dir.bound, f.dir.heter, f.dir.reg)
f.rhs <- c(f.rhs.hit, f.rhs.bound, f.rhs.heter, f.rhs.reg)

#save.image('temp.rda')

system.time(lps <- lp('max', f.obj, f.con, f.dir, f.rhs, int.vec = which(f.type == 'i'), binary.vec = which(f.type == 'b')))
solution <- lps$solution[alleles]
names(solution) <- allele.names
solution <- solution[order(-solution)]
solution <- solution[solution > 0]
solution <- names(solution)
print(solution)

explained <- sum(apply(mat[, solution], 1, max))
delta <- explained - sapply(seq_along(solution), function(i) sum(apply(mat[, solution[-i]], 1, max)))
importance <- delta / explained * length(solution)
solution <- data.frame(solution, importance)
print(solution)
write.table(solution, row = F, col = F, sep = '\t', quo = F, file = args[2])

#save.image(file = sprintf('%s.temp.rda', args[2]))
