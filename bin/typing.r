#!/usr/bin/env Rscript

args <- commandArgs()
code.source <- sub('--file=', '', args[4])
if(length(args) != 7) {
	cat("usage:", code.source, "input.tsv output.tsv\n")
	q()
}
data.dir <- dirname(code.source)
align.path <- args[6]
out.path <- args[7]
library(parallel)
options(mc.cores = detectCores())
library(data.table)
library(lpSolve)

all <- fread(align.path)
setnames(all, c('q', 'qpos0', 'qpos', 't', 'tlen', 'ts', 'te', 'mis', 'type', 'msa', 'exon', 'specific', 'left', 'right', 'start', 'end'))
all[, specific := as.double(specific)]
all <- all[mis == 0]

# for HLA alleles with frame shift variants, we require reads span over the frame shift site
frame.shift <- fread(sprintf('%s/../data/hla.shift', data.dir))
setnames(frame.shift, c('t', 'EXON', 'shift'))
frame.shift[, type := sub('-E.+', '', t)]
frame.shift[, MSA := sub('\\*.+', '', t)]
frame.shift[, MSA := ifelse(MSA %in% c('A', 'B', 'C'), 'ClassI', MSA)]
setkey(all, msa, exon, ts)
frame.shift[, mapped := all[msa == MSA & exon == EXON & ts < shift-1 & te > shift+1, .N], by = .(MSA, EXON, shift)]
frame.shift <- frame.shift[mapped > 0]
frame.shift[, t := NULL]
setkey(frame.shift, type, EXON)
frame.shift <- unique(frame.shift)
setkey(all, type, exon)
all <- frame.shift[all]
spanned <- all[ts < shift-1 & te > shift+1]
all <- all[type %in% spanned$type | !(type %in% frame.shift$type)]

# exon data avalability in IMGT
exons <- data.table(read.delim(sprintf('%s/../data/hla.tsv', data.dir), h = F, stringsAsFactor = F)[, c('V2', 'V3')])
setnames(exons, c('type', 'EXON'))

# filter non-specific matching
#all <- all[specific==1 & left==0 & right==0]
all <- all[left==0 & right==0]

# filter pair end matches
all[, qp := sub('/\\d$', '', q)]
nr <- all[, .(pair1 = length(unique(q))), keyby = qp]
setkey(all, qp)
all <- nr[all]
nr <- all[, .(pair2 = length(unique(q))), keyby = .(qp, type)]
setkey(all, qp, type)
all <- nr[all]
all <- all[pair1 == pair2]
specific.pairs <- all[pair1 == 2 & specific == 1, qp]
all[qp %in% specific.pairs, specific := 1]

# only keep Class I (A, B, C) and DRB1, DQB1, and DPB1
keep <- c('ClassI', 'DRB1', 'DQB1', 'DPB1')
ignore <- all[! msa %in% keep, q]
all <- all[msa %in% keep]
all[q %in% ignore, specific := 0]
#all[q %in% ignore, specific := 0.1]
#all[specific == 0, specific := 0.1]
all <- all[specific > 0]


# separate core and non-core exons
print(nrow(all))
non.core <- all[!((msa == 'ClassI' & EXON %in% c('E2', 'E3')) | (msa != 'ClassI' & EXON == 'E2'))]
print(nrow(non.core))
core <- all[(msa == 'ClassI' & EXON %in% c('E2', 'E3')) | (msa != 'ClassI' & EXON == 'E2')]
print(nrow(core))

# matrix for core exon alignments
mat <- dcast(core, q ~ type, value.var = 'specific', fun.aggregate = max, fill = 0)
qs <- mat$q
mat[, q := NULL]
mat <- as.matrix(mat)
weight <- apply(mat, 1, max)
mat2 <- mat
mat[mat > 0] <- 1

# matrix for non-core exon alignments
mat.noncore <- dcast(non.core, q ~ type, value.var = 'specific', fun.aggregate = max, fill = 0)
mat.noncore[, q := NULL]
mat.noncore <- as.matrix(mat.noncore)

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
weight <- weight[cand]


# preparing to setup Integer Linear Programming
allele.names <- colnames(mat)
allele.genes <- unique(sub('\\*.+', '', allele.names))
n.genes <- length(allele.genes)
alleles <- 1:ncol(mat)
reads <- 1:nrow(mat)
na <- length(alleles)
nr <- length(reads)
all.zero <- c(rep(0, na), rep(0, 2 * nr), 0)
heter <- length(all.zero)
zero.m <- t(matrix(all.zero, nrow = heter, ncol = nr))
yindex <- matrix(c(1:nr, na + 1:nr), ncol = 2)
gindex <- matrix(c(1:nr, na + nr + 1:nr), ncol = 2)
gamma <- 0.01
beta <- 0.009

# ILP objective
f.obj <-  c(rep(-gamma, na), weight,      -beta * weight, 0  )
f.type <- c(rep('b', na),    rep('b',nr), rep('i',   nr), 'i')

# constraints for number of chromosomes
f.con.bound <- t(matrix(all.zero, nrow = heter, ncol = n.genes))
for(g in seq_along(allele.genes)){
	this.gene <- grep(sprintf("^%s", allele.genes[g]), allele.names)
	f.con.bound[g, this.gene] <- 1
}
f.dir.bound <- rep(c('>=', '<='), each = n.genes)
f.rhs.bound <- rep(c( 1,    2  ), each = n.genes)


# constraints for hit incidence matrix
f.con.hit <- zero.m
f.con.hit[yindex] <- -1
f.con.hit[, 1:na] <- mat
f.dir.hit <- rep('>=', nr)
f.rhs.hit <- rep(0, nr)

# num of heterozygous genes
f.con.heter <- all.zero
f.con.heter[1:na] <- 1
f.con.heter[heter] <- -1
f.dir.heter <- '=='
f.rhs.heter <- n.genes

# regularization for heterozygous genes
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

# 1: initial solution candidate set via integer linear programming
system.time(lps <- lp('max', f.obj, f.con, f.dir, f.rhs, int.vec = which(f.type == 'i'), binary.vec = which(f.type == 'b')))
solution <- lps$solution[alleles]
names(solution) <- allele.names
solution <- solution[order(-solution)]
solution <- solution[solution > 0]
solution <- names(solution)
print(solution)

# 2: round 0 for better candidate searching
max.hit <- sum(apply(mat2[, solution], 1, max))
more <- do.call(rbind, mclapply(solution, function(s){
	minus1 <- solution[solution != s]
	minus1.hit <- apply(mat2[, minus1], 1, max)
	gene <- sub('(.+?)\\*.+', '^\\1', s)
	others <- allele.names[grepl(gene, allele.names)]
	other.hit <- sapply(others, function(i) sum(pmax(minus1.hit, mat2[, i])))
	cand <- data.table('rank' = 0, 'solution' = '', 'competitor' = as.character(others), 'missing' = max.hit - other.hit)
    cand <- cand[order(missing)]
	cand <- cand[1:max(sum(missing <= 2), min(50, length(others)))]

	ambig <- cand[missing <= 2, competitor]
	missing <- cand[missing <= 2, missing]
	names(missing) <- ambig
	sol <- s
	if(length(ambig) > 1){
		bests <- sort(ambig)
		x <- as.integer(sub('.+?\\*(\\d+):.+', '\\1', bests))
		y <- as.integer(sub('.+?\\*\\d+:(\\d+).*', '\\1', bests))
		bests <- bests[order(x * 1e5 + y)]
		total <- colSums(mat2[, bests, drop = F])
		in.noncore <- bests %in% colnames(mat.noncore)
		names(in.noncore) <- bests
		total.noncore <- rep(0, length(bests))
		total.noncore[in.noncore] <- colSums(mat.noncore[, bests[in.noncore], drop = F])
		missing <- missing[bests]
		bests <- bests[order(missing * 5 - total - total.noncore / 3)]

		sol <- bests[1]
		all <- unique(c(bests, cand$competitor))
		cand <- cand[match(all, competitor)]
	}
	cand[, rank := 1:nrow(cand)]
	cand[, solution := sol]
	return(copy(cand))
}))

# functions for downstream ranking
shared.exon <- function(a, b){
	exon1 <- exons[type == a, EXON]
	exon2 <- exons[type == b, EXON]
	return(exon1[exon1 %in% exon2])
}
all.exons <- exons[, unique(EXON)]
specific.reads <- function(a, superset, EXONS = all.exons, dt = all){
	superset <- superset[superset != a]
	EXONS <- EXONS[EXONS %in% dt$EXON]
	if(length(EXONS) >= 1){
		my.reads <- dt[type == a & EXON %in% EXONS]
		other.reads <- dt[type %in% superset & EXON %in% EXONS, q]
		return(my.reads[!q %in% other.reads])
	}else{
		return(dt[0])
	}
}
diff.count <- function(a, b, superset = NULL, dt = all){
	exon.shared <- shared.exon(a, b)
	superset <- superset[!superset %in% c(a, b)]
	a.reads <- specific.reads(a, superset, exon.shared, dt)
	b.reads <- specific.reads(b, superset, exon.shared, dt)
	diff1 <- a.reads[!q %in% b.reads$q, .N]
	diff2 <- b.reads[!q %in% a.reads$q, .N]
	return(c(nrow(a.reads), nrow(b.reads), diff1, diff2))
}
get.diff <- function(sol, comp, superset){
	superset <- superset[!superset %in% c(sol, comp)]
	diff.all <- diff.count(sol, comp)
	diff.core <- diff.count(sol, comp, superset, core)
	diff.noncore <- diff.count(sol, comp, superset, non.core)
	return(c(
	  'my.total' = diff.all[1], 
	  'comp.total' = diff.all[2], 
	  'my.alone' = diff.all[3],
	  'comp.alone' = diff.all[4],
	  'my.core' = diff.core[1],
	  'comp.core' = diff.core[2],
	  'missing.core' = diff.core[3],
	  'extra.core' = diff.core[4],
	  'my.noncore' = diff.noncore[1],
	  'comp.noncore' = diff.noncore[2],
	  'missing.noncore' = diff.noncore[3],
	  'extra.noncore' = diff.noncore[4]
	))
}

all.candidates <- unique(c(more$solution, more$competitor))
all <- all[type %in% all.candidates]
core <- core[type %in% all.candidates]
non.core <- non.core[type %in% all.candidates]

# temporaly switch all solutions with N/Q/L/etc suffix to non-suffixed version
non.suffix <- more[grepl('\\D$', solution) & !grepl('\\D$', competitor), .(solution, competitor)][!duplicated(solution)]
to.change <- non.suffix$competitor
names(to.change) <- non.suffix$solution
to.change <- to.change[!to.change %in% more$solution]
print(to.change)
more[solution %in% names(to.change), solution := to.change[solution]]

# 3: better candidate search iterations
bad <- 1
# TODO: it might keep running, ie: A -> B -> C -> A -> ...
while(length(bad) > 0){
	solution <- more[rank == 1, solution]
	max.hit <- sum(apply(mat2[, solution], 1, max))
	comp.info <- data.table(do.call(rbind, mclapply(1:nrow(more), function(x){
		sol <- more[x, solution]
		comp <- more[x, competitor]
		get.diff(sol, comp, solution)
	})))
	print(summary(comp.info))
	bad <- which(comp.info[, missing.core * 5 + missing.noncore < 0 & ((comp.noncore > 15 & missing.noncore < -5) | (missing.core < -2))])
	bad <- bad[!more[bad, competitor] %in% solution]
	bad <- bad[!duplicated(more[bad, solution])]
	bad <- bad[which.min(comp.info[bad, missing.noncore])]
	print(cbind(more[bad], comp.info[bad]))
	better <- more[bad, competitor]
	names(better) <- more[bad, solution]
	print(better)
	more[solution %in% names(better), solution := better[solution]]
}

# 4: generate some diagnositic numbers
solution <- more[rank == 1, solution]
max.hit <- sum(apply(mat2[, solution], 1, max))
comp.info <- do.call(rbind, mclapply(1:nrow(more), function(x){
	sol <- more[x, solution]
	comp <- more[x, competitor]
	get.diff(sol, comp, solution)
}))
more <- cbind(more, comp.info)
more[, missing := missing.core + missing.noncore - extra.core - extra.noncore]

#reads <- rbind(core[type %in% solution], non.core[type %in% solution])
#het.rate <- function(SD){
#	if(nrow(SD) == 2){
#		a <- SD[1, solution]
#		b <- SD[2, solution]
#		print(a)
#		print(b)
#		exon.shared <- shared.exon(a, b)
#		print(exon.shared)
#		a.rest <- solution[solution != a]
#		b.rest <- solution[solution != b]
#		a.reads <- reads[type == a & EXON %in% exon.shared] 
#		a.reads.others <- reads[type %in% a.rest]
#		a.reads <- a.reads[!q %in% a.reads.others$q]
#		b.reads <- reads[type == b & EXON %in% exon.shared] 
#		b.reads.others <- reads[type %in% b.rest]
#		b.reads <- b.reads[!q %in% b.reads.others$q]
#		print(c(nrow(a.reads), nrow(b.reads)))
##		print(reads[q %in% a.reads$q & type %in% c(a, b)])
#		if(b == 'A*24:03'){
#			print(reads[q %in% b.reads$q & type %in% c(a, b)])
#		}
#	}else{
#	}
#}
#
#het <- more[rank == 1]
#het[, gene := sub('\\*.+', '', solution)]
#het[, het.rate(.SD), by = gene]
#het[, sum := my.alone + my.noncore.sp]
#het.ratio <- het[, .(ratio = max(sum) / min(sum), min = solution[which.min(sum)]), by = gene]
#print(het.ratio)
#het.ratio <- het.ratio[ratio > 10]
##more[solution %in% het.ratio$min, rank := 1000L + rank]

important <- function(solution){
	explained <- sum(apply(mat[, solution], 1, max))
	delta <- explained - sapply(seq_along(solution), function(i) sum(apply(mat[, solution[-i]], 1, max)))
	gene <- sub('\\*.+', '', solution)
	gene.unique <- unique(gene)
	weight <- ifelse(gene.unique %in% c('A', 'B', 'C'), 89+91, 89)
	weight <- weight / sum(weight) / 2
	names(weight) <- gene.unique
	delta / explained / weight[gene]
}
more[, importance := 0]
more[rank == 1, importance := important(solution)]
more <- more[order(rank)]

print(more[rank == 1])
print(more[rank == 2])
write.table(more, row = F, sep = '\t', quo = F, file = out.path)

#wrong <- c('C*04:01', 'C*04:09N')
#key.match <- core[type %in% c(more[rank == 1, solution], wrong)]
#noncore.match <- non.core[type %in% c(more[rank == 1, solution], wrong)]
#save(key.match, noncore.match, file = 'temp.rda')
