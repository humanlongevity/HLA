#!/usr/bin/env Rscript
args <- commandArgs(T)
library(parallel)
options(mc.cores = detectCores())
library(data.table)

dt <- fread(args[1])
setnames(dt, c('q', 'qpos0', 'qpos', 't', 'tlen', 'ts', 'te', 'mis', 'type', 'msa', 'exon', 'specific', 'left', 'right', 'start', 'end'))
dt[, specific := as.double(specific)]
dt <- dt[mis == 0]

# for HLA alleles with frame shift variants, we require reads span over the frame shift site
frame.shift <- fread('data/hla.shift')
setnames(frame.shift, c('t', 'EXON', 'shift'))
frame.shift[, type := sub('-E.+', '', t)]
frame.shift[, MSA := sub('\\*.+', '', t)]
frame.shift[, MSA := ifelse(MSA %in% c('A', 'B', 'C'), 'ClassI', MSA)]
setkey(dt, msa, exon, ts)
frame.shift[, mapped := dt[msa == MSA & exon == EXON & ts < shift-1 & te > shift+1, .N], by = .(MSA, EXON, shift)]
frame.shift <- frame.shift[mapped > 0]
frame.shift[, t := NULL]
setkey(frame.shift, type, EXON)
frame.shift <- unique(frame.shift)
setkey(dt, type, exon)
dt <- frame.shift[dt]
spanned <- dt[ts < shift-1 & te > shift+1]
dt <- dt[type %in% spanned$type | !(type %in% frame.shift$type)]

# filter non-specific matching
#dt <- dt[specific==1 & left==0 & right==0]
dt <- dt[left==0 & right==0]
# TODO, filter core exons

# filter pair end matches
dt[, qp := sub('/\\d$', '', q)]
nr <- dt[, .(pair1 = length(unique(q))), keyby = qp]
setkey(dt, qp)
dt <- nr[dt]
nr <- dt[, .(pair2 = length(unique(q))), keyby = .(qp, type)]
setkey(dt, qp, type)
dt <- nr[dt]
dt <- dt[pair1 == pair2]
specific.pairs <- dt[pair1 == 2 & specific == 1, qp]
dt[qp %in% specific.pairs, specific := 1]

# only keep Class I (A, B, C) and DRB1, DQB1, and DPB1
# There are other genes very similar to the above:
# H, Y (27% and 11% with A), and 
# DRB3, DRB5, DRB7, DRB4, DRB2 (50%, 23%, 16%, 9%, and 4% with DRB1)
# However, including those makes performance drop. It is much easier to 
# just ignore any reads that mapped to those similar genes
keep <- c('ClassI', 'DRB1', 'DQB1', 'DPB1')
ignore <- dt[! msa %in% keep, q]
dt <- dt[msa %in% keep]
dt[q %in% ignore, specific := 0]
#dt[q %in% ignore, specific := 0.1]
#dt[specific == 0, specific := 0.1]


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
weight <- apply(mat, 1, max)
mat2 <- mat
mat[mat > 0] <- 1

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
#f.obj <-  c(rep(-gamma, na), rep(1,  nr), rep(-beta, nr), 0  )
f.obj <-  c(rep(-gamma, na), weight,      -beta * weight, 0  )
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
#f.rhs.bound[f.rhs.bound == 1 & grepl('^DRB[345]', rep(allele.genes, 2))] <- 0
#f.rhs.bound[f.rhs.bound == 2 & grepl('^DRB[345]', rep(allele.genes, 2))] <- 1

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

get.diff <- function(x, y) {
	diff.reads <- qs[which(apply(mat[, c(x, y)], 1, diff) != 0)]
	diff.match <- dt[q %in% diff.reads & type %in% c(solution, x, y)]
	by.others <- diff.match[!type %in% c(x, y)]
	#return(copy(diff.match))
	return(copy(diff.match[!q %in% by.others$q]))
}

max.hit <- sum(apply(mat2[, solution], 1, max))
empty.df <- data.frame(
	'rank' = 1,
	'solution' = '',
	'missing' = 0,
	'missing2' = 0,
	'tier1' = '',
	'tier2' = '',
	'tier3' = '',
	'best.sp' = max.hit,
	'best.nonsp' = 0,
	'comp.sp' = 0,
	'comp.nonsp' = 0,
	'competitor' = ''
)

more <- do.call(rbind, mclapply(solution, function(s){
	minus1 <- solution[solution != s]
	minus1.hit <- apply(mat2[, minus1], 1, max)
	gene <- sub('(.+?)\\*.+', '^\\1', s)
	others <- allele.names[grepl(gene, allele.names)]
	other.hit <- sapply(others, function(i) sum(pmax(minus1.hit, mat2[, i])))
	cand <- data.table('competitor' = as.character(others), 'missing' = max.hit - other.hit)
    cand <- cand[order(missing)]
	cand <- cand[1:min(50, length(others))]

	ambig <- cand[missing <= 2, competitor]
	missing <- cand[missing <= 2, missing]
	names(missing) <- ambig
	sol <- s
	if(length(ambig) > 0){
		bests <- sort(ambig)
		x <- as.integer(sub('.+?\\*(\\d+):.+', '\\1', bests))
		y <- as.integer(sub('.+?\\*\\d+:(\\d+).*', '\\1', bests))
		bests <- bests[order(x * 1e5 + y)]
		total <- sapply(bests, function(sol) sum(mat2[, sol]))
		missing <- missing[bests]
		bests <- bests[order(missing * 5 - total)]
		ambig <- bests[-1]
		sol <- bests[1]
		cand[, rank := 1:nrow(cand)]
	}

	competition <- do.call(rbind, lapply(cand$competitor, function(comp){
		diff.match <- get.diff(sol, comp)
		c(
		  	'rank' = 0,
		  	'solution' = 0,
		  	'my.total' = 0,
		  	'comp.total' = 0,
		  	'my.alone' = 0,
		  	'comp.alone' = 0,
		  	'missing' = cand[competitor == comp, missing],
		  	'missing2' = 0,
			'tier1' = 0,
		  	'tier2' = 0,
		  	'tier3' = 0,
			'best.sp' = nrow(diff.match[type == s & specific == 1]),
			'best.nonsp' = nrow(diff.match[type == s & specific <  1]),
			'comp.sp' = nrow(diff.match[type == comp & specific == 1]),
			'comp.nonsp' = nrow(diff.match[type == comp & specific <  1])
		)
	}))
	competition <- data.frame(competition)
	competition$competitor <- cand$competitor
	competition$solution <- sol
#	if(length(others) > 1){
#		competition <- subset(competition, solution != competitor)
#	}
	competition$rank <- 1:nrow(competition)
	competition$tier1 <- paste(ambig, collapse = ';')
	competition$tier2 <- paste(subset(competition, best.sp == 0 & best.nonsp > 0)$competitor, collapse = ';')
	competition$tier3 <- paste(subset(competition, best.sp > 0 & comp.sp > 0 & comp.sp * 5 >= best.sp)$competitor, collapse = ';')

	competition
}))
more <- data.table(more)
more[, solution := as.character(solution)]
more[, competitor := as.character(competitor)]
solution <- more[rank == 1, solution]
max.hit <- sum(apply(mat2[, solution], 1, max))
hit.counts <- do.call(rbind, mclapply(1:nrow(more), function(x){
	sol <- more[x, solution]
	comp <- more[x, competitor]
	minus1 <- solution[solution != sol]
	minus1.hit <- apply(mat2[, minus1], 1, max)
	other1.hit <- sum(pmax(minus1.hit, mat2[, comp]))
	missing1 <- max.hit - other1.hit
	my.total <- sum(mat2[, sol])
	comp.total <- sum(mat2[, comp])
	my.alone <- sum(mat2[minus1.hit == 0, sol])
	comp.alone <- sum(mat2[minus1.hit == 0, comp])
	gene <- sub('\\*.+', '', sol)
    minus2 <- solution[-grep(gene, solution)]
	minus2.hit <- apply(mat2[, minus2], 1, max)
	other2.hit <- sum(pmax(minus2.hit, mat2[, comp]))
	missing2 <- max.hit - other2.hit
	c(my.total, comp.total, my.alone, comp.alone, missing1, missing2)
}))
more[, my.total := hit.counts[, 1]]
more[, comp.total := hit.counts[, 2]]
more[, my.alone := hit.counts[, 3]]
more[, comp.alone := hit.counts[, 4]]
more[, missing := hit.counts[, 5]]
more[, missing2 := hit.counts[, 6]]
#    cand <- cand[order(missing * 1e8 + missing2)]

important <- function(sol){
	sol <- sub(';.+', '', as.character(sol))
	explained <- sum(apply(mat[, sol], 1, max))
	delta <- explained - sapply(seq_along(sol), function(i) sum(apply(mat[, sol[-i]], 1, max)))
	delta / explained * length(sol)
}
more[, importance := 0]
solution <- unique(sub(';.+', '', more$solution))
more[rank == 1, importance := important(solution)]
more <- more[order(rank)]
print(more[rank == 1])
write.table(more, row = F, col = F, sep = '\t', quo = F, file = args[2])

#wrong <- c('C*04:01', 'C*01:09N')
#print(get.diff(wrong[1], wrong[2]))
#key.match <- dt[type %in% c(more[rank == 1, solution], wrong)]
#save(key.match, file = 'temp.rda')
#save.image(file = sprintf('%s.temp.rda', args[2]))
