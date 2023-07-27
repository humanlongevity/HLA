#!/usr/bin/env Rscript
# Author: Xie Chao

args <- commandArgs()
code.source <- sub('--file=', '', args[4]) # code.source is the name of the file
if(length(args) != 8) {
	cat("usage:", code.source, "input.tsv output.tsv cores\n")
	q()
}
set.seed(131)
data.dir <- dirname(code.source)
align.path <- args[6]
out.path <- args[7]
cores <- args[8]
print('----->> Cores in R script',cores)
print(cores)
library(parallel)
options(mc.cores = cores)
library(data.table)
library(lpSolve)

#all <- fread(sprintf("gzip -dc %s", align.path))
all <- fread(align.path)
setnames(all, c('q', 'qpos0', 'qpos', 't', 'tlen', 'ts', 'te', 'mis', 'type', 'msa', 'exon', 'specific', 'left', 'right', 'start', 'end'))
all[, specific := as.double(specific)]
all <- all[mis == 0]
print(nrow(all))

if(length(unique(all[, q])) > 20000){
	reads <- all[, .(q, type)]
	reads[, gene := sub('\\*.+', '', type)]
	reads[, q := sub('/[12]$', '', q)]
	reads <- unique(reads[, .(q, gene)])
	setkey(reads, gene)
	keep <- reads[, .(q = sample(q, min(3000, .N))), by = gene]
	all <- all[sub('/[12]$', '', q) %in% keep[, q]]
}
print(nrow(all))
print(length(unique(all[,q])))

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

# remove alleles with bad spread coverage on core exons
n.core.reads <- all[msa %in% keep & (EXON == 'E2' | (EXON == 'E3' & msa == 'ClassI')), length(unique(q))]
print(n.core.reads)
filter.by.spread <- F
if(n.core.reads >= 1000)
{
	core.cover <- all[
		msa %in% keep & (EXON == 'E2' | (EXON == 'E3' & msa == 'ClassI')), 
		.(
			tlen[1], 
			cover = sum(IRanges::width(IRanges::reduce(IRanges::IRanges(ts, te))))
		), 
		by = .(type, EXON)
	]
	core.cover[, gene := sub('\\*.+', '', type)]
	core.cover[, max.cover := max(cover/V1), by = .(gene, EXON)]
	bad.alleles <- core.cover[cover/V1 < 0.9 * max.cover[1], type, by = .(type, EXON)][, type]
	print(length(bad.alleles))
	all <- all[!type %in% bad.alleles]
	filter.by.spread <- T
}
print(nrow(all))

# only keep Class I (A, B, C) and DRB1, DQB1, and DPB1
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
cat('dcast\n')
mat <- dcast(core, q ~ type, value.var = 'specific', fun.aggregate = max, fill = 0)
cat('done dcast\n')
qs <- mat$q
mat[, q := NULL]
mat <- as.matrix(mat)
cat('weight\n')
weight <- apply(mat, 1, max)
cat('done weight\n')
mat2 <- mat
cat('set to 1\n')
mat[mat > 0 & mat < 1] <- 1
cat('done set to 1\n')

# matrix for non-core exon alignments
cat('dcast2\n')
mat.noncore <- dcast(non.core, q ~ type, value.var = 'specific', fun.aggregate = max, fill = 0)
mat.noncore[, q := NULL]
mat.noncore <- as.matrix(mat.noncore)
cat('done dcast2\n')

# filter out types with too few reads
if(!filter.by.spread)
{
	counts <- colSums(mat)
	type.counts <- data.table(type = names(counts), counts)
	type.counts[, gene := sub('\\*.+', '', type)]
#	print(type.counts[, summary(counts), by = gene])
	cand <- type.counts[, .(type = type[counts > quantile(counts, 0.3)]), by = gene][, type]
#	print(type.counts[type %in% cand, summary(counts), by = gene])
	mat <- mat[, colnames(mat) %in% cand]
}
## filter out reads with no alleles mapped to
counts <- rowSums(mat)
#summary(counts)
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
beta <- 0.002
#beta <- 0.0001

# ILP objective
f.obj <-  c(rep(0, na),   weight,      -beta * weight, 0  )
f.type <- c(rep('i', na), rep('b',nr), rep('i',   nr), 'i')

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
cat('lpsolve\n')
lps <- lp('max', f.obj, f.con, f.dir, f.rhs, int.vec = which(f.type == 'i'), binary.vec = which(f.type == 'b'))
cat('done lpsolve\n')
solution <- lps$solution[alleles]
names(solution) <- allele.names
solution <- solution[order(-solution)]
solution <- solution[solution > 0]
solution <- names(solution)
print(solution)

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
	missing.core <- diff.core[3] - diff.core[4]
	missing.noncore <- diff.noncore[3] - diff.noncore[4]
	return(c(
	  'missing' = missing.core + missing.noncore,
	  'core' = missing.core,
	  'noncore' = missing.noncore,
	  'my.total' = diff.all[1], 
	  'comp.total' = diff.all[2], 
	  'my.alone' = diff.all[3],
	  'comp.alone' = diff.all[4],
	  'my.core' = diff.core[1],
	  'comp.core' = diff.core[2],
	  'my.core.sp' = diff.core[3],
	  'comp.core.sp' = diff.core[4],
	  'my.noncore' = diff.noncore[1],
	  'comp.noncore' = diff.noncore[2],
	  'my.noncore.sp' = diff.noncore[3],
	  'comp.noncore.sp' = diff.noncore[4]
	))
}

get.comp.info <- function(sols, comps, superset){
	data.table(do.call(rbind, 
		mclapply(1:length(sols), function(x) get.diff(sols[x], comps[x], superset))
	))
}

get.better <- function(sols, comps, superset){
	comp.info <- get.comp.info(sols, comps, superset)
#	bad <- which(comp.info[, core * 5 + noncore < 0 & ((comp.noncore > 10 & noncore < -5) | (core < -2))])
	bad <- which(ifelse(grepl('[A-Z]$', comps), 
		comp.info[, noncore < -10 & core <= 0],
		comp.info[, core * 5 + noncore < 0 & (noncore < -2  | core <= -2)]))
	bad <- bad[!comps[bad] %in% superset]
#	bad <- bad[!duplicated(more[bad, solution])]
	if(length(bad) > 0){
		bad <- bad[order(comp.info[bad, missing])]
		better <- comps[bad]
		names(better) <- sols[bad]
		print(better)
		return(better)
	}else{
		return(NULL)
	}
}

# 2: round 0 for better candidate searching
max.hit <- sum(apply(mat2[, solution], 1, max))
cat("pulling non-core exons in\n")
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
	ambig <- ambig[ambig == s | (! ambig %in% solution)]
#	missing <- cand[missing <= 2, missing]
#	names(missing) <- ambig
	sol <- s
	if(length(ambig) > 1){
		bests <- sort(ambig)
		x <- as.integer(sub('.+?\\*(\\d+):.+', '\\1', bests))
		y <- as.integer(sub('.+?\\*\\d+:(\\d+).*', '\\1', bests))
		bests <- bests[order(y * 5 + x)]
		sol <- bests[1]

		history <- sol
		better <- 1
		while(length(better) >= 1){
			better <- get.better(rep(sol, length(ambig)), ambig, solution)
			if(length(better) > 0){
				sol <- better[1]
				if(sol %in% history){
					break
				}
				history <- c(history, sol)
				all <- unique(c(better, cand$competitor))
				cand <- cand[match(all, competitor)]
			}
		}
	}
	cand[, rank := 1:nrow(cand)]
	cand[, solution := sol]
	return(copy(cand))
}))

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
cat("refining solution\n")
better <- 1
# TODO: it might keep running, ie: A -> B -> C -> A -> ...
history <- more[rank == 1, solution]
while(length(better) > 0){
	solution <- more[rank == 1, solution]
	better <- get.better(more$solution, more$competitor, solution)
	if(length(better) >= 1)
	{
		better <- better[1]
		more[solution %in% names(better), solution := better[solution]]
		if(better %in% history)
		{
			break
		}
		history <- c(history, better)
	}
}

# 4: generate some diagnositic numbers
solution <- more[rank == 1, solution]
comp.info <- get.comp.info(more$solution, more$competitor, solution)
more <- cbind(more[, .(rank, solution, competitor)], comp.info)
more[, heter.reads := 0]
more[, importance := 0]

# detect ambigous calls
freq <- fread(sprintf('%s/../data/hla.freq', data.dir))
sols <- more[rank == 1, unique(solution)]
pop <- freq[allele %in% sols, .(r = median(rank)), by = population][r == min(r), population]
freq <- freq[population %in% pop][, .(rank = median(rank)), by = allele]
setkey(freq, allele)
print(freq[sols])
more[, ambig := 'NA']
ambig.n.reads <- 2
ambig.core.ratio <- 0.7
extra <- do.call(rbind, lapply(sols, function(sol){
	final <- more[rank == 1 & solution == sol]
	core.cutoff <- more[solution == sol, ambig.core.ratio * max(my.core)]
	cands <- more[
			solution == sol & 
			missing <= ambig.n.reads & 
			core <= ambig.n.reads & 
			comp.alone >= my.alone &
			comp.core >= core.cutoff, 
		competitor]
	cands <- cands[!cands %in% sols]
	if(length(cands) >= 1){
		cands <- c(sol, cands)
		ord <- freq[cands]
		ord[, field1 := as.numeric(sub('.+\\*(\\d+):(\\d+)\\D*?$', '\\1', allele))]
		ord[, field2 := as.numeric(sub('.+\\*(\\d+):(\\d+)\\D*?$', '\\2', allele))]
		setorder(ord, rank, field2, field1, na.last = T)
		print(ord)
		best <- ord[1, allele]
		if(ord[1, is.na(rank) | rank >= 100]){
			best <- sol
		}
		my.ambig = paste(ord[allele != best, allele], collapse = ';')
		if(best == sol){
			final[, rank := 0]
			final[, ambig := my.ambig]
		}else{
			final <- more[solution == sol & competitor == best]
			final <- final[, .(
				rank = 0, 
				solution = best, 
				competitor = sol, 
				missing = -1 * missing, 
				core = -1 * core, 
				noncore = -1 * noncore, 
				my.total = comp.total, 
				comp.total = my.total, 
				my.alone = comp.alone, 
				comp.alone = my.alone, 
				my.core = comp.core, 
				comp.core = my.core, 
				my.core.sp = comp.core.sp, 
				comp.core.sp = my.core.sp, 
				my.noncore = comp.noncore, 
				comp.noncore = my.noncore, 
				my.noncore.sp = comp.noncore.sp, 
				comp.noncore.sp = my.noncore.sp, 
				heter.reads = heter.reads, 
				importance = importance,
				ambig = my.ambig
			)]
		}
	}else
	{
		final[, rank := 0]
	}
	copy(final)
}))
more <- rbind(extra, more)
setorder(more, rank)

het.reads <- function(sols){
	if(length(sols) == 2){
		diff <- get.diff(sols[1], sols[2], solution)
		return(as.double(c(diff['my.alone'], diff['comp.alone'])))
	}else{
		return(as.double(1))
	}
}
het.cover <- function(sols){
	if(length(sols) == 2){
		cat(sols, fill = T)
		exon.shared <- shared.exon(sols[1], sols[2])
		a <- specific.reads(sols[1], solution, exon.shared)
		b <- specific.reads(sols[2], solution, exon.shared)
		a.cover <- 0
		b.cover <- 0
		for(e in unique(c(a$EXON, b$EXON))){
			cat(sols, e, fill = T)
			a.range <- IRanges::IRanges(a[EXON == e, ts], a[EXON == e, te])
			a.cover <- a.cover + sum(IRanges::width(IRanges::reduce(a.range)))
			b.range <- IRanges::IRanges(b[EXON == e, ts], b[EXON == e, te])
			b.cover <- b.cover + sum(IRanges::width(IRanges::reduce(b.range)))
		}
		return(as.double(c(a.cover, b.cover)))
	}else{
		return(as.double(1))
	}
}
more[rank == 0, heter.reads := het.reads(solution), by = .(sub('\\*.+', '', solution))]

important <- function(solution){
	explained <- sum(apply(mat[, solution], 1, max))
	delta <- explained - sapply(seq_along(solution), function(i) sum(apply(mat[, solution[-i]], 1, max)))
	gene <- sub('\\*.+', '', solution)
	gene.unique <- unique(gene)
	weight <- ifelse(gene.unique %in% c('A', 'B', 'C'), 89+91, 89)
	weight <- weight / sum(weight) / 2
	names(weight) <- gene.unique
	round(delta / explained / weight[gene], 3)
}
more[rank == 0, importance := important(solution)]
more <- more[order(rank)]

het <- more[rank == 0]
het[, gene := sub('\\*.+', '', solution)]
het.ratio <- het[, 
	.(
		ratio = max(heter.reads) / min(heter.reads), 
		min.imp = importance[which.min(heter.reads)],
		ratio.imp = importance[which.max(heter.reads)] / importance[which.min(heter.reads)],
		min = solution[which.min(heter.reads)]
	), by = gene]
het.ratio <- het.ratio[ratio >= 5 & (min.imp < 0.1 | ratio.imp >= 10)]
more[solution %in% het.ratio$min, rank := 1000L + rank]
setorder(more, rank)

print(more[rank == 0])
print(freq[more[rank == 0, solution]])
write.table(more, row = F, sep = '\t', quo = F, file = out.path)

