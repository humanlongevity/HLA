library(parallel)
library(IRanges)
library(data.table)
options(mc.cores = detectCores())

dt <- fread('187521910.m8')
setnames(dt, c('q', 't', 'iden', 'len', 'mis', 'gap', 'qf', 'qt', 'tf', 'tt', 'e', 'score'))
#dt <- dt[iden == 100]
dt[, gene := sub('\\*.+', '', t)]
setkey(dt, q)

nexpr <- dt[grepl('\\D-', t)]
dt <- dt[grepl('\\d-', t)]
dt[, c('LEN', 'GOOD', 'GOOD.N', 'BAD.N', 'N') := .(
    max(len), 
    paste(unique(gene[len == max(len)]), collapse = ','),
    sum(len == max(len) & iden == 100),
    sum(len < max(len) | iden < 100),
    .N)
, by = q]

good <- dt[len == LEN & iden == 100]
dt <- dt[q %in% good$q]
bad <- dt[gene == GOOD & (len < LEN | iden < 100)]
setkey(good, t)
setkey(bad, t)

score.cov <- good[, .(
    nreads = .N,
    cov = sum(width(reduce(IRanges(tf, tt))))
), by = t]
score.cov[, total := as.integer(sub('.+-', '', t))]
score.bad <- bad[, .(
    nreads.bad = .N,
    cov.bad = sum(width(reduce(IRanges(tf, tt))))
), keyby = t]
score.cov <- score.bad[score.cov]

rescue <- function(to.rescue, rescue.by, BAD = bad, GOOD = good){
    explained <- unique(GOOD[t %in% c(rescue.by, to.rescue), q])
    g <- sub('\\*.+', '', to.rescue)
    unique(BAD[gene == g & !(q %in% explained), q])
}

score.cov[, mapped := cov / total]
score.cov <- score.cov[order(-cov)]
head(score.cov)

pairs <- function(G){
    scoreA <- score.cov[grepl(sprintf('^%s', G), t)][mapped/max(mapped) >= 0.95 | cov/max(cov) >= 0.9][order(-cov)]
	print(nrow(scoreA))
    goodA <- good[t %in% scoreA$t]
    badA <- bad[gene == G & q %in% goodA$q]
    print(length(unique(badA$q)))
    todo <- combn(scoreA$t, 2, simplify = F)
    cand <- do.call(rbind, mclapply(todo, function(p){
        h1 <- rescue(p[1], p[2], badA, goodA)
        h2 <- rescue(p[2], p[1], badA, goodA)   
		reads <- length(unique(goodA[t == p[1] | t == p[2], q]))
        data.table(t1 = p[1], t2 = p[2], H1 = length(h1), H2 = length(h2), hopeless = length(unique(c(h1, h2))), reads = reads)
    }))
    setkey(scoreA, t)
    setkey(cand, t1)
    cand <- cbind(cand, scoreA[cand, .(nreads, cov, mapped, nreads.bad, cov.bad)])
    setnames(cand, c('nreads', 'cov', 'mapped', 'nreads.bad', 'cov.bad'), c('nreads1', 'cov1', 'mapped1', 'nreads.bad1', 'cov.bad1'))
    setkey(cand, t2)
    cand <- cbind(cand, scoreA[cand, .(nreads, cov, mapped, nreads.bad, cov.bad)])
    setnames(cand, c('nreads', 'cov', 'mapped', 'nreads.bad', 'cov.bad'), c('nreads2', 'cov2', 'mapped2', 'nreads.bad2', 'cov.bad2'))
    cand[, score := 2 * reads - 10 * hopeless]
    return(cand[order(-score)])
}

head(A <- pairs('A'))
head(B <- pairs('B'))
head(C <- pairs('C'))
head(DRB1 <- pairs('DRB1'))
head(DQB1 <- pairs('DQB1'))
head(DPB1 <- pairs('DPB1'))

save.image('session.rda')


