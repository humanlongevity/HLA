#!/usr/bin/env Rscript
args <- commandArgs(T)
library(data.table)
dt <- data.frame(fread(args[1]))
dt <- dt[, seq(2, ncol(dt), by = 2)]
cn <- read.table(args[2], stringsAsFactor = F, sep = '\t')$V1
colnames(dt) <- c('Allele', cn)
dt <- dt[, grep('Allele|NMDP', colnames(dt))]
summary(dt)
subset(dt, `USA NMDP Middle Eastern or North Coast of Africa (n=70890)` > 0.1)
write.table(dt, file = args[3], row = F, sep = '\t', quo = F)
