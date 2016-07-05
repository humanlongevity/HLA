library(data.table)
data <- do.call(rbind, lapply(list.files('long', '.tsv', full = T), fread))
setnames(data, c('allele', 'population', 'freq'))
data[, rank := rank(-freq, na.last = T), by = population]
setorder(data, rank)
write.table(data, row = F, quo = F, sep = '\t', file = 'hla.freq')
