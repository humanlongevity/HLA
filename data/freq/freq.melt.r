library(data.table)

dt <- fread('wide/freq-A.tsv')
dt <- melt(dt, id = 'Allele')
setorder(dt, -value, na.last = T)
write.table(dt, file = 'A.tsv', row = F, sep = '\t', quo = F)

dt <- fread('wide/freq-B.tsv')
dt <- melt(dt, id = 'Allele')
setorder(dt, -value, na.last = T)
write.table(dt, file = 'B.tsv', row = F, sep = '\t', quo = F)

dt <- fread('wide/freq-C.tsv')
dt <- melt(dt, id = 'Allele')
setorder(dt, -value, na.last = T)
write.table(dt, file = 'C.tsv', row = F, sep = '\t', quo = F)

dt <- fread('wide/freq-DRB1.tsv')
dt <- melt(dt, id = 'Allele')
setorder(dt, -value, na.last = T)
write.table(dt, file = 'DRB1.tsv', row = F, sep = '\t', quo = F)

dt <- fread('wide/freq-DQB1.tsv')
dt <- melt(dt, id = 'Allele')
setorder(dt, -value, na.last = T)
write.table(dt, file = 'DQB1.tsv', row = F, sep = '\t', quo = F)
