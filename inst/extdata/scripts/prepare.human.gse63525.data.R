library(InteractionSet)
gr <- GRanges('21:28208606-30257695')
GM12878 <- read.delim('~/Downloads/GSE63525/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt.gz')
GM12878_gi <- with(GM12878, GInteractions(GRanges(chr1, IRanges(x1, x2)),
                                          GRanges(chr2, IRanges(y1, y2)),
                                          score = log2(o) - log2(e_donut),
                                          id = seq.int(nrow(GM12878))))
IMR90 <- read.delim('~/Downloads/GSE63525/GSE63525_IMR90_HiCCUPS_looplist.txt.gz')
IMR90_gi <- with(IMR90, GInteractions(GRanges(chr1, IRanges(x1, x2)),
                                      GRanges(chr2, IRanges(y1, y2)),
                                      score = log2(o) - log2(e_donut),
                                      id = seq.int(nrow(IMR90))))

IMR90_gi_sub <- subsetByOverlaps(IMR90_gi, gr)
GM12878_gi_sub <- subsetByOverlaps(GM12878_gi, gr)

IMR90_sub <- IMR90[IMR90_gi_sub$id, ]
GM12878_sub <- GM12878[GM12878_gi_sub$id, ]
write.table(IMR90_sub, file.path('inst', 'extdata', 'GSE63525', 'GSE63525_IMR90_HiCCUPS_looplist.txt'), sep='\t', 
            quote=FALSE, row.names = FALSE)
write.table(GM12878_sub, file.path('inst', 'extdata', 'GSE63525', 'GSE63525_GM12878_HiCCUPS_looplist.txt'), sep='\t', 
            quote=FALSE, row.names = FALSE)
system(paste('gzip', file.path('inst', 'extdata', 'GSE63525', 'GSE63525_IMR90_HiCCUPS_looplist.txt')))
system(paste('gzip', file.path('inst', 'extdata', 'GSE63525', 'GSE63525_GM12878_HiCCUPS_looplist.txt')))
