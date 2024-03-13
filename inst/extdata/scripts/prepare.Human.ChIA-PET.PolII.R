library(InteractionSet)
library(readxl)
f <- '~/Downloads/GSE33664/mmc2.xls'
d <- lapply(c('Table S3g. K562_saturated IXN_1', 'Table S3g. K562_saturated IXN_2', 'Table S3h. MCF7_saturated IXN'), function(.ele) read_xls(f, sheet=.ele, skip=3))
d <- list('K562'=rbind(d[[1]], d[[2]]), 'MCF7'=d[[3]])
gi <- lapply(d, function(.ele){
  colnames(.ele)[1:7] <- paste0('V', 1:7)
  with(.ele, GInteractions(GRanges(V1, IRanges(V2, V3)),
                           GRanges(V4, IRanges(V5, V6)),
                           score = log2(V7)))
})
gr <- GRanges('chr8:127187814-129141059')
gi.s <- lapply(gi, function(.ele) subsetByOverlaps(.ele, gr))
gi.s
extfolder <- file.path('inst/extdata/GSE33664')
dir.create(extfolder)
mapply(gi.s, names(gi.s), FUN = function(.ele, .name){
  rtracklayer::export(.ele, file.path(extfolder, paste0(.name, '.bedpe')), format='BEDPE')
})
system('gzip inst/extdata/GSE33664/*')

