library(rtracklayer)
library(InteractionSet)
library(trackViewer)
bws <- dir("~/Downloads/kazu-bedpe/HiCAR_R2_bigwig", ".bigWig",
           full.names = TRUE)
range.ext <- range <- GRanges("chr5", IRanges(21100000, 22100000))
start(range.ext) <- start(range.ext) - 10000000
end(range.ext) <- end(range.ext) + 10000000
sigs <- lapply(bws, import.bw, selection=range)
mapply(function(.s, .n){
  export(.s, file.path('inst', 'extdata', .n), format = 'BigWig')
}, sigs, basename(bws))

interactions <- dir("~/Downloads/kazu-bedpe",
                    "10k.2.sig3Dinteractions.bedpe",
                    full.names = TRUE)
loops_loopBouquet <- lapply(interactions, function(.ele){
  .loop <- read.delim(.ele, header = FALSE)
  .loop <- with(.loop, GInteractions(
    GRanges(V1, IRanges(V2, V3)),
    GRanges(V4, IRanges(V5, V6)),
    score = V13
  ))
  .loop <- subsetByOverlaps(.loop, range)
  .loop <- .loop[seqnames(first(.loop))==seqnames(range)[1] &
                   seqnames(second(.loop))==seqnames(range)[1]]
})
mapply(function(.s, .n){
  export(.s, file.path('inst', 'extdata', .n), format = 'BEDPE')
}, loops_loopBouquet, basename(interactions))

hics <- dir("~/Downloads/kazu-bedpe/MAPS_signals",
            "10000.hic", full.names = TRUE)
gis <- lapply(hics, function(.ele){
  importGInteractions(.ele, ranges = range.ext, format = "hic",
                      resolution = 10000,
                      out = "GInteractions",
                      normalization = "NONE",
                      matrixType = "observed")
})
mapply(function(.s, .n){
  saveRDS(.s, file.path('inst', 'extdata',
                        paste0(sub('hic', 'ginteractions',
                                   sub('PEAK', 'signals', .n)),
                               '.rds')))
}, gis, basename(hics))
