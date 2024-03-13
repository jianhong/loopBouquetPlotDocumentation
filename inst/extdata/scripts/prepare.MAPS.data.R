url <- 'https://doi.org/10.1371/journal.pcbi.1006982.s027'
tmpf <- tempfile()
download.file(url=url, destfile = tmpf)
tmpdir <- tempdir()
unzip(tmpf, exdir = tmpdir)
d <- file.path(tmpdir, 'Supplementary data')
MAPS_files <- dir(d, 'MAPS')
(names(MAPS_files) <- sub('_combine_MAPS_call', '',
                          sub('.bedpe', '', MAPS_files)))
hichipper_files <- dir(d, 'hichipper')
(names(hichipper_files) <- sub('_combine_hichipper_converted.bedpe', 
                               '', hichipper_files))
maps_loops <- lapply(MAPS_files, function(.ele){
  .ele <- read.delim(file.path(d, .ele))
  with(.ele, GInteractions(GRanges(chr1, IRanges(start1+1, end1)),
                           GRanges(chr2, IRanges(start2+1, end2)),
                           score = -log10(fdr)))
})
hichipper_loops <- lapply(hichipper_files, function(.ele){
  .ele <- read.delim(file.path(d, .ele))
  with(.ele, GInteractions(GRanges(chr1, IRanges(start1+1, end1)),
                           GRanges(chr2, IRanges(start2+1, end2)),
                           score = -log10(fdr)))
})

mm_gi <- list(
  maps_CTCF         = maps_loops[["mESC_CTCF"]],
  maps_H3K4me3      = maps_loops[["mESC_H3K4me3"]],
  hichipper_CTCF    = hichipper_loops[["mESC_CTCF"]],
  hichipper_H3K4me3 = hichipper_loops[["mESC_H3K4me3"]]
)

mm_range <- GRanges('chr10', IRanges(60600001, 61660000)) 

mm_gi <- lapply(mm_gi, function(.ele){
  .ele <- subsetByOverlaps(.ele, mm_range)
  .ele[seqnames(first(.ele))=='chr10' & seqnames(second(.ele))=='chr10']
})
mapply(mm_gi, names(mm_gi), FUN=function(.d, .n){
  rtracklayer::export(.d, file.path('inst/extdata/GSE119663', paste0(.n, '.bedpe')))
})

(hics <- dir('~/github/loopBouquetPlotDocumentation1/GSE119663',
             ".hic",
             full.names = TRUE))
## set the sample names
(names(hics) <- sub("GSE119663_F123_(.*?)_PLACseq_combine.hic", "\\1",
                    basename(hics)))
mm_range <- GRanges('chr10', IRanges(60600001, 61660000)) 
range.ext <- mm_range
seqlevelsStyle(range.ext) <- 'ensembl'
start(range.ext) <- start(range.ext) - 10000000
end(range.ext) <- end(range.ext) + 10000000
gis <- lapply(hics, function(.ele){
  .ele <- importGInteractions(.ele, ranges = range.ext,
                              format = "hic",
                              resolution = 10000,
                              out = "track",
                              normalization = "KR",
                              matrixType = "observed")
  seqlevelsStyle(.ele$dat) <- 'UCSC'
  seqlevelsStyle(.ele$dat2) <- 'UCSC'
  .ele$dat$score <- log2(.ele$dat$score)
  .ele
})
null <- mapply(gis, names(gis), FUN=function(.ele, .name){
  saveRDS(.ele, file.path('inst/extdata/GSE119663', paste0(.name, '_PLACseq.hic')))
})
