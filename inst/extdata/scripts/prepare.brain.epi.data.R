library(InteractionSet)
library(rtracklayer)
human_loops <- read.delim('~/Downloads/BrainComparativeEpigenome/Supplementary Tables 1-34/Supplementary_Table_22.tsv')
mouse_loops <- read.delim('~/Downloads/BrainComparativeEpigenome/Supplementary Tables 1-34/Supplementary_Table_25.tsv')
hs_gi <- GInteractions(GRanges(sub('-', ':', human_loops$hg38_coord_bin1)), GRanges(sub('-', ':', human_loops$hg38_coord_bin2)), celltype=human_loops$celltype)
names(hs_gi) <- human_loops$loop_id
mm_gi <- GInteractions(GRanges(sub('-', ':', mouse_loops$bin1)), GRanges(sub('-', ':', mouse_loops$bin2)), celltype=mouse_loops$celltype)
names(mm_gi) <- mouse_loops$loop_id
mm_range <- GRanges("chr12:35216062-39306229")
hs_range <- GRanges('chr7:13552503-17745392')

hs_gi <- subsetByOverlaps(hs_gi, hs_range)
mm_gi <- subsetByOverlaps(mm_gi, mm_range)

extdata <- file.path('inst', 'extdata', 'BrainData')
write.table(human_loops[human_loops$loop_id %in% names(hs_gi), ],file = file.path(extdata, 'Supplementary_Table_22.subset.tsv'), sep='\t', quote = FALSE, row.names = FALSE)
write.table(mouse_loops[mouse_loops$loop_id %in% names(mm_gi), ],file = file.path(extdata, 'Supplementary_Table_25.subset.tsv'), sep='\t', quote = FALSE, row.names = FALSE)

hs_bws <- dir('~/Downloads/BrainComparativeEpigenome/ATAC/human_m1_atac',
              '*.bw', full.names = TRUE)
bws <- lapply(hs_bws, import.bw, selection=hs_range)
dir.create(file.path(extdata, 'human_m1_atac'))
null <- mapply(bws, basename(hs_bws), FUN=function(bw, n){
  export(bw, file.path(extdata, 'human_m1_atac', n))
})

mm_bws <- dir('~/Downloads/BrainComparativeEpigenome/ATAC/mouse_mop_atac',
              '*.bw', full.names = TRUE)
bws <- lapply(mm_bws, import.bw, selection=mm_range)
dir.create(file.path(extdata, 'mouse_mop_atac'))
null <- mapply(bws, basename(mm_bws), FUN=function(bw, n){
  export(bw, file.path(extdata, 'mouse_mop_atac', n))
})
