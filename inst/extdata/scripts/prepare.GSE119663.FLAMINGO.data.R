library(FLAMINGOrLite)
res = flamingo_main(hic_data='ignore/GSE119663/GSE119663_F123_CTCF_PLACseq_combine.hic',
                    file_format='hic',
                    domain_res=1e6,
                    frag_res=5e3,
                    chr_name='chr10',
                    normalization='KR',
                    nThread=6)
gr <- with(res, GRanges(seqnames = chr, ranges=IRanges(start+1, end),
                        x=x, y=y, z=z))
saveRDS(gr, 'ignore/GSE119663/GSE119663_F123_CTCF.chr10.flamingo.rds')
res4 = flamingo_main(hic_data='ignore/GSE119663/GSE119663_F123_H3K4me3_PLACseq_combine.hic',
                     file_format='hic',
                     domain_res=1e6,
                     frag_res=5e3,
                     chr_name='chr10',
                     normalization='KR',
                     nThread=6)
gr4 <- with(res4, GRanges(seqnames = chr, ranges=IRanges(start+1, end),
                          x=x, y=y, z=z))
saveRDS(gr4, 'ignore/GSE119663/GSE119663_F123_H3K4me3.chr10.flamingo.rds')

mm_range <- GRanges('chr10', IRanges(60600001, 61660000))
gr.s <- subsetByOverlaps(gr, mm_range)
saveRDS(gr.s, 'inst/extdata/GSE119663/CTCF.chr10.flamingo.rds')
gr4.s <- subsetByOverlaps(gr4, mm_range)
saveRDS(gr4.s, 'inst/extdata/GSE119663/H3K4me3.chr10.flamingo.rds')

