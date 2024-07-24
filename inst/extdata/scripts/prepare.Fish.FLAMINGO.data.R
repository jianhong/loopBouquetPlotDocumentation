library(FLAMINGOrLite)
res = flamingo_main(hic_data='/Users/ouj/Downloads/kazu-bedpe/danRer11/cooler/mcool/0dpa5000.mcool',
                    file_format='mcool',
                    domain_res=1e6,
                    frag_res=5e3,
                    chr_name='chr5',
                    normalization='KR',
                    nThread=6)
gr <- with(res, GRanges(seqnames = chr, ranges=IRanges(start+1, end),
                        x=x, y=y, z=z))
saveRDS(gr, 'danRer11.fishfin.0dpa.chr5.flamingo.rds')
res4 = flamingo_main(hic_data='/Users/ouj/Downloads/kazu-bedpe/danRer11/cooler/mcool/4dpa5000.mcool',
                    file_format='mcool',
                    domain_res=1e6,
                    frag_res=5e3,
                    chr_name='chr5',
                    normalization='KR',
                    nThread=6)
gr4 <- with(res4, GRanges(seqnames = chr, ranges=IRanges(start+1, end),
                        x=x, y=y, z=z))
saveRDS(gr4, 'danRer11.fishfin.4dpa.chr5.flamingo.rds')
range <- GRanges("chr5", IRanges(21100000, 22100000))

## not work, got error.
# res = flamingo_main(hic_data='/Users/ouj/Downloads/kazu-bedpe/MAPS_signals/MAPS_PEAK_0dpa5000.mcool',
#                     file_format='mcool',
#                     domain_res=1e6,
#                     frag_res=5e3,
#                     chr_name='chr5',
#                     normalization='NONE',
#                     nThread=6)
# gr <- with(res, GRanges(seqnames = chr, ranges=IRanges(start+1, end),
#                         x=x, y=y, z=z))
# saveRDS(gr, 'danRer11.fishfin.maps.0dpa.chr5.flamingo.rds')
# res4 = flamingo_main(hic_data='/Users/ouj/Downloads/kazu-bedpe/MAPS_signals/MAPS_PEAK_0dpa5000.mcool',
#                      file_format='mcool',
#                      domain_res=1e6,
#                      frag_res=5e3,
#                      chr_name='chr5',
#                      normalization='NONE',
#                      nThread=6)
# gr4 <- with(res4, GRanges(seqnames = chr, ranges=IRanges(start+1, end),
#                           x=x, y=y, z=z))
# saveRDS(gr4, 'danRer11.fishfin.maps.4dpa.chr5.flamingo.rds')
# range <- GRanges("chr5", IRanges(21100000, 22100000))