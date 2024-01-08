#Rscript to generate dotplots from outputs of downstream isoseq pipeline.
#load libraries
library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(glue)
#cur_nhp = snakemake.wildcards[["SPRPOP"]]
#cur_sample = snakemake.wildcards[["SMP"]]
cur_tbl_path = snakemake@input[["tbl"]]
regions_path <- snakemake@input[["rgns"]]
cur_nhp = snakemake@params[['cur_nhp']]
cur_sample = snakemake@params[['cur_sample']]

tbl = read_tsv(cur_tbl_path) %>% distinct()
colnames(tbl)[1] = "reference_name"
regions = read_table(regions_path, col_names = c("seqnames", "start", "stop", "name", ".", "strand"))


#turn ranges into GRanges
regions = makeGRangesFromDataFrame(regions, keep.extra.columns=TRUE)
regions$gene = regions$name #add gene column for logic I pulled from previous code

#figure out primary alignments by matches column
tbl = tbl %>% arrange( desc(matches) )
tbl = tbl %>% group_by(query_name) %>% mutate(alignment_n =  row_number() ) %>% mutate(is_primary = ifelse(alignment_flag < 20, yes = TRUE, no = FALSE) ) %>% ungroup() #order alignments by matches and number them as primary, secondary, etc.
tbl_ranges = makeGRangesFromDataFrame(tbl, keep.extra.columns = T, 
                                      seqnames.field = "reference_name", 
                                      start.field = "reference_start",end.field = "reference_end", strand.field = "strand")



#get overlaps
overlaps = GenomicRanges::findOverlaps(tbl_ranges, regions)
if( length(overlaps@from) > nrow(tbl) ){
    stop(glue("reads are overlapping multiple regions. Not sure which region to annotate them as. Check table and regions for {cur_sample}.") )
}
#remove reads that don't overlap with : regions.bed
tbl$gene = NA
tbl$gene[overlaps@from] = regions$gene[overlaps@to] #add gene name
tbl = tbl[overlaps@from,] #remove any genes that were not overlaps

#generate first plot
#get plot with all alignments
MIN_p_THRESH = snakemake@config[["min_perID"]]

#plot 1: percentage identity plot for all reads
p1 = ggplot(data = tbl %>% filter(perID_by_events >= MIN_p_THRESH), mapping = aes(x = gene, y = perID_by_events)) + 
  geom_jitter(width = 0.1, height = 0, mapping = aes(color = is_primary), alpha = 0.3) +
  ggtitle(glue("{cur_nhp} {cur_sample} isoseq read primary alignments to nhp reference paralogs.")) +
  xlab(glue("{cur_nhp} TBC1D3 paralogs identified by mapping")) +
  ylab("percent identity by events") +
  scale_color_discrete(name = "is primary alignment") +
  theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1.2))


#plot 2: primary vs secondary plot
#get stats of secondary vs. primary alignment
prim_vs_secondary = tbl %>% group_by(gene, query_name) %>% arrange(desc(perID_by_events)) %>% filter(row_number() == 1 ) %>% ungroup() %>%
  arrange(query_name, desc(is_primary), perID_by_events)


prim_vs_secondary = prim_vs_secondary %>% 
  group_by(query_name) %>% 
  arrange(desc(is_primary), desc(perID_by_events)) %>% 
  filter(row_number() %in% c(1,2) ) %>%
  summarize(difference_with_secondary = diff(perID_by_events) ) %>%
  ungroup()

tbl2 = left_join(tbl, prim_vs_secondary, by = "query_name" )
tbl2 = tbl2 %>% mutate(difference_with_secondary = ifelse(is.na(difference_with_secondary), 1, difference_with_secondary)) #make NAs 1
p2 = ggplot(data = tbl2 , mapping = aes(x = gene, y = difference_with_secondary)) + 
  geom_jitter(width = 0.1, height = 0,  mapping = aes(color = is_primary)) +
  ggtitle(glue("{cur_nhp} {cur_sample} isoseq: primary vs secondary alignments %-ID")) +
  xlab(glue("{cur_nhp} TBC1D3 paralogs identified by mapping")) +
  ylab("percent identity difference")


#save plots
plots_list = list(p1, p2)
pdf( snakemake@output[['plt']] , onefile = TRUE, width = 9, height = 6 )
for(i in 1:length(plots_list) ){
  plot(plots_list[[i]])
}
dev.off()
