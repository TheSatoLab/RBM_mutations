#!/usr/bin/env R

library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
#sessionInfo()
#setwd('/Users/yuuka.m/Dropbox/study/ggtree/')


metadata.name <- commandArgs()[5]
mut_data.name <- commandArgs()[6]
out.name <-  commandArgs()[7]
out.metadata.name <- commandArgs()[8]
mut_list <- commandArgs()[9:length(commandArgs())]


# setwd("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project")
# metadata.name <- "data/gisaid/metadata_2021-09-07_13-20.tsv"
# mut_data.name <- "data/gisaid/mutation.aa.all_protein.long.2021-09-07_13-20.txt"
# out.name <- "data/output/mut_rate_pango.metadata_2021-09-07_13-20.modified_mut_rate.tsv"
# out.metadata.name <- "data/output/metadata_2021-09-07_13-20.mut_RBM.tsv"
# mut_list <- c("L452", "T478", "E484", "N501")


metadata <- read_tsv(metadata.name)
mut_data <- read_tsv(mut_data.name)
#metadata$virus %>% table()
#mut_data$locus %>% unique()

#metadata
#metadata_filtered <- metadata %>% filter(host=="Human") %>% select(gisaid_epi_isl, date, region_exposure, pango_lineage)
# for(i in 1:length(mut_list)){
#   mut_data_filtered <- filter(mut_data, str_detect(mut,mut_list[i]), gene=="S") %>% select(-gene)
#   colnames(mut_data_filtered)[1] <- mut_list[i]
#   metadata_filtered <- left_join(metadata_filtered, mut_data_filtered, by=c("gisaid_epi_isl"="Id"))
# }
# metadata_f_long <- gather(select(metadata_filtered,gisaid_epi_isl,pango_lineage, mut_list),
#                           key = locus,
#                           value = mut,
#                           -gisaid_epi_isl, -pango_lineage)


###select gisaid_epi_isl, pango_lineage, date###
metadata_filtered <- metadata %>% filter(Host=="Human") %>% select(3,12,4)
colnames(metadata_filtered) <- c("gisaid_epi_isl", "pango_lineage", "date")
#metadata_filtered
for(i in 1:length(mut_list)){
  mut_data_filtered <- filter(mut_data, str_detect(mut,mut_list[i]))
  colnames(mut_data_filtered)[2] <- mut_list[i]
  metadata_filtered <- left_join(metadata_filtered, mut_data_filtered, by=c("gisaid_epi_isl"="Id"))
}
# metadata_filtered$date <- as.Date(metadata_filtered$date)
# metadata_filtered %>% arrange(date) %>% head()
# metadata_filtered %>% arrange(desc(date)) %>% head()

write_tsv(metadata_filtered, out.metadata.name)

metadata_f_long <- gather(select(metadata_filtered,gisaid_epi_isl,pango_lineage, mut_list),
                          key = locus,
                          value = mut,
                          -gisaid_epi_isl, -pango_lineage)

res_df_mut_rate <- metadata_f_long %>% group_by(pango_lineage, locus) %>% mutate(total_count=n()) %>% ungroup()
res_df_mut_rate <- res_df_mut_rate %>% group_by(pango_lineage, locus, mut) %>%
  summarise(total_count=total_count, mut_rate=n()/total_count) %>% unique()
#res_df_mut_rate$locus <- gsub("[A-Z|\\*|\\-]+", "", res_df_mut_rate$locus)
res_df_mut_rate$locus <- str_extract(res_df_mut_rate$locus, pattern="[0-9]+")
write_tsv(res_df_mut_rate, out.name)





