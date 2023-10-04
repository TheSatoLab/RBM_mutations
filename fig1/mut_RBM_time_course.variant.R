#! usr/bin/env R

library(tidyverse)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(scales)
options(scipen=1000)

### dir="/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/"
###
# out_fig_time_course_path="data/output/${input_fasta_prefix}/histgram_mut_rate.metadata_${input_fasta_prefix}.mut_RBM.pdf"
# ${R_path} --vanilla --slave --args \
#     ${out_metadata_path} \
#     ${out_fig_time_course_path} \
#     Spike_L452 Spike_T478 Spike_E484 Spike_N501 < \
#     script/mut_RBM_time_course.R
###

metadata_mut_path <- commandArgs()[5]
out_plot_path <- commandArgs()[6]
#locus_list

# metadata_mut_path <- "/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/output/2022_03_23/metadata_2022_03_23.mut_RBM.tsv"
# locus_list <- c("Spike_L452", "Spike_T478", "Spike_E484", "Spike_N501")


metadata_mut <- read_tsv(metadata_mut_path)
metadata_mut$date[nchar(metadata_mut$date)!=10] <- NA
metadata_mut$date <- as.Date(metadata_mut$date)
first_date  <-  min(na.omit(as.numeric(metadata_mut$date)))
metadata_mut  <- metadata_mut %>% mutate(date.num = as.numeric(date) - first_date)
metadata_mut <- filter(metadata_mut , !is.na(date.num), !is.na(pango_lineage))
#metadata_mut <- select(metadata_mut, -date, -region_exposure)
#metadata_mut

# lineage_list <- metadata_mut$pango_lineage %>% unique() %>% sort()

metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="B.1.1.7"|str_detect(metadata_mut$pango_lineage, "Q\\."), "Alpha", "others")
metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="B.1.617.2"|str_detect(metadata_mut$pango_lineage, "AY\\."), "Delta",metadata_mut$pango_lineage2)
metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="B.1.1.529"|str_detect(metadata_mut$pango_lineage, "BA\\."), "Omicron", metadata_mut$pango_lineage2)
metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="P.1"|str_detect(metadata_mut$pango_lineage, "P\\.1\\."), "Gamma", metadata_mut$pango_lineage2)
metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="B.1.351"|str_detect(metadata_mut$pango_lineage, "B\\.1\\.351"), "Beta", metadata_mut$pango_lineage2)
metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="B.1.621" | metadata_mut$pango_lineage=="BB.2" |str_detect(metadata_mut$pango_lineage, "B\\.1\\.621\\."), "Mu", metadata_mut$pango_lineage2)

metadata_mut %>% group_by(pango_lineage2) %>% summarise(count=n()) %>% arrange(desc(count))
# strain_order.v <- metadata_mut %>% group_by(pango_lineage2) %>% summarise(count=n()) %>% arrange(desc(count)) %>% pull(pango_lineage2)
# strain_order.v  <- strain_order.v[-which(strain_order.v=="others" | is.na(strain_order.v))]
# strain_order.v  <- c(strain_order.v, "others")
strain_order.v <- metadata_mut$pango_lineage2 %>% unique()  %>% sort()

metadata_mut$pango_lineage2 <- factor(metadata_mut$pango_lineage2, levels = strain_order.v)

my_pal <- c(brewer.pal(length(strain_order.v)-1, "Paired"), "grey50")
g <- metadata_mut %>% ggplot() +
  geom_histogram(aes_string(x="date", fill="pango_lineage2"), binwidth = 10)+
  theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))+
  theme(legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.title =  element_text(colour="black", size=12),
        axis.line.x = element_line(color="black", size=0.5, lineend = "square"),
        axis.line.y = element_line(color="black", size=0.5, lineend = "square"),
        axis.text.y = element_text(colour="black", size=12, angle = 45),
        axis.text.x = element_text(colour="black", size=12),
        axis.ticks = element_line(color="black", size=0.5, lineend = "square"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right")+
  scale_fill_manual(values=my_pal) +
  scale_y_continuous(labels=scientific)
g
ggsave(filename=out_plot_path, plot = g, width = 12, height = 8)

