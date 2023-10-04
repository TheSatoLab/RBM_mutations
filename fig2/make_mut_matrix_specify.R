#! usr/bin/env R

library(tidyverse)

#####path#####
mut_data.name <- commandArgs()[5]
out.name <- commandArgs()[6]
mut_list <- commandArgs()[7:length(commandArgs())]

# mut_data.name <- '/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/2021_12_22/mutation.aa.all_protein.long.2021_12_22.txt'
# out.name <- '/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/2021_12_22/mutation.df.2021_12_22.RBM.txt'
# mut_list <- c("Spike_L452R", "Spike_T478K", "Spike_E484K", "Spike_N501Y")

mut_data <- read_tsv(mut_data.name)
mut_data_label.v <- mut_data$Id %>% unique()
mut_target_df <- tibble(Id=mut_data_label.v)

for (mut_target in mut_list){
    #mut_hold_label.v <- filter(mut_data, mut==mut_target)$Id
    mut_hold_label.v <- filter(mut_data, str_detect(mut, pattern=mut_target)) %>% pull(Id)
    mut_target_df[,mut_target] <- ifelse(mut_target_df$Id %in% mut_hold_label.v, 1, 0)
}
mut_target_df
write_tsv(mut_target_df, out.name)
