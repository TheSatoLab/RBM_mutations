#! usr/bin/env R

library(tidyverse)
library(ape)

#####path#####
tree.name <- commandArgs()[5]
tree_data.name <- commandArgs()[6]
metadata.name <- commandArgs()[7]
out.name <- commandArgs()[8]

# working_dr_path <- 'data/GISAID-hCoV-19-phylogeny-2021-12-13'
# setwd("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project")
# setwd(working_dr_path)
# 
# tree.name <- 'global.tree'
# tree_data.name <- 'metadata.csv'
# metadata.name <- 'metadata_tsv_2021_12_22/metadata.tsv'
# #output
# out.name <- 'tree_info.with_num_children.global.tree.txt'

#####command#####
### read data ###
tree_data <- read_csv(tree_data.name)

metadata <- read_tsv(metadata.name)
metadata_label.v <- metadata$`Accession ID`

##read tree and filter##
tree <- read.tree(tree.name)
tree.tip.df <- data.frame(tip_Id = 1:length(tree$tip.label), gisaid_epi_isl = tree$tip.label)
drop_id_list <- setdiff(tree.tip.df$gisaid_epi_isl, metadata_label.v)
tree.tip.df <- filter(tree.tip.df, gisaid_epi_isl %in% metadata_label.v)
tree_filtered <- drop.tip(tree, drop_id_list)

##make tree_info_df##
edge_df <- tree_filtered$edge  %>% data.frame() %>% tibble()
edge_df <- edge_df %>% arrange(X2)
N_tip <- length(tree_filtered$tip.label)
tree.info.df <- tibble(node = edge_df$X2, parent= edge_df$X1)
new_row <- c(N_tip+1, NA)
tree.info.df <- rbind(tree.info.df , new_row)
tree.info.df$is.Tip <- ifelse(tree.info.df$node <= N_tip, TRUE, FALSE)
tree.info.df <- arrange(tree.info.df, node)
tree.info.df$label <- c(tree_filtered$tip.label, tree_filtered$node.label)

tree.info.df$num.descendants <- ifelse(tree.info.df$is.Tip, 1, NA)
index_list <- nrow(tree.info.df):(N_tip+1)

for(i in index_list){
  Node <- tree.info.df$node[i]
  total_decendants <- filter(tree.info.df, parent==Node)$num.descendants %>% sum()
  tree.info.df$num.descendants[i] <- total_decendants
}
#head(tree.info.df)
write_tsv(x=tree.info.df, file=out.name)