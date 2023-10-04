#!/usr/bin/env R

library(tidyverse)
library(RColorBrewer)
library(ape)
library(castor)
library(phangorn)

working_dr <-'Dropbox/Research/cov-phylo/RBM_project/'
setwd(working_dr)

tree.name <- 'data/gisaid/2022_03_23/GISAID-hCoV-19-phylogeny-2022-02-21.global.tree'
tree <- read.tree(tree.name)

res_asr_data.name <- 'data/output/2022_03_23/res_asr.global.tree_only_first_node_merged.common_200.100.mut_gain_modified.txt'
res_asr_data <- read_tsv(res_asr_data.name, col_types = 'ddc')
mut_common_list <- res_asr_data$mut %>% unique()

tree_info.name <- 'data/gisaid/2022_03_23/GISAID-hCoV-19-phylogeny-2022-02-21.global.tree.info.txt'
tree_info <- read_tsv(tree_info.name)

tree_filtered <- keep.tip(tree, filter(tree_info, is.Tip)$label)
tree_filtered$edge.length %>% sum() #5373743

mut_gain_df  <- read.table('data/output/2022_03_23/res_asr_state_table_internal_node.filtered.txt')
mut_list <- c("Spike_L452R","Spike_T478K", "Spike_E484K", "Spike_N501Y")
mut_gain_df <- select(mut_gain_df, mut_list)
rownames(mut_gain_df) <- mut_list 

all_mut_gain_node.v <- res_asr_data$node %>% unique()
length(all_mut_gain_node.v)

mut_gain_node.L452R <- filter(res_asr_data, mut=="Spike_L452R")$node
mut_gain_node.T478K <- filter(res_asr_data, mut=="Spike_T478K")$node
mut_gain_node.E484K <- filter(res_asr_data, mut=="Spike_E484K")$node
mut_gain_node.N501Y <- filter(res_asr_data, mut=="Spike_N501Y")$node
descendant.L452R <- Descendants(tree_filtered, mut_gain_node.L452R, type="all") %>% unlist() 
descendant.T478K <- Descendants(tree_filtered, mut_gain_node.T478K, type="all") %>% unlist() 
descendant.E484K <- Descendants(tree_filtered, mut_gain_node.E484K, type="all") %>% unlist() 
descendant.N501Y <- Descendants(tree_filtered, mut_gain_node.N501Y, type="all") %>% unlist() 

##only RBM
RBM_mut_gain_node.v <- res_asr_data %>% filter(mut %in% mut_list )  %>% pull(node) %>% unique()
length(RBM_mut_gain_node.v)

test.df <- tibble(node=RBM_mut_gain_node.v)
test.df$L452R_desc <- ifelse(test.df$node %in% descendant.L452R, 1, 0)
test.df$T478K_desc <- ifelse(test.df$node %in% descendant.T478K, 1, 0)
test.df$E484K_desc <- ifelse(test.df$node %in% descendant.E484K, 1, 0)
test.df$N501Y_desc <- ifelse(test.df$node %in% descendant.N501Y, 1, 0)

for(i in 1:length(mut_list)){
  mut1 <- mut_list[i]
  tmp.row <- ifelse(test.df$node %in% filter(res_asr_data, mut==mut1)$node, 1, 0)
  test.df <- test.df %>% bind_cols(tibble(tmp.row))
  colnames(test.df)[ncol(test.df)] <- mut1
}
test.df

mut_RBM.v <- colnames(test.df)[2:5]
fisher.df <- data.frame(matrix(rep(NA, 8), nrow=1))[numeric(0), ]
colnames(fisher.df) <- c("RBM", "target", "pp", "pn", "np", "nn", "odds.ratio", "p.value")
for(i in 1:4){
  for(j in 1:4){
    if(i!=j){
      mut_RBM <- mut_RBM.v[i]
      mut_target <- mut_list[j]
      mt <- table(pull(test.df, mut_RBM), pull(test.df, mut_target))
      res <- fisher.test(mt)
      odds.ratio <- res$estimate[[1]]
      p.value <- res$p.value
      tmp.row <- data.frame(mut_RBM, mut_target, mt[2,2], mt[2,1], mt[1,2], mt[1,1], odds.ratio, p.value)
      colnames(tmp.row) <- c("RBM", "target", "pp", "pn", "np", "nn", "odds.ratio", "p.value")
      fisher.df <- rbind2(fisher.df, tmp.row)
    }
  }
}

# for(i in 1:4){
#   for(j in 1:4){
#     if(i!=j){
#       mut_A <- mut_list[i]
#       #mut_RBM <- mut_RBM.v[i]
#       mut_B <- mut_list[j]
#       mut_A_desc <- str_c(gsub("Spike_", "", mut_A), "_desc")
#       test.df_tmp <- test.df[pull(test.df,mut_A)==0, ]
#       mt <- table(pull(test.df_tmp, mut_A_desc), pull(test.df_tmp, mut_B))
#       res <- fisher.test(mt)
#       odds.ratio <- res$estimate[[1]]
#       p.value <- res$p.value
#       tmp.row <- data.frame(mut_A, mut_B, mt[2,2], mt[2,1], mt[1,2], mt[1,1], odds.ratio, p.value)
#       colnames(tmp.row) <- c("RBM", "target", "pp", "pn", "np", "nn", "odds.ratio", "p.value")
#       fisher.df <- rbind2(fisher.df, tmp.row)
#     }
#   }
# }


fisher.df
fisher.df <- arrange(fisher.df, p.value)
fisher.df$fdr  <- p.adjust(fisher.df$p.value, method="fdr")
fisher.df 
write.table(fisher.df, 'data/output/2022_03_23/res_fisher_mut_gain_RBM.txt', quote = F, row.names = F, col.names = T, sep = "\t")

