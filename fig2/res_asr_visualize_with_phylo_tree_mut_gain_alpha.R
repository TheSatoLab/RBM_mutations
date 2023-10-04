#!/usr/bin/env R

library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggtree)
#library(ggnewscale)
library(RColorBrewer)
library(circlize)
library(ape)
library(ggnewscale)

##path##
metadata.name <- commandArgs()[5]
tree.name <- commandArgs()[6]
res_asr_df.name <- commandArgs()[7]
out.dr <- commandArgs()[8]

# input_fasta_prefix="2022_03_23"
# metadata.name <-  "/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/2022_03_23/metadata_2022_03_23.tsv"
# tree.name <- "/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/2022_03_23/GISAID-hCoV-19-phylogeny-2022-02-21.global.tree"
# res_asr_df.name <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/output/", input_fasta_prefix, "/res_asr.global.tree_only_first_node_merged.RBM.100.txt")

out.name <- str_c(out.dr, "/phylo_tree_alpha_with_state.GISAID.pdf")


##read metadata GISAID##
metadata <- read_tsv(metadata.name)
# metadata$date[nchar(metadata$date)!=10] <- NA
# metadata$date <- as.Date(metadata$date)
# first_date  <-  min(na.omit(as.numeric(metadata$date)))
# metadata  <- metadata %>% mutate(date.num = as.numeric(date) - first_date)
#colnames(metadata)
# nrow(metadata) #[1] 9574301
# metadata$`Is high coverage?` %>% table(useNA = "always")
# metadata$`Is low coverage?` %>% table(useNA = "always")
# metadata$`Is complete?` %>% table(useNA = "always")
# metadata$`N-Content` %>% hist()

##read tree##
tree <- read.tree(tree.name)
# tip.label.v <- tree$tip.label
# tip.label.v %>% length() #6295310

###read result of asr###
res_asr_df <- read_tsv(res_asr_df.name, col_types="ddlcdcdddd")
res_asr_df$state %>% table()
res_asr_df <- res_asr_df %>% select(label, num.descendants, state, Spike_L452R, Spike_T478K, Spike_E484K, Spike_N501Y)
mut_gain_node.v <- res_asr_df %>% filter(str_detect(state, "Spike_")) %>% pull(label)
length(mut_gain_node.v)


tree.merged <- full_join(tree, res_asr_df)
tree.merged@data %>% filter(state=="Spike_N501Y") %>% arrange(desc(num.descendants))
N501Y_node.v <- tree.merged@data %>% filter(state=="Spike_N501Y") %>% arrange(desc(num.descendants)) %>% pull(node)
alpha.node <-  N501Y_node.v[2] #6476473

###extract tree of alpha strain###
tree_alpha <- extract.clade(tree.merged@phylo, alpha.node)
alpha.tip.v <- tree_alpha$tip.label
alpha.tip.v %>% length()

####alpha####
metadata_alpha <- metadata %>% filter(`Pango lineage` == "B.1.1.7" | str_detect(`Pango lineage`,"^Q\\."))
metadata_alpha <- metadata_alpha %>% filter(is.na(`Is low coverage?`),
                                            Host=="Human",
                                            `Is complete?`,
                                            !(`N-Content`>0.02),
                                            str_length(`Collection date`) == 10)

metadata_alpha_filtered <- metadata_alpha %>% filter(`Accession ID` %in% alpha.tip.v)
nrow(metadata_alpha_filtered)

metadata_alpha_filtered <- metadata_alpha_filtered %>%
  select("gisaid_epi_isl"=`Accession ID`,
         "pango_lineage"=`Pango lineage`,
         "date"=`Collection date`)

#merge tree.tip.df & metadata
tree_alpha <- keep.tip(tree_alpha, metadata_alpha_filtered$gisaid_epi_isl)
#tree_alpha <- full_join(tree_alpha, res_asr_df)



# tree.tip.df <- tibble(tip_Id = 1:length(tree_filtered$tip.label), gisaid_epi_isl = tree_filtered$tip.label)
# tree.tip.df.merged <- inner_join(tree.tip.df,metadata_alpha_filtered)
# tree.tip.df.merged <- tree.tip.df.merged %>% arrange(tip_Id)


# ####delta####
# metadata_delta <- metadata %>% filter(`Pango lineage` =="B.1.1.529" | str_detect(`Pango lineage`,"AY\\."))
# metadata_delta <- metadata_delta %>% filter(is.na(`Is low coverage?`),
#                                             Host=="Human",
#                                             `Is complete?`,
#                                             !(`N-Content`>0.02),
#                                             str_length(`Collection date`) == 10)
#
# metadata_delta_filtered <- metadata_delta %>% filter(`Accession ID` %in% tip.label.v)
# nrow(metadata_delta_filtered)
#
# metadata_delta_filtered <- metadata_delta_filtered %>%
#   select("gisaid_epi_isl"=`Accession ID`,
#          "pango_lineage"=`Pango lineage`,
#          "date"=`Collection date`)
#
# #merge tree.tip.df & metadata
# tree_filtered <- keep.tip(tree, metadata_delta_filtered$gisaid_epi_isl)
# tree.tip.df <- tibble(tip_Id = 1:length(tree_filtered$tip.label), gisaid_epi_isl = tree_filtered$tip.label)
# tree.tip.df.merged <- inner_join(tree.tip.df,metadata_delta_filtered)
# tree.tip.df.merged <- tree.tip.df.merged %>% arrange(tip_Id)



# #MRCA
# delta_node <- getMRCA(tree_filtered@phylo, metadata_alpha_filtered$gisaid_epi_isl)
# tree_filtered@data %>% filter(node==delta_node)


###visualize###
res_asr_df$Spike_L452R <- ifelse(res_asr_df$Spike_L452R==1, "L452R", "-")
res_asr_df$Spike_T478K <- ifelse(res_asr_df$Spike_T478K==1,"T478K", "-")
res_asr_df$Spike_E484K <- ifelse(res_asr_df$Spike_E484K==1,"E484K", "-")
res_asr_df$Spike_N501Y <- ifelse(res_asr_df$Spike_N501Y==1,"N501Y",  "-")
res_asr_df$mut_set <- str_c(res_asr_df$Spike_L452R,
                                 res_asr_df$Spike_T478K,
                                 res_asr_df$Spike_E484K,
                                 res_asr_df$Spike_N501Y, sep="/")
mut_set_list <- c("-/-/-/-","L452R/-/-/-", "-/T478K/-/-", "-/-/E484K/-", "-/-/-/N501Y",
                  "L452R/T478K/-/-", "L452R/-/E484K/-","L452R/-/-/N501Y",
                  "-/T478K/E484K/-", "-/T478K/-/N501Y", "-/-/E484K/N501Y",
                  "L452R/T478K/E484K/-", "L452R/T478K/-/N501Y", "L452R/-/E484K/N501Y", "-/T478K/E484K/N501Y",
                  "L452R/T478K/E484K/N501Y")
res_asr_df$mut_set <- factor(res_asr_df$mut_set, levels=c(mut_set_list))
res_asr_df %>% pull(mut_set) %>% table()
head(res_asr_df)
anno_df <- res_asr_df %>% select(label, mut_set, state)


pal_list <- c("grey70", brewer.pal(4,"Pastel2"), brewer.pal(6,"Set2"), brewer.pal(4,"Dark2"), brewer.pal(9,"Greens")[9])

g <- ggtree(tree_alpha)
g$data <- g$data %>% left_join(anno_df)
g <- g + aes(color=mut_set) +
  scale_color_manual(values=pal_list) +
  new_scale_color()+
  geom_point2(aes(subset=(label %in% mut_gain_node.v), color=state), shape=16, size=8)+
  scale_color_manual(values=c(brewer.pal(5, "RdBu")[1],brewer.pal(5, "RdBu")[5]))+
  theme(legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(colour = "black", size = 12),
        legend.text =element_text(colour = "black", size = 12))
ggsave(filename=out.name,  plot = g, width = 10, height = 30)




###for test###
# tree_ex <- extract.clade(tree_alpha,  350561)
# g <- ggtree(tree_ex)
# anno_df <- tree_alpha@data %>% bind_cols(label=c(tree_alpha@phylo$tip.label, tree_alpha@phylo$node.label))  %>% select(label, mut_set)
# g$data <- left_join(g$data, anno_df)
# g$data
#
# g + aes(color=mut_set) +
#   scale_color_manual(values=pal_list)
#
# ggtree(tree_ex)$data
