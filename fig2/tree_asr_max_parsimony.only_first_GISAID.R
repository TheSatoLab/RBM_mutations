#! usr/bin/env R

library(tidyverse)
library(ape)
library(castor)
library(igraph)
library(phangorn)


###############path###############
tree.name <- commandArgs()[5]
tree.info.name <- commandArgs()[6]
mut_mat.name <- commandArgs()[7]
output.name <- commandArgs()[8]
min.num.descendants <- commandArgs()[9]
#min.num.descendants <- 100

#####tenporary#####
#setwd("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project")
#working_dr_path <- 'data/GISAID-hCoV-19-phylogeny-2021-12-13'
#setwd(working_dr_path)

##input##
# tree.name <- 'data/gisaid/2022_03_23/GISAID-hCoV-19-phylogeny-2022-02-21.global.tree'
# tree.info.name <- 'data/gisaid/2022_03_23/GISAID-hCoV-19-phylogeny-2022-02-21.global.tree.info.txt'
## metadata.name <-'data/gisaid/2022_03_23/metadata_2022_03_23.tsv'
## metadata <-read_tsv(metadata.name)
## metadata_filtered <- metadata %>% filter(`Accession ID` %in% tree.tip.df$label)
# tree.name <- 'data/gisaid/2021_12_22/GISAID-hCoV-19-phylogeny-2021-12-13.global.tree'
# tree.info.name <- 'data/gisaid/2021_12_22/GISAID-hCoV-19-phylogeny-2021-12-13.info.txt'
# mut_mat.name <- 'data/gisaid/2021_12_22/mutation.df.2021_12_22.RBM.txt'

##output##
#output.name  <- 'output/res_asr.global.tree_only_first_node_merged.100.txt'


###############command###############
### read data ###
##tree info
tree.info.df <- read_tsv(tree.info.name)
tree.tip.df <- tree.info.df %>% filter(is.Tip)

##read tree
tree <- read.tree(tree.name)
tree_filtered <- keep.tip(tree, tree.tip.df$label)
# N_node <- tree_filtered$Nnode
N_tip <- length(tree_filtered$tip.label)

##mut data
mut_data <- read_tsv(mut_mat.name)
mut.v <- colnames(mut_data)[2:ncol(mut_data)]

mut.mat <- mut_data %>% data.frame(row.names = 1)
mut.mat <- mut.mat[tree.tip.df$label,]


state.list <- list()
node.0to1.list <- list()
state.array <- numeric(nrow(mut.mat)+tree_filtered$Nnode)
#ncol(mut.mat)
for(i in 1:ncol(mut.mat)){
  mut.name <- colnames(mut.mat)[i]
  state.v <- mut.mat[,i]
  #fit.Mk <- asr_mk_model(tree_filtered, state.v + 1, rate_model = "ARD")
  fit.Mk <- asr_max_parsimony(tree_filtered, state.v + 1, transition_costs="all_equal")

  state.mat <- fit.Mk$ancestral_likelihoods
  state.node.v <- state.mat[,1]
  state.node.v <- ifelse(state.node.v>0.5,0,1)
  state.list[[mut.name]] <- c(state.v,state.node.v)
  state.array <- state.array + c(state.v,state.node.v)
  state.color <- c(state.v,state.node.v)
  state.color <- factor(state.color,levels=c(0,1))
  tree.info.df.interest <- tree.info.df %>% select(parent,node,num.descendants) %>% mutate(state = state.color)
  tree.info.df.interest.parent <- tree.info.df.interest %>% select(node,state) %>% rename(parent = node, state.parent = state)
  #tree.info.df.interest.merged <- merge(tree.info.df.interest,tree.info.df.interest.parent,by="parent")
  tree.info.df.interest.merged <- left_join(tree.info.df.interest,tree.info.df.interest.parent)
  tree.info.df.interest.merged.0to1 <- tree.info.df.interest.merged %>% filter(state == 1, state.parent == 0, num.descendants >= min.num.descendants)
  node.0to1.list[mut.name] <- NA
  if(nrow(tree.info.df.interest.merged.0to1) > 0){
    node.0to1.list[[mut.name]] <- tree.info.df.interest.merged.0to1$node
  }
  node.0to1.v <- node.0to1.list[[mut.name]]
  if(!is.na(node.0to1.list[mut.name])){
    if(length(node.0to1.v)>1){
      ancestral.nodes.k.l <- list()
      for(j in 1:length(node.0to1.v)){
        node.0to1 <- node.0to1.v[j]
        ancestral.nodes.v <- Ancestors(tree_filtered, node.0to1, type = "parent")
        #ancestral.nodes.v <- na.omit(ancestral.nodes.5.l[[node.0to1]]) %>% as.vector()
        #ancestral.nodes.v <- ancestral.nodes.v[1:min(2, length(ancestral.nodes.v))]
        ancestral.nodes.k.l[[j]] <- ancestral.nodes.v
      }
      link.df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
      for(j1 in 1:length(ancestral.nodes.k.l)){
        for(j2 in 2:length(ancestral.nodes.k.l)){
          if(j1 < j2){
            overlap.v <- intersect(ancestral.nodes.k.l[[j1]],ancestral.nodes.k.l[[j2]])
            if(length(overlap.v)>0){
              temp.df <- data.frame(j1 = j1, j2 = j2)
              link.df <- rbind(link.df,temp.df)
            }
          }
        }
      }
      if(nrow(link.df)>0){
        g <- graph.data.frame(link.df,directed=F)
        dg.l <- decompose(g, min.vertices = 1)
        network_node.v.new <- c()
        for(j in 1:length(dg.l)){
          dg <- dg.l[[j]]
          dg.info <- as.data.frame(degree(dg))
          network_node.v <- as.numeric(rownames(dg.info))
          node.0to1.v.clustered <- node.0to1.v[network_node.v]
          node.0to1.new <- getMRCA(tree_filtered, node.0to1.v.clustered)
          network_node.v.new <- c(network_node.v.new,node.0to1.new)
        }
        node.0to1.v <- node.0to1.v[-(network_node.v)]
        node.0to1.v <- c(node.0to1.v, network_node.v.new)
        node.0to1.list[[mut.name]] <- node.0to1.v
      }
    }
    node.0to1.v <- node.0to1.list[[mut.name]]
    for(j in 1:length(node.0to1.v)){
      #node.0to1 <- node.0to1.v[j] %>% as.character()
      node.0to1 <- node.0to1.v[j]
      ancestral.nodes.v <- Ancestors(tree_filtered, node.0to1, type = "all")
      #ancestral.nodes.v <- na.omit(ancestral.nodes.5.l[[node.0to1]]) %>% as.vector()
      #tree.info.df.interest.ancestors_1 <- tree.info.df.interest %>% filter(node %in% ancestral.nodes.v, state==1)
      tree.info.df.interest.ancestors_1 <- intersect(node.0to1.v , ancestral.nodes.v)
      if(length(tree.info.df.interest.ancestors_1) > 0){
        node.0to1.v[j] <- NA
      }
    }
    node.0to1.v <- as.vector(na.omit(node.0to1.v))
    node.0to1.list[[mut.name]]<- node.0to1.v
  }
}
node.0to1.list

mut.0to1.list <- list()
for(i  in 1:length(node.0to1.list)){
  mut.name <- names(node.0to1.list)[i]
  node_list <- node.0to1.list[[i]]
  for(node in node_list){
    mut.0to1.list[[as.character(node)]] <- c(mut.0to1.list[[as.character(node)]], mut.name)
  }
}
for(idx in names(mut.0to1.list)){
  idx_i = as.numeric(idx)
  state.array[idx_i] <- str_c(mut.0to1.list[[idx]],collapse = "+")
}
tree.info.df$state <- state.array
for(mut.name in names(state.list)){
  tree.info.df[,mut.name] <- state.list[[mut.name]]
}
tree.info.df
tree.info.df$state %>% table()
write_tsv(x=tree.info.df, file=output.name)



