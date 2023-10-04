#! usr/bin/env R

library(tidyverse)
library(ggtree)
library(ape)
library(castor)
library(igraph)
library(phangorn)


#####bash script#####
##make mutaion matrix##
##common mutation##
# ${R_path} --vanilla --slave --args \
#     ${mut_data_path} \
#     ${mut_mat_common_path} < \
#     script/make_mut_matrix_common.R

# ##common mutation##
# ${R_path} --vanilla --slave --args \
#     ${tree_path} \
#     ${tree_info_path} \
#     ${mut_mat_common_path} \
#     ${res_asr_common_mut_path}  < \
#     script/tree_asr_max_parsimony.only_first_GISAID.R
######################

#####path#####
mut_data.name <- commandArgs()[5]
tree.name <-  commandArgs()[6]
tree.info.name <- commandArgs()[7]
out.name <- commandArgs()[8]
min.num.descendants <- commandArgs()[9]
#min.num.descendants <- 100


#####command#####
#######make mutaion matrix############
mut_data <- read_tsv(mut_data.name)
mut_data_label.v <- mut_data$Id %>% unique()
mut_target_df <- tibble(Id=mut_data_label.v)
mut_data_top200 <- mut_data %>% group_by(mut) %>% summarize(count=n()) %>% top_n(200, count)
mut_top200.v <- mut_data_top200$mut

for (mut_target in mut_top200.v){
    mut_hold_label.v <- filter(mut_data, mut==mut_target)$Id
    mut_target_df[,mut_target] <- ifelse(mut_target_df$Id %in% mut_hold_label.v, 1, 0)
}
mut_target_df


#################asr#################

tree.info.df <- read_tsv(tree.info.name)
tree.tip.df <- tree.info.df %>% filter(is.Tip)

##read tree
tree <- read.tree(tree.name)
tree_filtered <- keep.tip(tree, tree.tip.df$label)
N_tip <- length(tree_filtered$tip.label)

##mut data
mut.v <- colnames(mut_target_df)[2:ncol(mut_target_df)]
mut.mat <- mut_target_df %>% data.frame(row.names = 1)
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

mut_gain_df <- filter(tree.info.df,str_detect(state, "\\_"))
mut_gain_df <- select(mut_gain_df, 1,5,6) %>% arrange(state)
write_tsv(x=mut_gain_df, file=out.name)


