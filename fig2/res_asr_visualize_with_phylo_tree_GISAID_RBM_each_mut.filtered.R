library(tidyverse)
library(ggtree)
library(ape)
library(phangorn)
library(circlize)

################path################
tree.name <- commandArgs()[5]
res_asr_data.name <- commandArgs()[6]
out.dr <- commandArgs()[7]
mut_target.v <- c("Spike_L452R", "Spike_T478K", "Spike_E484K", "Spike_N501Y")


#tree.name <- "/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/2022_03_23/GISAID-hCoV-19-phylogeny-2022-02-21.global.tree"
#res_asr_data.name <- "/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/output/2022_03_23/res_asr.global.tree_only_first_node_merged.RBM.100.txt"
#out.dr <- "/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/output/2022_03_23/phylotree"

# res_asr_data <- read_tsv(res_asr_data.name, col_types="ddlcdcdddd")
# res_asr_data$state %>% table()
# res_asr_tip <- filter(res_asr_data, is.Tip)
#
# tree <- read.tree(tree.name)
# tree_filtered <- keep.tip(tree, res_asr_tip$label)
#
# res_asr.mut_gain <- res_asr_data %>% filter(str_detect(state, "^Spike_"))
# gain.node.v <- res_asr.mut_gain %>% pull(node)
# print(gain.node.v)
# gain.node.desc.list <- Descendants(tree_filtered, gain.node.v, type = "tips")
# #head(gain.node.desc.list)
# gain.node.desc.all <- gain.node.desc.list %>% unlist() %>% unique()
# length(gain.node.desc.all)
#
# for (i in 1:length(gain.node.v)){
#   tmp.v <- gain.node.desc.list[[i]]
#   tmp.n <- length(tmp.v)
#   rand.v <- sample(tmp.v, tmp.n %/% 10 + 1)
#   gain.node.desc.list[[i]] <- rand.v
# }
#
# tip.selected <-  gain.node.desc.list %>% unlist() %>% unique()
# length(tip.selected)
# tip.others <- setdiff(res_asr_tip$node,tip.selected)
# tip.others.selected <- sample(tip.others, size=length(tip.others) %/% 10)
# keep.tip.label.v <- res_asr_data %>% filter(node %in% union(tip.selected, tip.others.selected)) %>% pull(label)
# length(keep.tip.label.v)
# #keep.tip.label.v %>% head()
#
# tree_filtered_selected <- keep.tip(tree_filtered, keep.tip.label.v)
# keep.node.label.v <- c(tree_filtered_selected$node.label, tree_filtered_selected$tip.label)
# res_asr_data_filtered <- res_asr_data %>% filter(label %in% keep.node.label.v)

# print(length(c(tip.selected, tip.others.selected)))
# print(tree_filtered_selected$Nnode)
# print(nrow(res_asr_data_filtered))
# length(c(tip.selected, tip.others.selected))+tree_filtered_selected$Nnode == nrow(res_asr_data_filtered)

out.tree.name <- str_c(out.dr, "/GISAID-hCoV-19-phylogeny-2022-02-21.global.tree.filtered")
out.data.name <- str_c(out.dr, "/res_asr.global.tree_only_first_node_merged.RBM.100.filtered.txt")
# write_tsv(x=res_asr_data_filtered, file=out.data.name)
# write.tree(tree_filtered_selected, file=out.tree.name)
res_asr_data_filtered <- read_tsv(out.data.name, col_types="ddlcdcdddd")
tree_filtered_selected <- read.tree(out.tree.name)


###plot###
res_asr_data_filtered$Spike_L452R <- factor(res_asr_data_filtered$Spike_L452R)
res_asr_data_filtered$Spike_T478K <- factor(res_asr_data_filtered$Spike_T478K)
res_asr_data_filtered$Spike_E484K <- factor(res_asr_data_filtered$Spike_E484K)
res_asr_data_filtered$Spike_N501Y <- factor(res_asr_data_filtered$Spike_N501Y)
tree_filtered_selected.merged <- full_join(tree_filtered_selected,
                                           select(res_asr_data_filtered,label, state, Spike_L452R, Spike_T478K, Spike_E484K, Spike_N501Y),
                                           by='label')

pal_tree <- colorRamp2(c(0,5), c("white", brewer.pal(3,"Dark2")[2]))
pallete.v <- c("grey70", pal_tree(3))
for (mut_target in mut_target.v){
  gain.node.v <- tree_filtered_selected.merged@data  %>% filter(str_detect(state, mut_target)) %>% pull(node) %>% unique()
  g <- ggtree(tree_filtered_selected.merged) +
    aes_(color=as.name(mut_target)) +
    scale_color_manual(values = pallete.v) +
    geom_point2(aes(subset=(node %in% gain.node.v)), shape=16, size=10, color =pal_tree(5))+
    theme(legend.key.size = unit(0.5, 'cm'),
          legend.title = element_text(colour = "black", size = 12),
          legend.text =element_text(colour = "black", size = 12)) +
    coord_cartesian(clip = 'off')

  out.name <- str_c(out.dr, "/phylo_tree_all_with_state.GISAID.", mut_target, ".re.pdf")
  ggsave(filename=out.name,  plot = g, width = 10, height = 40)
}

