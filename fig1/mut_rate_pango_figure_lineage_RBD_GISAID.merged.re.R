#!/usr/bin/env R

library(tidyverse)
library(ggtree)
library(reshape2)
library(ggnewscale)
library(dummies)
library(phangorn)
library(ape)
library(Biostrings)
library(outliers)
library(RColorBrewer)
library(ggstance)
library(patchwork)
library(circlize)

setwd("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project")
mut_rate_df.name <- "data/output/2022_03_23/mut_rate_pango.metadata_2022_03_23.RBM.txt"
tree.name <-'data/gisaid/timetree.pango_lineage_2022_01_19.nwk'
# out_fig.name <-'data/gisaid/barplot_mut_rate.metadata_2021-09-07_13-20.RBM.with_mutation_type.phylotree_GISAID.pdf'

#read tree
tree <- read.tree(tree.name)
#ggtree(tree) + geom_tiplab()
pango_lineage.v <- tree$tip.label
length(pango_lineage.v)

####read pango_lineage mut data and modify####
mut_rate_df <- read_tsv(mut_rate_df.name)
mut_rate_df$mut_rate  <- ifelse(is.na(mut_rate_df$mut), 0, mut_rate_df$mut_rate)
pango.v <- mut_rate_df %>% pull(pango_lineage) %>% unique()
pango_lineage.v <- pango_lineage.v[pango_lineage.v  %in% pango.v]
length(pango_lineage.v)
mut_rate_df <- filter(mut_rate_df, pango_lineage %in% pango_lineage.v)
mut_rate_df <- filter(mut_rate_df, total_count>200)
mut_rate_df %>% pull(pango_lineage) %>% unique() %>% length()

pango_count_df <- mut_rate_df %>% select(pango_lineage, total_count) %>% distinct()

delta.pango.v <- mut_rate_df %>% filter(pango_lineage=="B.1.617.2" | str_detect(pango_lineage,"^AY\\.")) %>% pull(pango_lineage) %>% as.character() %>% unique()
alpha.pango.v <- mut_rate_df %>% filter(pango_lineage=="B.1.1.7" | str_detect(pango_lineage,"^Q\\.")) %>% pull(pango_lineage) %>% as.character() %>% unique()
gamma.pango.v <- mut_rate_df %>% filter(pango_lineage=="P.1" | str_detect(pango_lineage,"^P\\.1\\.")) %>% pull(pango_lineage) %>% as.character() %>% unique()
beta.pango.v <- mut_rate_df %>% filter(pango_lineage=="B.1.351" | str_detect(pango_lineage,"^B\\.1\\.351\\.")) %>% pull(pango_lineage) %>% as.character() %>% unique()
mu.pango.v <- mut_rate_df %>% filter(pango_lineage=="B.1.621" | pango_lineage=="BB.2" | str_detect(pango_lineage,"B\\.1\\.621\\.")) %>% pull(pango_lineage) %>% as.character() %>% unique()
omicron.pango.v <- mut_rate_df %>% filter(pango_lineage=="B.1.1.529" | str_detect(pango_lineage,"^BA\\.")) %>% pull(pango_lineage) %>% as.character() %>% unique()

mut_rate_df$strain <- ifelse(mut_rate_df$pango_lineage %in% delta.pango.v, "Delta", mut_rate_df$pango_lineage)
mut_rate_df$strain <- ifelse(mut_rate_df$pango_lineage %in% alpha.pango.v, "Alpha", mut_rate_df$strain)
mut_rate_df$strain <- ifelse(mut_rate_df$pango_lineage %in% gamma.pango.v, "Gamma", mut_rate_df$strain)
mut_rate_df$strain <- ifelse(mut_rate_df$pango_lineage %in% beta.pango.v, "Beta", mut_rate_df$strain)
mut_rate_df$strain <- ifelse(mut_rate_df$pango_lineage %in% mu.pango.v, "Mu", mut_rate_df$strain)
mut_rate_df$strain <- ifelse(mut_rate_df$pango_lineage %in% omicron.pango.v, "Omicron", mut_rate_df$strain)

pango_strain_df <- mut_rate_df %>% select(pango_lineage, strain) %>% distinct()
pango_strain_df <- pango_strain_df  %>% left_join(pango_count_df)
strain_count_df <- pango_strain_df %>% group_by(strain) %>% summarise(total_count.strain=sum(total_count))

mut_rate_df.merged <- mut_rate_df %>% group_by(strain,locus, mut) %>%
  summarize(mut_rate.mean=mean(mut_rate)) %>% ungroup()
mut_rate_df.merged <- mut_rate_df.merged %>% left_join(strain_count_df)

mut_major_df <- mut_rate_df.merged %>% mutate(mut_count=total_count.strain*mut_rate.mean) %>%
  group_by(locus, mut) %>% summarise(mut_total_count=sum(mut_count)) %>%
  filter(!is.na(mut)) %>% top_n(3, mut_total_count) %>% arrange(locus, desc(mut_total_count))

mut_rate_df.merged.modified <- mut_rate_df.merged %>%
  mutate(mut_p = ifelse(mut %in% mut_major_df$mut, mut, ifelse(is.na(mut), mut, "others"))) %>%
  group_by(strain, locus, mut_p ,total_count.strain) %>% summarise(mut_rate=sum(mut_rate.mean))


####Modify Tree####
tree_filtered <- keep.tip(tree, mut_rate_df$pango_lineage)
tree.tip.df <- tibble(tip_Id = 1:length(tree_filtered$tip.label), pango_lineage = tree_filtered$tip.label)

MRCA_delta <- getMRCA(tree_filtered,delta.pango.v)
MRCA_alpha <- getMRCA(tree_filtered,alpha.pango.v)
MRCA_gamma <- getMRCA(tree_filtered,gamma.pango.v)
MRCA_beta <- getMRCA(tree_filtered,beta.pango.v)
MRCA_mu <- getMRCA(tree_filtered,mu.pango.v)
MRCA_omicron <- getMRCA(tree_filtered,omicron.pango.v)


tree_joined <- full_join(tree_filtered, pango_count_df, by=c('label'='pango_lineage'))
# ggtree(tree_joined) +
#   geom_tippoint(aes(colour=total_count),size=5) +
#   scale_colour_gradient(low='blue', high='red')

####tree visualize####
my_pal <- c(brewer.pal(3,"Dark2"), "grey70")
p <- ggtree(tree_joined)
p1 <- ggtree::collapse(p, node=MRCA_delta, mode='none', clade_name="Delta")
p2 <- ggtree::collapse(p1, node=MRCA_alpha, mode='none', clade_name="Alpha")
p3 <- ggtree::collapse(p2, node=MRCA_gamma, mode='none', clade_name="Gamma")
p4 <- ggtree::collapse(p3, node=MRCA_beta, mode='none', clade_name="Beta")
p5 <- ggtree::collapse(p4, node=MRCA_mu, mode='none', clade_name="Mu")
p6 <- ggtree::collapse(p5, node=MRCA_omicron, mode='none', clade_name="Omicron")

p6 + geom_tiplab(size=2) +
  geom_cladelabel(node=MRCA_delta, "Delta") +
  geom_cladelabel(node=MRCA_alpha, "Alpha") +
  geom_cladelabel(node=MRCA_gamma, "Gamma") +
  geom_cladelabel(node=MRCA_beta, "Beta") +
  geom_cladelabel(node=MRCA_mu, "Mu")


major_strain_count_df <- strain_count_df %>%
  filter(!str_detect(strain, "[0-9]")) %>%
  filter(str_detect(strain, "[a-z]")) %>%
  rename("strain"="label")

p6$data <- p6$data %>% left_join(major_strain_count_df)
p6$data$total_count <- ifelse(is.na(p6$data$total_count), p6$data$total_count.strain ,p6$data$total_count)
p6$data$total_count_log <- log10(p6$data$total_count)

color_tip <- brewer.pal(5 , "RdBu")[1]
color_anno <- brewer.pal(5 , "RdBu")[5]
anno_size <- 7
p_tree <- p6 +
  geom_tippoint(aes(colour=total_count_log),size=4) +
  geom_nodepoint(aes(color=total_count_log),size=4) +
  scale_colour_gradient(low='grey70', high=color_tip,na.value=NA) +
  geom_point2(aes(subset=(node==MRCA_delta)),shape=22, size=anno_size , color=color_anno) +
  geom_point2(aes(subset=(node==MRCA_alpha)),shape=22, size=anno_size , color=color_anno) +
  geom_point2(aes(subset=(node==MRCA_gamma)),shape=22, size=anno_size , color=color_anno) +
  geom_point2(aes(subset=(node==MRCA_beta)),shape=22, size=anno_size , color=color_anno) +
  geom_point2(aes(subset=(node==MRCA_mu)),shape=22, size=anno_size , color=color_anno) +
  geom_point2(aes(subset=(node==MRCA_omicron)),shape=22, size=anno_size , color=color_anno) +
  geom_cladelabel(node=MRCA_delta, "Delta", fontsize = 4, align=T) +
  geom_cladelabel(node=MRCA_alpha, "Alpha", fontsize = 4, align=T) +
  geom_cladelabel(node=MRCA_gamma, "Gamma", fontsize = 4, align=T) +
  geom_cladelabel(node=MRCA_beta, "Beta", fontsize = 4, align=T) +
  geom_cladelabel(node=MRCA_mu, "Mu", fontsize = 4, align=T)  +
  geom_cladelabel(node=MRCA_omicron, "Omicron", fontsize = 4, align=T)
p_tree
####merge tree and mut rate####
strain.order.v <- p6$data %>% arrange(y) %>% filter(!is.na(total_count), !is.na(x)) %>% pull(label)
length(strain.order.v)
mut_rate_df.merged.modified$strain <- factor(mut_rate_df.merged.modified$strain, levels=strain.order.v)
locus_list <- mut_rate_df.merged.modified$locus %>% unique() %>% sort()

###bar plot for each locus###
##colorpal
pal.base <- brewer.pal(4, "Dark2")
pal1 <- colorRamp2(c(0, 5), c('white', pal.base[1]))
pal2 <- colorRamp2(c(0, 5), c('white', pal.base[2]))
pal3 <- colorRamp2(c(0, 5), c('white', pal.base[3]))
pal4 <- colorRamp2(c(0, 5), c('white', pal.base[4]))

# mut_major_df$color <- c(rev(brewer.pal(3,"Greens")), rev(brewer.pal(3,"Oranges")), rev(brewer.pal(3,"Purples")),rev(brewer.pal(3,"RdPu")))
mut_major_df$color <- c(pal1(c(5, 3, 1)), pal2(c(5, 3, 1)), pal3(c(5, 3, 1)),pal4(c(5, 3, 1)))

#L452
locus_n <- locus_list[1]
mut_rate_merged_df.tmp <- mut_rate_df.merged.modified %>% filter(locus==locus_n) %>% rename("mut_p"=str_c("mut_", locus_n))
mut_order.v <- mut_major_df %>% filter(locus==locus_n) %>% pull(mut)
mut_order.v  <- c(mut_order.v , "others")
mut_rate_merged_df.tmp[,3] <- factor(mut_rate_merged_df.tmp[[3]], levels=mut_order.v)
my_pal <- mut_major_df %>% filter(locus==locus_n) %>% pull(color) %>% c("grey50")
g_452 <-mut_rate_merged_df.tmp %>% ggplot() +
  geom_barh(aes(x=mut_rate,y=strain, fill=mut_452), stat="identity", size=3)+
  scale_fill_manual(values=my_pal,na.translate=FALSE) +
  scale_x_continuous(breaks=seq(0,1))+
  theme_classic(base_size = 12,base_family = "Helvetica") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(color="black", size=0.5, lineend = "square"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="black", size=0.5),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(colour = "black", size = 12),
        legend.text =element_text(colour = "black", size = 12),
        plot.title = element_text(hjust = 0.5))+
  labs(title="452")


#T478
locus_n <- locus_list[2]
my_pal <- mut_major_df %>% filter(locus==locus_n) %>% pull(color) %>% c("grey50")
mut_rate_merged_df.tmp <- mut_rate_df.merged.modified %>%
  filter(locus==locus_n) %>%
  rename("mut_p"=str_c("mut_", locus_n))
mut_order.v <- mut_major_df %>% filter(locus==locus_n) %>% pull(mut)
mut_order.v  <- c(mut_order.v , "others")
mut_rate_merged_df.tmp[,3] <- factor(mut_rate_merged_df.tmp[[3]], levels=mut_order.v)
mut_rate_merged_df.tmp

g_478 <-mut_rate_merged_df.tmp %>% ggplot() +
  geom_barh(aes(x=mut_rate,y=strain, fill=mut_478), stat="identity", size=3)+
  scale_fill_manual(values=my_pal,na.translate=FALSE) +
  scale_x_continuous(breaks=seq(0,1))+
  theme_classic(base_size = 12,base_family = "Helvetica") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(color="black", size=0.5, lineend = "square"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="black", size=0.5),
        axis.text = element_text(colour = "black", size = 12),
        axis.text.y = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(colour = "black", size = 12),
        legend.text =element_text(colour = "black", size = 12),
        plot.title = element_text(hjust = 0.5))+
  labs(title="478")


#484
locus_n <- locus_list[3]
my_pal <- mut_major_df %>% filter(locus==locus_n) %>% pull(color) %>% c("grey50")
mut_rate_merged_df.tmp <- mut_rate_df.merged.modified %>%
  filter(locus==locus_n) %>%
  rename("mut_p"=str_c("mut_", locus_n))
mut_order.v <- mut_major_df %>% filter(locus==locus_n) %>% pull(mut)
mut_order.v  <- c(mut_order.v , "others")
mut_rate_merged_df.tmp[,3] <- factor(mut_rate_merged_df.tmp[[3]], levels=mut_order.v)
mut_rate_merged_df.tmp
g_484 <- mut_rate_merged_df.tmp %>% ggplot() +
  geom_barh(aes(x=mut_rate,y=strain, fill=mut_484), stat="identity", size=3)+
  scale_fill_manual(values=my_pal,na.translate=FALSE) +
  scale_x_continuous(breaks=seq(0,1))+
  theme_classic(base_size = 12,base_family = "Helvetica") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(color="black", size=0.5, lineend = "square"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="black", size=0.5),
        axis.text = element_text(colour = "black", size = 12),
        axis.text.y = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(colour = "black", size = 12),
        legend.text =element_text(colour = "black", size = 12),
        plot.title = element_text(hjust = 0.5))+
  labs(title="484")


#501
locus_n <- locus_list[4]
my_pal <- mut_major_df %>% filter(locus==locus_n) %>% pull(color) %>% c("grey50")
mut_rate_merged_df.tmp <- mut_rate_df.merged.modified %>%
  filter(locus==locus_n) %>%
  rename("mut_p"=str_c("mut_", locus_n))
mut_order.v <- mut_major_df %>% filter(locus==locus_n) %>% pull(mut)
mut_order.v  <- c(mut_order.v , "others")
mut_rate_merged_df.tmp[,3] <- factor(mut_rate_merged_df.tmp[[3]], levels=mut_order.v)
mut_rate_merged_df.tmp
g_501 <- mut_rate_merged_df.tmp %>% ggplot() +
  geom_barh(aes(x=mut_rate,y=strain, fill=mut_501), stat="identity", size=3)+
  scale_fill_manual(values=my_pal,na.translate=FALSE) +
  scale_x_continuous(breaks=seq(0,1))+
  theme_classic(base_size = 12,base_family = "Helvetica") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(color="black", size=0.5, lineend = "square"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color="black", size=0.5),
        axis.text = element_text(colour = "black", size = 12),
        axis.text.y = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(colour = "black", size = 12),
        legend.text =element_text(colour = "black", size = 12),
        plot.title = element_text(hjust = 0.5))+
  labs(title="501")

p_tree_modified <- p_tree +
  coord_cartesian(clip = 'off')+
  #ylim(0,705) +
  theme_classic(base_size = 12,base_family = "Helvetica") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(colour = "black", size = 12),
        legend.text =element_text(colour = "black", size = 12),
        plot.margin = margin(0,0,0,3,unit = "pt")) +
  xlim(0, 2.4)+ labs(fill="log(count)")


g_final <- p_tree_modified + g_452 + g_478 + g_484 + g_501 +
  plot_layout(nrow=1, guides = "collect", widths = c(4.5, 2, 2, 2, 2))
g_final
out_fig.name2 <-"data/output/2022_03_23/barplot_mut_rate_pango.metadata_2022_03_23.RBM.with_mutation_type.phylotree_GISAID.filtered.merged.modified.pdf"
ggsave(out_fig.name2,g_final,width=12, height=25)

