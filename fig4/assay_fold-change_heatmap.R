#! usr/bin/env R

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

setwd('/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/assay')

assay.data.name <- '210924_assay_data.csv'
data <- read_csv(assay.data.name)
data$L452R <- factor(data$L452R, levels = c("0", "1"))
data$T478K <- factor(data$T478K, levels = c("0", "1"))
data$E484K <- factor(data$E484K, levels = c("0", "1"))
data$N501Y <- factor(data$N501Y, levels = c("0", "1"))

assay_list <- data$assay %>% unique()
mut_list <- colnames(data)[1:4]

#####heatmap###
data.h <- data
data.h$mut_set <- ifelse(data.h$L452R==1, "L452R", "-")
data.h$mut_set <- str_c(data.h$mut_set, "/", ifelse(data.h$T478K==1, "T478K", "-"))
data.h$mut_set <- str_c(data.h$mut_set, "/", ifelse(data.h$E484K==1, "E484K", "-"))
data.h$mut_set <- str_c(data.h$mut_set, "/", ifelse(data.h$N501Y==1, "N501Y", "-"))

mut_set_list <- data.h$mut_set %>% unique()
data.h$mut_set <- factor(data.h$mut_set, levels=mut_set_list)
#data.h <- select(data.h, assay, mut_set, value)
data.m <- data.h %>% group_by(assay, mut_set) %>% summarize(mean=mean(value))
data.m <- spread(data.m, key = assay, value = mean) %>% as.data.frame()
#variance <- function(x) var(x)*(length(x)-1)/length(x)

WT_data <- data.m[data.m$mut_set=="-/-/-/-",]

for (i in 2:4){
  data.m[,i] <- log2(data.m[,i]/WT_data[,i])
  data.m[,i] <- scale(data.m[,i], center=F, scale=T)
}
data.m
data.m[,4]  <- data.m[,4]* (-1)
mat <- data.m %>% select(-mut_set)  %>% as.matrix()
rownames(mat)  <-  data.m$mut_set
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
g <- Heatmap(mat, name = "fold_change", row_title = "mut_set", col=col_fun)
g
# pdf_name <- 'heatmap_assay_fold_change.pdf'
# pdf(file = pdf_name, width = 8, height=6)
# g
# dev.off()
#


####moderna assay####
moderna.data.name <- '/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/assay/moderna_NT50_VOC_fold_change.csv'
moderna.data <- read_csv(moderna.data.name)

#moderna.data$mut_set <- mut_set_list
#moderna.data <- select(moderna.data, -mut)
moderna.data  <- gather(data=moderna.data, -sample, key=mut_number, value = fold_change)
mutset_number <- moderna.data$mut_number %>% unique()
mut_set_df <- tibble(mut_number=mutset_number, mut_set=mut_set_list)

moderna.data <- moderna.data %>% left_join(mut_set_df) %>% select(-mut_number)
moderna.data_mean <- moderna.data %>% group_by(mut_set) %>% summarise(mean_fold_change= mean(fold_change))
moderna.data_mean$mean_fold_change <- log2(moderna.data_mean$mean_fold_change)
moderna.data_mean$mean_fold_change <- scale(moderna.data_mean$mean_fold_change, center=F, scale=T)
data.m$moderna <- moderna.data_mean$mean_fold_change



####pfizer assay####
pfizer.data.name <- '/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/assay/pfizer_NT50_VOC_fold_change.csv'
pfizer.data <- read_csv(pfizer.data.name)
pfizer.data  <- gather(data=pfizer.data, -sample, key=mut_number, value = fold_change)

pfizer.data <- pfizer.data %>% left_join(mut_set_df) %>% select(-mut_number)
pfizer.data_mean <- pfizer.data %>% group_by(mut_set) %>% summarise(mean_fold_change= mean(fold_change))
pfizer.data_mean$mean_fold_change <- log2(pfizer.data_mean$mean_fold_change)
pfizer.data_mean$mean_fold_change <- scale(pfizer.data_mean$mean_fold_change, center=F, scale=T)
data.m$pfizer <- pfizer.data_mean$mean_fold_change

####mean
data.m$neutralization <- data.m %>% select(moderna, pfizer) %>% apply(1, mean)
mat <- data.m %>% select(-mut_set,-moderna, -pfizer)  %>% as.matrix()
rownames(mat)  <-  data.m$mut_set
# col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
# g <- Heatmap(mat, name = "log fold-change", row_title = "mut_set", col=col_fun)
# g
# pdf_name <- 'heatmap_assay_fold_change_ver5_mean.pdf'
# pdf(file = pdf_name, width = 8, height=6)
# g
# dev.off()


###
mut_set_mat <- unique(select(data.h, -assay, -value, -mut_set)) %>% as.data.frame() %>% as.matrix()
rownames(mut_set_mat)  <-  data.m$mut_set
mat <- data.m %>% select(-mut_set,-moderna, -pfizer)  %>% as.matrix()
rownames(mat)  <-  data.m$mut_set
g <- Heatmap(mat, name = "log fold-change", row_title = "mut_set", col=col_fun)+
Heatmap(mut_set_mat[,"L452R"], name = "L452R", col = c("grey", "navyblue"), width = unit(5, "mm"), show_heatmap_legend = F)+
Heatmap(mut_set_mat[,2], name = "T478K", col = c("grey", "navyblue"), width = unit(5, "mm"), show_heatmap_legend = F)+
Heatmap(mut_set_mat[,3], name = "E484K", col = c("grey", "navyblue"), width = unit(5, "mm"), show_heatmap_legend = F)+
Heatmap(mut_set_mat[,4], name = "N501Y", col = c("grey", "navyblue"), width = unit(5, "mm"), show_heatmap_legend = F)
g
# pdf_name <- 'heatmap_assay_fold_change_ver5_mean_with_mut_heatmap.pdf'
# pdf(file = pdf_name, width = 8, height=6)
# g
# dev.off()


###single mutation####
single_mut_set.v <- c("L452R/-/-/-", "-/T478K/-/-", "-/-/E484K/-", "-/-/-/N501Y")
data.h.single <- data.h %>% filter(mut_set %in% single_mut_set.v)
data.m.single <- data.m %>% filter(mut_set %in% single_mut_set.v)
mut_names <- c("L452R", "T478K", "E484K","N501Y")
mut_set_mat <- unique(select(data.h.single, -assay, -value, -mut_set)) %>% as.data.frame() %>% as.matrix()
rownames(mut_set_mat)  <-  mut_names
col_fun = colorRamp2(c(-2, 0, 2), brewer.pal(5, "RdBu")[c(5,3,1)])
mat <- data.m.single %>% select(-mut_set,-moderna, -pfizer)  %>% as.matrix()
colnames(mat) <- c("fusion activity", "transmissibility", "ACE2-binding", "neutralization")
rownames(mat)  <- mut_names
g <- Heatmap(mat,
             name = "log fold-change",
             row_title = "mut", column_title = "assay",
             col=col_fun,
             column_names_gp = grid::gpar(fontsize = 12),
             column_title_gp = grid::gpar(fontsize = 12),
             row_names_gp = grid::gpar(fontsize = 12),
             row_title_gp = grid::gpar(fontsize = 12),
             heatmap_legend_param = list(labels_gp = gpar(font = 12),
                                         title_gp = grid::gpar(fontsize = 12))
             )
g
pdf_name <- 'heatmap_assay_fold_change_ver5_mean.single_mut.pdf'
pdf(file = pdf_name, width = 5, height=4)
g
dev.off()



####
mut_count_df.name <- '/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/2022_03_23/mut_set_df.2022_03_23.txt'
mut_count_df <- read_tsv(mut_count_df.name)
mut_count_df$mut_set <- gsub("Spike_", "", mut_count_df$mut_set)
mut_count_df$mut_set <- gsub("_", "-", mut_count_df$mut_set)
mut_count_df$`log(count)` <- log2(mut_count_df$count)


mat <- data.m %>% select(-mut_set,-moderna, -pfizer)  %>% as.matrix()
rownames(mat)  <-  data.m$mut_set
colnames(mat) <- c("fusion activity", "transmissibility", "ACE2-binding", "neutralization")

mut_set_count_df <- unique(select(data.h, -assay, -value)) %>% left_join(mut_count_df) %>%
  select(-Spike_L452, -Spike_T478, -Spike_E484, -Spike_N501, -count)
mut_set_mat <- unique(select(data.h, -assay, -value, -mut_set)) %>% as.data.frame() %>% as.matrix()
rownames(mut_set_mat)  <-  data.m$mut_set


count_mat <- mut_set_count_df %>% select(`log(count)`) %>% as.matrix()
rownames(count_mat)  <-  data.m$mut_set

col_binary <- c('0'='grey', '1'='navyblue')
ha <- rowAnnotation(L452R=mut_set_mat[,"L452R"], T478K=mut_set_mat[,"T478K"], E484K=mut_set_mat[,"E484K"], N501Y=mut_set_mat[,"N501Y"],
                   col=list(L452R=col_binary, T478K=col_binary, E484K=col_binary, N501Y=col_binary),
                   show_legend = c("L452R" = FALSE, "T478K" = FALSE, "E484K" = FALSE, "N501Y" = FALSE))

ha2  <- rowAnnotation(`log(count)`=count_mat, col=list(`log(count)`=colorRamp2(c(1, 20), c("white",brewer.pal(3,"Dark2")[1]))),
                    annotation_legend_param=list(`log(count)`=list(labels_gp = gpar(font = 12),title_gp = grid::gpar(fontsize = 12))))

g <- Heatmap(mat, name = "log fold-change",
             right_annotation = c(ha, ha2),
             row_title = "mut", column_title = "assay",
             col=col_fun,
             column_names_gp = grid::gpar(fontsize = 12),
             column_title_gp = grid::gpar(fontsize = 12),
             row_names_gp = grid::gpar(fontsize = 12),
             row_title_gp = grid::gpar(fontsize = 12),
             heatmap_legend_param = list(labels_gp = gpar(font = 12),
                                         title_gp = grid::gpar(fontsize = 12)))
g


pdf_name <- 'heatmap_assay_fold_change_ver5_mean_with_mut_heatmap_count_log.pdf'
pdf(file = pdf_name, width = 8, height=6)
g
dev.off()

