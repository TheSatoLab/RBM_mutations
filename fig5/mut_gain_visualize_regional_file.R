#!/usr/bin/env R

library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggrepel)
library(scales)
options(scipen=1000)


####input name####
# setwd("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/")
# date.info <-  "2022_03_23"
# start.date <- "2020-11-09"
# end.date <- "2021-09-04"
# region.interest <- "USA"
# metadata.out.name <- str_c("data/output/",date.info,"/transmissibility/regional/method3.spike_hap.",start.date,"_",end.date,".",region.interest,".metadata.txt",sep = "")
# res.growth_gain.name <- str_c("data/output/",date.info,"/transmissibility/student_t/method3.student_t.spike_hap.",start.date,"_",end.date,".",region.interest,".growth_gain.txt",sep = "")
# pdf.growth.gain.rank.name <- args[7]
# pdf.g.date.name <- args[8]
# pdf.g.count.name <- args[9]

args = commandArgs(trailingOnly=T)
date.info <-  args[1]
start.date <- args[2]
end.date <- args[3]
region.interest <- args[4]
metadata.out.name <- args[5]
res.growth_gain.name <- args[6]
pdf.growth.gain.rank.name <- args[7]
pdf.g.date.name <- args[8]
pdf.g.count.name <- args[9]


# mut.info.name <- str_c("data/gisaid/", date.info, "/mutation.aa.all_protein.long.", date.info, ".txt")
mut.info.name <- str_c("gisaid/", date.info, "/mutation.aa.all_protein.long.", date.info, ".txt")
# metadata.out.name <- str_c("output/", date.info, "/transmissibility/regional/method3.spike_hap.", start.date, "_", end.date, ".",region.interest,".metadata.txt")
# res.growth_gain.name <- str_c("output/", date.info, "/transmissibility/regional/method3.spike_hap.", start.date, "_", end.date, ".",region.interest,".growth_gain.95.ex_T478K_group.txt")

stat.info.growth_gain <- read_tsv(res.growth_gain.name)
metadata.analyzed <- read_tsv(metadata.out.name)
mut.info <- read_tsv(mut.info.name)


###data modify###
mut.info <- mut.info %>% filter(Id %in% metadata.analyzed$Accession.ID)
mut.info <- mut.info %>% left_join(metadata.analyzed %>% select(-hap_Id), by=c("Id"="Accession.ID"))

count.mut.df <- mut.info %>% group_by(mut) %>% summarise(count.mut=n())
count.mut.df <- count.mut.df %>% mutate(rate.mut=count.mut/nrow(metadata.analyzed))
count.mut.df %>% arrange(desc(rate.mut))
count.mut.df$feature <- 'no'

for (i in 1:nrow(count.mut.df)){
  mut <- count.mut.df$mut[i]
  if (any(str_detect(stat.info.growth_gain$feature, mut))){
    if(stat.info.growth_gain$feature[str_detect(stat.info.growth_gain$feature, mut)] %>% length==1){
      count.mut.df$feature[i] <- stat.info.growth_gain$feature[str_detect(stat.info.growth_gain$feature, mut)]
    }
  }
}
count.mut.df <- count.mut.df %>% filter(feature!= "no")
count.mut.df <- count.mut.df %>% group_by(feature) %>% mutate(count.feature=mean(count.mut))

mut.l <- c()
date.l <- c()
for(mut.target in count.mut.df$mut){
  mut.info.filtered <- filter(mut.info, mut %in% mut.target)
  mut.info.filtered <- mut.info.filtered %>% arrange(Collection.date) %>% head(1000)
  date.1000 <- mut.info.filtered$Collection.date[1000]
  mut.l <- c(mut.l, mut.target)
  date.l <- c(date.l, as.character(date.1000))
}

mut_date_df <- data.frame(mut = mut.l, date=as.Date(date.l))
count.mut.df <- count.mut.df %>% left_join(mut_date_df)
count.mut.df <- count.mut.df %>% mutate(date.mean=mean(date))

stat.info.growth_gain <- stat.info.growth_gain %>% left_join(count.mut.df %>% select(feature, count.feature, date.mean) %>% unique())

####growth_gain plot####
stat.info.growth_gain.mut <- stat.info.growth_gain %>% arrange(mean) %>%
  mutate(feature = gsub("Spike_", "S_", feature),
         #feature = ifelse(str_length(feature)>50, str_c(substr(feature, 1, 50), "...", sep=""), feature),
         feature = factor(feature, levels=feature))

hap_Id.reference.mut.v <- "S_D614G+NSP12_P323L"

#make label
mut.target.v <- c("S_L452R", "S_T478K", "S_E484K", "S_N501Y")
stat.info.growth_gain.mut <- stat.info.growth_gain.mut %>%
  mutate(label = ifelse(signif != "not signif",
                        as.character(feature),
                        ifelse(feature %in% mut.target.v,
                               as.character(feature),
                               NA)),
         lbl.wrp =  gsub(":", "\n", label))

##rank plot
ylabel <-　paste('growth gain (per ',hap_Id.reference.mut.v,')',sep="")
# g.rank <- ggplot(stat.info.growth_gain.mut %>% tail(100))
g.rank <- ggplot(stat.info.growth_gain.mut, aes(y=mean, x=reorder(feature, -mean)))
g.rank <- g.rank + geom_point(aes(color=signif))
g.rank <- g.rank + geom_label_repel(aes(label=lbl.wrp), size=2, max.overlaps = 25)
g.rank <- g.rank + geom_hline(yintercept=1, linetype="dashed",color="black")
g.rank <- g.rank + theme_set(theme_classic(base_size = 8, base_family = "Helvetica"))
g.rank <- g.rank + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       axis.title.y = element_text(colour="black", size=8),
                       axis.text.x = element_blank(),
                       axis.text.y = element_text(colour="black", size=8),
                       axis.line = element_line(color="black", size=0.5, lineend = "square"),
                       axis.ticks.y = element_line(color="black", size=0.5, lineend = "square"),
                       axis.ticks.x = element_blank(),
                       strip.text = element_text(colour="black", size=8),
                       legend.title = element_text(colour="black", size=8),
                       legend.text = element_text(colour="black", size=8))
#g.rank <- g.rank + scale_color_manual(breaks=c("positive","negative",'not signif'),values=c("lightsalmon","lightskyblue","gray70"))
g.rank <- g.rank + scale_color_manual(breaks=c("positive","negative",'not signif'),
                                      values=c(brewer.pal(3, "RdBu")[1],brewer.pal(3, "RdBu")[3],"gray70"))
g.rank <- g.rank + xlab('') + ylab(ylabel) + ggtitle(region.interest)
g.rank <- g.rank  + coord_cartesian(xlim=c(-5, nrow(stat.info.growth_gain.mut)+5))
g.rank

# pdf.growth.gain.rank.name <- str_c("output/", date.info, "/transmissibility/regional/method3.spike_hap.", start.date, "_", end.date, ".",region.interest,".growth_gain.95.ex_T478K_group.rank.pdf")
pdf(pdf.growth.gain.rank.name,width=5,height=7)
plot(g.rank)
dev.off()


####count####
ylabel <-paste('growth gain (per ',hap_Id.reference.mut.v,')',sep="")
g.count <- ggplot(stat.info.growth_gain.mut, aes(x=count.feature, y=mean)) + geom_jitter(aes(color=signif))
g.count <- g.count + geom_label_repel(aes(label=lbl.wrp), size=2, max.overlaps = 25)
g.count <- g.count + geom_hline(yintercept=1, linetype="dashed",color="gray30")
g.count <- g.count + theme_set(theme_classic(base_size = 8, base_family = "Helvetica"))
g.count <- g.count + theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.title = element_text(colour="black", size=8),
                         axis.text = element_text(colour="black", size=8),
                         axis.line = element_line(color="black", size=0.5, lineend = "square"),
                         axis.ticks = element_line(color="black", size=0.5, lineend = "square"),
                         strip.text = element_text(colour="black", size=8),
                         legend.title = element_text(colour="black", size=8),
                         legend.text = element_text(colour="black", size=8))
g.count <- g.count + scale_color_manual(breaks=c("positive","negative",'not signif'),
                                        values=c(brewer.pal(3, "RdBu")[1],brewer.pal(3, "RdBu")[3],"gray70"))
g.count <- g.count  + xlab('mut count') + ylab(ylabel) + ggtitle(region.interest)
g.count
# pdf.g.count.name <- str_c("output/", date.info, "/transmissibility/regional/method3.spike_hap.", start.date, "_", end.date, ".",region.interest,".growth_gain.95.ex_T478K_group.count.pdf")
pdf(pdf.g.count.name,width=7,height=5)
plot(g.count)
dev.off()


####date####
g.date <- ggplot(stat.info.growth_gain.mut, aes(x=date.mean, y=mean)) + geom_jitter(aes(color=signif))
# g.date <- g.date + geom_text_repel(aes(x=date.mean, y=mean, label=label), size = 3, vjust=0)
g.date <- g.date + geom_label_repel(aes(label=lbl.wrp), size=2,  max.overlaps = 25)
g.date <- g.date + theme_set(theme_classic(base_size =8, base_family = "Helvetica"))
g.date <- g.date + theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.title = element_text(colour="black", size=8),
                         axis.text = element_text(colour="black", size=8),
                         axis.line = element_line(color="black", size=0.5, lineend = "square"),
                         axis.ticks = element_line(color="black", size=0.5, lineend = "square"),
                         strip.text = element_text(colour="black", size=8),
                         legend.title = element_text(colour="black", size=8),
                         legend.text = element_text(colour="black", size=8))
g.date <- g.date + scale_color_manual(breaks=c("positive","negative",'not signif'),
                                      values=c(brewer.pal(3, "RdBu")[1],brewer.pal(3, "RdBu")[3],"gray70"))
g.date <- g.date + xlab('Date') + ylab(ylabel) + ggtitle(region.interest)
g.date
# pdf.g.date.name <- str_c("output/", date.info, "/transmissibility/regional/method3.spike_hap.", start.date, "_", end.date, ".",region.interest,".growth_gain.95.ex_T478K_group.date.pdf")
pdf(pdf.g.date.name,width=7,height=5)
plot(g.date)
dev.off()



