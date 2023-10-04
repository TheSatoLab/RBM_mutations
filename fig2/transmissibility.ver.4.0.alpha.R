#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(ggplot2)
library(rbin)
library(cmdstanr)
library(patchwork)
library(RColorBrewer)

download.date <- "2022_03_23"
stan_f.name <- "/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/script/multinomial_independent.stan"
metadata.name <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/", download.date,  "/metadata_", download.date, ".tsv")
mut.info.name <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/", download.date, "/mutation.aa.all_protein.long.", download.date, ".txt")
mut.interest <- "Spike_L452"
mut.interest <- "Spike_E484"
out.prefix <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/output/", download.date, "/transmissibility")

download.date <- gsub( "_", "-", download.date)
download.date <- as.Date(download.date)

##########parameters##########
#general
core.num <- 1

#min numbers
limit.count.analyzed <- 50

#Transmissibility
bin.size <- 7
generation_time <- 5.0
#generation_time <- 2.1

country.interest.v <- c('USA','United Kingdom')

#model
#set.seed(1)
#cmdstanpath <- "/usr/local/cmdstan"
#set_cmdstan_path("/Users/yuuka.m/.cmdstan/cmdstan-2.29.2")
multi_nomial_model <- cmdstan_model(stan_f.name)

metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)
# colnames(metadata)

metadata <- metadata %>% distinct(Accession.ID,.keep_all=T) %>%
  filter(Host == "Human",
         str_length(Collection.date) == 10,
         !N.Content >　0.02,
         Pango.lineage != "",
         Pango.lineage != "None",
         str_detect(Additional.location.information,"[Qq]uarantine", negate=T)) %>%
  select(Accession.ID, Collection.date, Pango.lineage, Location)

metadata <- metadata %>%
  mutate(Collection.date = as.Date(Collection.date),
         region = str_split(Location," / ",simplify = T)[,1],
         country = str_split(Location," / ",simplify = T)[,2]) %>%
  select(-Location)

metadata <- metadata %>%
  mutate(region_analyzed = ifelse(country %in% country.interest.v,
                                  as.character(country),
                                  as.character(region))) %>%
  select(-country, -region)

###extract alpha###
metadata <- metadata %>% mutate(Pango.lineage = ifelse(Pango.lineage == "B.1.1.7" | str_detect(Pango.lineage,"Q\\."),"Alpha",as.character(Pango.lineage)))
metadata <- metadata %>% filter(Pango.lineage == "Alpha")

###load mut info and extract target loci###
mut.info <- read_tsv(mut.info.name)
mut.info <- mut.info %>% filter(Id %in% metadata$Accession.ID)
alpha_N501Y_Id.v <- mut.info %>% filter(mut=="Spike_N501Y") %>% pull(Id)
mut.info <- mut.info %>% filter(str_detect(mut,mut.interest))

metadata <- metadata %>% filter(Accession.ID %in% alpha_N501Y_Id.v)


########extract data of target duration########
#period to be analyzed
#day.delay <- 150
#day.analyzed <- 300
# date.end <- download.date - day.delay
# date.start <- date.end - day.analyzed + 1
# metadata %>% arrange(Collection.date) %>% head(20)
date.end <- metadata %>% arrange(Collection.date) %>% pull(Collection.date) %>% tail(1)
date.start <- metadata %>% arrange(Collection.date) %>% pull(Collection.date) %>% head(1)

metadata.analyzed <- metadata %>% filter(Collection.date >= date.start, Collection.date <= date.end)
Id.analyzed.v <- as.character(metadata$Accession.ID)

metadata.analyzed <- metadata %>% filter(Collection.date >= date.start, Collection.date <= date.end)
# Id.analyzed.v <- as.character(metadata$Accession.ID)

metadata.analyzed <- metadata.analyzed %>%
  left_join(mut.info, by=c("Accession.ID" = "Id"))

metadata.analyzed <- metadata.analyzed %>%
  mutate(hap = ifelse(is.na(mut), "hap_0", mut))
metadata.analyzed$hap %>% table()

region_analyzed.v <- unique(metadata$region_analyzed)
metadata.analyzed  %>% group_by(region_analyzed, mut) %>% summarise(count=n()) %>% filter(count>50, !is.na(mut))
region_analyzed.v <- c('Europe', 'United Kingdom', 'USA', 'Asia')
# metadata.analyzed  %>% group_by(region_analyzed, mut) %>% summarise(count=n()) %>% arrange(desc(count)) %>% head(20)

hap.ref <- 'hap_0'

stat.info.sum <- data.frame()
plot_observed.l <- list()
plot.growth.l <- list()
for(region.interest in region_analyzed.v){
  print(region.interest)
  metadata.analyzed.region <- metadata.analyzed %>% filter(region_analyzed == region.interest)
  num.variant.ref <- metadata.analyzed.region %>% filter(hap == hap.ref) %>% nrow()

  total.analyzed <- nrow(metadata.analyzed.region)

  count.analyzed.region.df <- metadata.analyzed.region %>% group_by(hap) %>% summarize(count.analyzed = n())
  count.analyzed.region.df.filtered <- count.analyzed.region.df %>% filter(count.analyzed > limit.count.analyzed) %>% ungroup()
  print(count.analyzed.region.df)
  variant.interest.v <- count.analyzed.region.df.filtered$hap %>% as.character()

  metadata.analyzed.region.selected <- metadata.analyzed.region %>% filter(hap %in% variant.interest.v)

  ######Transmissibility######
  metadata.filtered.interest <- metadata.analyzed.region.selected %>% mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1, date.bin = cut(date.num,seq(0,max(date.num),bin.size)), date.bin.num = as.numeric(date.bin))
  metadata.filtered.interest <- metadata.filtered.interest %>% filter(!is.na(date.bin))

  metadata.filtered.interest.bin <- metadata.filtered.interest %>% group_by(date.bin.num,hap) %>% summarize(count = n()) %>% ungroup()

  metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin %>% spread(key=hap,value = count)
  metadata.filtered.interest.bin.spread[is.na(metadata.filtered.interest.bin.spread)] <- 0

  X <- as.matrix(data.frame(X0 = 1, X1 = metadata.filtered.interest.bin.spread$date.bin.num))

  Y <- metadata.filtered.interest.bin.spread %>% select(- date.bin.num) %>% select(all_of(variant.interest.v))

  count.group <- apply(Y,2,sum)
  count.total <- sum(count.group)
  prop.group <- count.group / count.total

  Y <- Y %>% as.matrix()

  group.df <- data.frame(group_Id = 1:ncol(Y), group = colnames(Y))

  Y_sum.v <- apply(Y,1,sum)
  dim(Y)

  data.stan <- list(K = ncol(Y),
                    N = nrow(Y),
                    D = 2,
                    X = X,
                    Y = Y,
                    generation_time = generation_time,
                    bin_size = bin.size,
                    Y_sum = Y_sum.v)

  fit.stan <- multi_nomial_model$sample(
    data=data.stan,
    iter_sampling=1000,
    iter_warmup=1900,
    seed=1234,
    parallel_chains = 2,
    #adapt_delta = 0.99,
    max_treedepth = 20,
    chains=4)


  #growth rate
  stat.info <- fit.stan$summary("growth_rate") %>% as.data.frame()
  stat.info$hap <- colnames(Y)[2:ncol(Y)]

  stat.info.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.005,0.995))) %>% as.data.frame() %>% rename(q0.5 = `0.5%`, q99.5 = `99.5%`)
  stat.info <- stat.info %>% inner_join(stat.info.q,by="variable")

  stat.info <- stat.info %>% mutate(signif = ifelse(q0.5 > 1,'higher','not higher'))
  stat.info <- stat.info %>% arrange(desc(hap))
  # stat.info <- stat.info %>% mutate(Pango.lineage = factor(Pango.lineage))
  #stat.info <- stat.info %>% mutate(Pango.lineage2 = factor(Pango.lineage2,levels=Pango.lineage2))

  stat.info.sum <- rbind(stat.info.sum,stat.info)

  #growth rate plot
  # variant.top10.v <- stat.info %>% top_n(10,mean) %>% arrange(desc(mean)) %>% pull(Pango.lineage2) %>% as.character()
  # variant.top10.v <- c(variant.ref,variant.top10.v)

  # ylabel <- paste('relative growth rate per generation (per ',variant.ref,')',sep="")
  ylabel <- str_c('relative growth rate per generation (per ',mut.interest ,')')

  # g <- ggplot(stat.info %>% filter(Pango.lineage2 %in% variant.top10.v),aes(x=Pango.lineage2,y=mean,fill=signif))
  g <- stat.info  %>% ggplot(aes(x=hap,y=mean,fill=signif))
  g <- g + geom_bar(stat = "identity")
  g <- g + geom_errorbar(aes(ymin= q5, ymax = q95), width = 0.2)
  g <- g + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
  g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
  g <- g + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line.x = element_line(color="black", size=0.5, lineend = "square"),
                 axis.line.y = element_line(color="black", size=0.5, lineend = "square"),
                 axis.ticks = element_line(color="black", size=0.5, lineend = "square"),
                 axis.text = element_text(colour = "black"),
                 strip.text = element_text(size=8))
  g <- g + coord_flip()
  g <- g + scale_fill_manual(breaks=c("higher","not higher"),values=c("lightsalmon","gray70"))
  g <- g + xlab('') + ylab(ylabel)
  g <- g + ggtitle(region.interest)

  plot.growth.l[[region.interest]] <-  g

  # #plot observed

  metadata.plot <- metadata.analyzed.region %>%
    mutate(group = factor(hap,levels=variant.interest.v))

  col.v <- brewer.pal(n=5, name = "Paired")
  col.selected.v <- col.v[1:length(variant.interest.v)]
  #col.selected.v <- c(col.v[1:(nrow(stat.info)+1)], "grey70")

  g <- ggplot(metadata.plot,aes(x=Collection.date,fill=group))
  g <- g + geom_bar(stat = 'count', position = "identity")
  g <- g + scale_x_date(date_labels = "%y-%m", date_breaks = "2 months", date_minor_breaks = "2 month")
  g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
  g <- g + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line.x = element_line(color="black", size=0.5, lineend = "square"),
                 axis.line.y = element_line(color="black", size=0.5, lineend = "square"),
                 axis.ticks = element_line(color="black", size=0.5, lineend = "square"),
                 axis.text = element_text(colour = "black"),
                 strip.text = element_text(size=8)
  )
  g <- g + scale_fill_manual(values=col.selected.v)
  g <- g + ggtitle(region.interest)
  g

  plot_observed.l[[region.interest]] <- g
}


#output
pdf.observed.name <- paste(out.prefix, paste("alpha.method1.observed",mut.interest, date.start, date.end,"pdf",sep="."), sep="/")
pdf.growth.rate.name <- paste(out.prefix, paste("alpha.method1.growth_rate",mut.interest, date.start, date.end,"pdf",sep="."), sep="/")
txt.growth.rate.name <- paste(out.prefix, paste("alpha.method1.growth_rate",mut.interest, date.start, date.end,"txt",sep="."), sep="/")
# pdf.observed.name <- paste(out.prefix, paste("N501Y.method1.observed",mut.interest, date.start, date.end,"pdf",sep="."), sep="/")
# pdf.growth.rate.name <- paste(out.prefix, paste("N501Y.method1.growth_rate",mut.interest, date.start, date.end,"pdf",sep="."), sep="/")
# txt.growth.rate.name <- paste(out.prefix, paste("N501Y.method1.growth_rate",mut.interest, date.start, date.end,"txt",sep="."), sep="/")


#write output
write.table(stat.info.sum,txt.growth.rate.name,col.names=T,row.names=F,sep="\t",quote=F)

#plot growth rate
# pdf(pdf.growth.rate.name,width=25,height=15)
# wrap_plots(plot.growth.l,ncol=3)
# dev.off()

pdf(pdf.growth.rate.name,width=10,height=4)
wrap_plots(plot.growth.l,ncol=1)
dev.off()

#plot observed
pdf(pdf.observed.name,width=15, height=8)
wrap_plots(plot_observed.l,ncol=2)
dev.off()



