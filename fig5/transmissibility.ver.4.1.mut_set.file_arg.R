#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(ggplot2)
library(rbin)
library(cmdstanr)
library(patchwork)
library(RColorBrewer)


##########args##########
#input
args = commandArgs(trailingOnly=T)

download.date <- args[1] #"2022-04-07"
stan_f.name <- args[2]
metadata.name <- args[3] #'/media/jampei/C6700DA7700D9F75/data/Sato_analysis/GISAID_data/2022_04_07/metadata.tsv' #args[4]
mut.info.name <- args[4] #'/media/jampei/C6700DA7700D9F75/data/Sato_analysis/GISAID_data/2022_04_07/metadata.mut_long.tsv'
out.prefix <- args[5] #'/media/jampei/C6700DA7700D9F75/data/Sato_analysis/transmissibility_auto/GitHub/next_variant_detection/method4/test/test.2022_04_07' #args[5]

# download.date <- "2022_09_26"
# stan_f.name <- "/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/script/multinomial_independent.stan"
# metadata.name <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/", download.date,  "/metadata_", download.date, ".tsv")
# mut.info.name <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/", download.date, "/mutation.aa.all_protein.long.", download.date, ".txt")
# out.prefix <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/output/", download.date, "/transmissibility")

mut.interest.v <- c("Spike_L452", "Spike_T478", "Spike_E484", "Spike_N501")
##convert data type
download.date <- gsub( "_", "-", download.date)
download.date <- as.Date(download.date)

##########parameters##########
#general

#period to be analyzed
day.delay <- 30
day.analyzed <- 600

date.end <- download.date - day.delay
date.start <- date.end - day.analyzed + 1

#min numbers
limit.count.analyzed <- 50

#Transmissibility
bin.size <- 7
generation_time <- 5.0
# generation_time <- 2.1


country.interest.v <- c('USA','United Kingdom')



#################output path##################
pdf.observed.name <- paste(out.prefix,"/method4.",date.start, "_", date.end, ".RBM_mut_set.observed.pdf",sep="")
pdf.growth.rate.name <- paste(out.prefix,"/method4.",date.start, "_", date.end, ".RBM_mut_set.growth_rate.pdf",sep="")
txt.growth.rate.name <- paste(out.prefix,"/method4.",date.start, "_", date.end, ".RBM_mut_set.growth_rate.txt",sep="")


#################set model#################
#set.seed(1)
#cmdstanpath <- "/usr/local/cmdstan"
#set_cmdstan_path("/Users/yuuka.m/.cmdstan/cmdstan-2.29.2")
multi_nomial_model <- cmdstan_model(stan_f.name)

#################data preprocessing & QC##############
metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)

metadata <- metadata %>%
    filter(Host == "Human",
           str_length(Collection.date) == 10,
           Pango.lineage != "",
           Pango.lineage != "None",
           !N.Content > 0.02,
           !str_detect(Additional.location.information,"[Qq]uarantine"))

metadata <- metadata %>% select(Accession.ID, Collection.date, Location)

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

metadata.analyzed <- metadata %>%
    filter(Collection.date >= date.start, Collection.date <= date.end)
nrow(metadata.analyzed)

#mutation data
# mut.info <- fread(mut.info.name,header=T,sep="\t",check.names=T)
mut.info <- read_tsv(mut.info.name)
D614G.Id.v <- mut.info %>% filter(mut == "Spike_D614G") %>% pull(Id)

metadata.analyzed <-  metadata.analyzed %>% filter(Accession.ID %in% D614G.Id.v)
rm(D614G.Id.v)
nrow(metadata.analyzed)
Id.analyzed.v <- as.character(metadata.analyzed$Accession.ID)
mut.info <- mut.info %>% filter(Id %in% Id.analyzed.v)

#mut.info.interest <- mut.info %>% filter(str_detect(mut,mut.interest))
for(name in mut.interest.v){
    mut.info.RBM <- mut.info %>% filter(str_detect(mut, name)) %>% rename(!!name := mut)
    metadata.analyzed <- metadata.analyzed %>% left_join(mut.info.RBM, by = c("Accession.ID"="Id"))
}
rm(mut.info.RBM)
# metadata.analyzed <- metadata %>%
#     left_join(mut.info.interest %>% rename(Accession.ID = Id), by = "Accession.ID")

metadata.analyzed$Spike_L452 <- ifelse(is.na(metadata.analyzed$Spike_L452), "_", metadata.analyzed$Spike_L452)
metadata.analyzed$Spike_T478 <- ifelse(is.na(metadata.analyzed$Spike_T478), "_", metadata.analyzed$Spike_T478)
metadata.analyzed$Spike_E484 <- ifelse(is.na(metadata.analyzed$Spike_E484), "_", metadata.analyzed$Spike_E484)
metadata.analyzed$Spike_N501 <- ifelse(is.na(metadata.analyzed$Spike_N501), "_", metadata.analyzed$Spike_N501)


metadata.analyzed <- metadata.analyzed %>%
    mutate(mut_set=paste(Spike_L452,Spike_T478,Spike_E484,Spike_N501, sep = "/"))
metadata.analyzed

##########transmissibility estimation##########

region_analyzed.v <- unique(metadata$region_analyzed)

stat.info.sum <- data.frame()
plot_observed.l <- list()
plot.growth.l <- list()

variant.ref <- "_/_/_/_"


for(region.interest in region_analyzed.v){
    # region.interest <- region_analyzed.v[2]
    print(region.interest)
    metadata.analyzed.region <- metadata.analyzed %>% filter(region_analyzed == region.interest)
    num.variant.ref <- metadata.analyzed.region %>% filter(mut_set == variant.ref) %>% nrow()
    print(num.variant.ref)
    total.analyzed <- nrow(metadata.analyzed.region)

    count.analyzed.region.df <- metadata.analyzed.region %>% group_by(mut_set) %>% summarize(count.analyzed = n(), prop.analyzed = count.analyzed / total.analyzed)
    count.analyzed.region.df.filtered <- count.analyzed.region.df %>% filter(count.analyzed > limit.count.analyzed & prop.analyzed > 0.0005) %>% arrange(desc(count.analyzed))

    variant.interest.v <- count.analyzed.region.df.filtered$mut_set %>% as.character()
    count.top10.v  <- variant.interest.v[1:10]

    if(variant.ref %in% variant.interest.v){
        variant.interest.v <- c(variant.ref,variant.interest.v[-which(variant.interest.v == variant.ref)])
    } else {
        print("not enough reference samples for the analysis")
        variant.interest.v <- c(variant.ref,variant.interest.v)
    }
    #kill
    if(length(variant.interest.v) < 2) next
    if(num.variant.ref < limit.count.analyzed) next

    metadata.analyzed.region.selected <- metadata.analyzed.region %>% filter(mut_set %in% variant.interest.v)

    ######Transmissibility
    bin.size <- 7
    metadata.filtered.interest <- metadata.analyzed.region.selected %>%
        mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1,
            date.bin = cut(date.num,seq(0,max(date.num),bin.size)),
            date.bin.num = as.numeric(date.bin))
    metadata.filtered.interest <- metadata.filtered.interest %>% filter(!is.na(date.bin))

    metadata.filtered.interest.bin <- metadata.filtered.interest %>% group_by(date.bin.num,mut_set) %>% summarize(count = n()) %>% ungroup()

    metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin %>% spread(key=mut_set,value = count)
    metadata.filtered.interest.bin.spread[is.na(metadata.filtered.interest.bin.spread)] <- 0

    X <- as.matrix(data.frame(X0 = 1, X1 = metadata.filtered.interest.bin.spread$date.bin.num))

    Y <- metadata.filtered.interest.bin.spread %>% select(- date.bin.num) %>% select(all_of(variant.interest.v))

    count.group <- apply(Y,2,sum)
    count.total <- sum(count.group)
    prop.group <- count.group / count.total

    Y <- Y %>% as.matrix()

    group.df <- data.frame(group_Id = 1:ncol(Y), group = colnames(Y))

    Y_sum.v <- apply(Y,1,sum)

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
        iter_warmup=1500,
        seed=1234,
        parallel_chains = 4,
        #adapt_delta = 0.99,
        max_treedepth = 20,
        chains=4)


    #growth rate
    stat.info <- fit.stan$summary("growth_rate") %>% as.data.frame()
    print(stat.info$rhat %>% sort() %>% tail())
    stat.info$mut_set <- colnames(Y)[2:ncol(Y)]

    stat.info.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.005,0.995))) %>% as.data.frame() %>% rename(q0.5 = `0.5%`, q99.5 = `99.5%`)
    stat.info <- stat.info %>% inner_join(stat.info.q,by="variable")

    stat.info <- stat.info %>% mutate(signif = ifelse(q0.5 > 1,'positive',ifelse(q99.5 < 1,'negative', 'not signif')))
    stat.info <- stat.info %>% arrange(mean)
    stat.info <- stat.info %>% mutate(mut_set = factor(mut_set,levels=mut_set))

    stat.info.sum <- rbind(stat.info.sum,stat.info)

    ylabel <- paste('relative growth rate per generation (per ',variant.ref,')',sep="")

    g <- ggplot(stat.info,aes(x=mut_set,y=mean,fill=signif))
    g <- g + geom_bar(stat = "identity")
    g <- g + geom_errorbar(aes(ymin= q5, ymax = q95), width = 0.2)
    g <- g + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
    g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
    g <- g + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              strip.text = element_text(size=8))
    g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
    g <- g + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              strip.text = element_text(size=8))
    g <- g + coord_flip()
    g <- g + scale_fill_manual(breaks=c("positive","negative",'not signif'),values=c("lightsalmon","lightskyblue","gray70"))
    g <- g + xlab('') + ylab(ylabel)
    g <- g + ggtitle(region.interest)

    plot.growth.l[[region.interest]] <- g


    #plot observed
    #growth rate plot
    variant.top10.v <- stat.info %>% filter(signif=="positive") %>%
        top_n(12,mean) %>% arrange(desc(mean)) %>% pull(mut_set) %>% as.character()
    variant.top10.v <- c(variant.top10.v, variant.ref)
    metadata.plot <- metadata.analyzed.region %>%
        mutate(group = ifelse(mut_set %in% variant.top10.v,as.character(mut_set),"others")) %>%
        mutate(group = factor(group,levels=c(variant.top10.v,'others')))

    col.v <- c(brewer.pal(length(variant.top10.v)-1, "Paired"), "grey30", "gray70")

    g <- ggplot(metadata.plot,aes(x=Collection.date,fill=group))
    g <- g + geom_bar(stat = 'count', position = "fill")
    g <- g + scale_x_date(date_labels = "%y-%m", date_breaks = "1 months", date_minor_breaks = "1 month")
    g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
    g <- g + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              strip.text = element_text(size=8)
        )
    g <- g + scale_fill_manual(values=col.v)
    g <- g + ggtitle(region.interest)

    plot_observed.l[[region.interest]] <- g

}

#write output
write.table(stat.info.sum,txt.growth.rate.name,col.names=T,row.names=F,sep="\t",quote=F)

#plot growth rate
pdf(pdf.growth.rate.name,width=20,height=20)
wrap_plots(plot.growth.l,ncol=2)
dev.off()

#plot observed
pdf(pdf.observed.name,width=30,height=20)
wrap_plots(plot_observed.l,ncol=2)
dev.off()

