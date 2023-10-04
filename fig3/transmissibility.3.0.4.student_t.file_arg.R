#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(cmdstanr)
library(tictoc)
library(circlize)
library(igraph)
library(makedummies)

#input
args = commandArgs(trailingOnly=T)
download.date <- args[1]
metadata.name <- args[2]
mut.info.name <- args[3]
out.prefix <- args[4]
day.delay <- as.numeric(args[5])
day.analyzed <- as.numeric(args[6])
region.interest <- args[7]
stan_f.name <- args[8]

# download.date <- "2022_03_23"
# metadata.name <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/", download.date,  "/metadata_", download.date, ".tsv")
# mut.info.name <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/gisaid/", download.date, "/mutation.aa.all_protein.long.", download.date, ".txt")
# out.prefix <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/output/", download.date, "/transmissibility/student_t")
# day.delay <- 200
# day.analyzed <- 300
# region.interest <- "USA"
# stan_f.name <- "/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/script/multinomial_mut_regression.mut_only.2.stan"




##convert data type
download.date <- gsub( "_", "-", download.date)
download.date <- as.Date(download.date)

##########parameters##########
#general
# core.num <- 1
S_w.value <- 0.5
S_w.value <- 1

#download.date <- as.Date(download.date)
date.end <- download.date - day.delay
date.start <- date.end - day.analyzed + 1

#transmissibility
#variant.ref <- "BA.1"
bin.size <- 3
# generation_time <- 2.1
generation_time <- 5.0
country.interest.v <- c('USA','United Kingdom')

tree.depth <- 15


#################output path##################
# pdf.growth.rate.name <- paste(out.prefix,"/method3.spike_hap.",date.start, "_", date.end, ".USA.growth_rate.pdf", sep="")
# pdf.growth.gain.name <- paste(out.prefix,"/method3.spike_hap.",date.start, "_", date.end, ".USA.growth_gain.pdf", sep="")
# pdf.plot_observed.name <- paste(out.prefix,"/method3.spike_hap.",date.start, "_", date.end, ".USA.variant_freq.observed.pdf", sep="")
# pdf.growth.gain_summary.name <- paste(out.prefix,"/method3.spike_hap.",date.start, "_", date.end, ".USA.growth_gain.summary_each_lineage.pdf", sep="")
txt.metadata_out.name <- paste(out.prefix, "/method3.spike_hap.student_t.",date.start, "_", date.end, ".",region.interest,".metadata.txt", sep="")
txt.growth.rate.name <- paste(out.prefix,"/method3.spike_hap.student_t.",date.start, "_", date.end, ".",region.interest,".tree_", tree.depth ,".bin_", bin.size, ".Sw_", S_w.value, ".growth_rate.txt", sep="")
txt.growth.gain.name <- paste(out.prefix,"/method3.spike_hap.student_t.",date.start, "_", date.end, ".",region.interest,".tree_", tree.depth ,".bin_", bin.size, ".Sw_", S_w.value, ".growth_gain.txt", sep="")

region.interest <- gsub("_", " ", region.interest)

# ##########data preprocessing & QC##########
#metadata
metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)

metadata <- metadata %>%
  filter(Host == "Human",
         str_length(Collection.date) == 10,
         !N.Content > 0.02,
         !str_detect(Additional.location.information,"[Qq]uarantine"))

metadata <- metadata %>% select(Accession.ID, Collection.date, Location)

metadata <- metadata %>%
  mutate(Collection.date = as.Date(Collection.date),
         region = str_split(Location," / ",simplify = T)[,1],
         country = str_split(Location," / ",simplify = T)[,2]) %>%
  select(-Location)


country.interest.v <- c('USA','United Kingdom')

metadata <- metadata %>%
  mutate(region_analyzed = ifelse(country %in% country.interest.v,
                                  as.character(country),
                                  as.character(region))) %>%
  select(-country, -region)

####extract USA data####
metadata.analyzed <- metadata %>%
  filter(region_analyzed == region.interest) %>%
  filter(Collection.date >= date.start, Collection.date <= date.end)  %>%
  select(Accession.ID, Collection.date)

#metadata.analyzed.Id.v <- metadata.analyzed$Accession.ID

#mutation data
mut.info <- read_tsv(mut.info.name)

#extract data which have D614G
Id.D614G.v <- mut.info %>% filter(mut=="Spike_D614G") %>% pull(Id)
Id.P323L.v <- mut.info %>% filter(mut=="NSP12_P323L") %>% pull(Id)
Id.D614G_P323L.v <- intersect(Id.D614G.v, Id.P323L.v)
metadata.analyzed <- metadata.analyzed %>% filter(Accession.ID %in% Id.D614G_P323L.v)
rm(Id.D614G.v)
rm(Id.P323L.v)
mut.info <- mut.info %>% filter(Id %in% metadata.analyzed$Accession.ID)

mut.count.df <- mut.info %>% group_by(mut) %>% summarize(count.mut = n())
mut.count.df  %>% arrange(desc(count.mut))
# mut.count.df.filtered <- mut.count.df %>% filter(count.mut>500)
# mut.count.df.filtered <- mut.count.df %>% filter(count.mut>1000)
mut.count.df.filtered <- mut.count.df %>% arrange(desc(count.mut)) %>% head(500)

mut.analyzed.v <- mut.count.df.filtered %>% pull(mut) %>% unique() %>% sort()
length(mut.analyzed.v)
mut.info <- mut.info %>% filter(mut %in% mut.analyzed.v)

mut.info.filtered.spread <- tibble(Id = metadata.analyzed$Accession.ID) %>%
  left_join(mut.info) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from=mut,values_from=value)
mut.info.filtered.spread[is.na(mut.info.filtered.spread)] <- 0
mut.info.filtered.spread <- mut.info.filtered.spread %>% select(Id,all_of(mut.analyzed.v))

#####make haplotype######
mut.sd.gt0.exclude.v <- apply(mut.info.filtered.spread[,2:ncol(mut.info.filtered.spread)],2,sd)
mut.sd.gt0.exclude.v <- mut.sd.gt0.exclude.v[mut.sd.gt0.exclude.v==0] %>% names()
mut.info.filtered.spread <- mut.info.filtered.spread %>% select(-all_of(mut.sd.gt0.exclude.v))
mut.analyzed.v <- mut.analyzed.v[!(mut.analyzed.v %in%  mut.sd.gt0.exclude.v)]
mut.analyzed.v   %>% length()

hap.v <- apply(mut.info.filtered.spread[,2:ncol(mut.info.filtered.spread)],1,paste,collapse="")
hap.v %>% unique() %>% length()
hap_num.v <- as.character(as.numeric(as.factor(hap.v)))

hap.df <- data.frame(Accession.ID = mut.info.filtered.spread$Id, hap = hap_num.v)
hap.df <- hap.df %>% mutate(hap_Id = str_c("hap_", hap))
hap_Id.count.df <- hap.df %>% group_by(hap_Id) %>% summarise(hap_Id.count=n()) %>% arrange(desc(hap_Id.count))
hap_Id.count.df %>% filter(hap_Id == "hap_1")

# hap_Id.interest.v <- hap_Id.count.df %>% filter(hap_Id.count > 5) %>% pull(hap_Id)
hap_Id.interest.v <- hap_Id.count.df %>% pull(hap_Id)
mut.info.filtered.spread <- mut.info.filtered.spread %>% filter(Id %in% (hap.df %>% filter(hap_Id %in% hap_Id.interest.v) %>% pull(Accession.ID)))
mut.info.filtered.spread <- mut.info.filtered.spread %>% left_join(hap.df %>% select(Id=Accession.ID,hap_Id))
mut.info.filtered.spread.hap <- mut.info.filtered.spread[!duplicated(mut.info.filtered.spread[,2:ncol(mut.info.filtered.spread)]),] %>% select(-Id)

nrow(mut.info.filtered.spread)
ncol(mut.info.filtered.spread)

hap.mut.diff.mat.cor <- cor(mut.info.filtered.spread.hap[,1:(ncol(mut.info.filtered.spread.hap)-1)])
hap.mut.diff.mat.cor.long <- hap.mut.diff.mat.cor %>%
  as.data.frame() %>%
  mutate(mut1 = rownames(hap.mut.diff.mat.cor)) %>%
  gather(key=mut2,value=cor,-mut1) %>%
  filter(mut1 != mut2,cor>0.90)
rm(hap.mut.diff.mat.cor)
g <- graph.data.frame(hap.mut.diff.mat.cor.long[,1:2],directed=F)
g <- simplify(g,remove.multiple=T,remove.loops=T)

# plot(g,vertex.size=4)

connected.l <- decompose(g)

mut_group.df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
for(i in 1:length(connected.l)){
  v <- V(connected.l[[i]])$name
  if ("Spike_L452R" %in% v){
    print(v)
    v <- v[-which(v == "Spike_L452R")]
  }
  if ("Spike_T478K" %in% v){
    print(v)
    v <- v[-which(v == "Spike_T478K")]
  }
  if ("Spike_E484K" %in% v){
    print(v)
    v <- v[-which(v == "Spike_E484K")]
  }
  if ("Spike_N501Y" %in% v){
    print(v)
    v <- v[-which(v == "Spike_N501Y")]
  }
  name <- paste(v,collapse=":")
  temp.df <- data.frame(mut = v, mut_group = name)
  mut_group.df <- rbind(mut_group.df,temp.df)
}

mut.info.filtered.long <- mut.info.filtered.spread.hap %>% gather(key=mut, value=value, -hap_Id)
mut.info.filtered.long <- mut.info.filtered.long %>% left_join(mut_group.df,by="mut")
mut.info.filtered.long <- mut.info.filtered.long %>%
  mutate(mut_group = ifelse(is.na(mut_group),as.character(mut),as.character(mut_group)))

mut.info.filtered.long.mut_group <- mut.info.filtered.long %>%
  group_by(hap_Id,mut_group) %>% summarize(value = mean(value)) %>% ungroup()
rm(mut.info.filtered.long)

mut.info.filtered.long.mut_group.mat <- mut.info.filtered.long.mut_group %>%
  spread(key=mut_group,value=value)

metadata.analyzed.merged <- merge(metadata.analyzed,hap.df %>% select(Accession.ID, hap_Id),by="Accession.ID")

##output metadata##
# write.table(metadata.analyzed.merged,txt.metadata_out.name,col.names=T,row.names=F,sep="\t",quote=F)

count.hap.df <- metadata.analyzed.merged %>% group_by(hap_Id) %>% summarise(hap_Id.count=n()) %>% arrange(desc(hap_Id.count))
print(region.interest)
count.hap.df %>% filter(hap_Id=="hap_1")
# count.hap.df %>% arrange(desc(hap_Id.count)) %>% head(500) %>% pull(hap_Id.count) %>% sum()
# count.hap.df %>% filter(hap_Id.count>50)
# count.hap.interest.df <- count.hap.df %>% arrange(desc(hap_Id.count)) %>% head(500)
count.hap_1 <- count.hap.df %>% filter(hap_Id=="hap_1") %>% pull(hap_Id.count)
if (count.hap_1>100){
  count.hap.interest.df <- count.hap.df %>% arrange(desc(hap_Id.count)) %>% filter(hap_Id.count>100)
}else{
  if(count.hap_1>50){
    count.hap.interest.df <- count.hap.df %>% arrange(desc(hap_Id.count)) %>% filter(hap_Id.count>50)
  }else{
    stop("not enough reference seq number")
  }
}
count.hap.interest.df  %>% pull(hap_Id.count) %>% sum()

hap_Id.interest.v <- as.character(count.hap.interest.df$hap_Id)

hap_Id.reference <- "hap_1"
# hap_Id.interest.v <- as.character(mut.info.filtered.long.mut_group.mat$hap_Id.2)

####sort columns####
hap_Id.interest.v <- c(hap_Id.reference,hap_Id.interest.v[hap_Id.interest.v!= hap_Id.reference])
mut.mat <- mut.info.filtered.long.mut_group.mat[match(hap_Id.interest.v, mut.info.filtered.long.mut_group.mat$hap_Id),]
rm(mut.info.filtered.long.mut_group.mat)

mut.sd.gt0.interest.v <- apply(mut.mat %>% select(-hap_Id),2,sd)
mut.sd.gt0.interest.v <- mut.sd.gt0.interest.v[mut.sd.gt0.interest.v>0] %>% names()
ncol(mut.mat)
mut.sd.gt0.interest.v %>% length()

mut.mat <- mut.mat %>% select(all_of(mut.sd.gt0.interest.v),hap_Id)
COEF <- mut.mat %>% select(-hap_Id)


# hap_Id.interest.v <- c(hap_Id.reference, hap_Id.interest.v[-which(hap_Id.interest.v == hap_Id.reference)])
metadata.analyzed.interest <- metadata.analyzed.merged %>% filter(hap_Id　%in% hap_Id.interest.v)
metadata.analyzed.interest <- metadata.analyzed.interest %>%
  mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))+ 1,
         date.bin = cut(date.num,seq(0,max(date.num),bin.size)),
         date.bin.num = as.numeric(date.bin)
  )

metadata.analyzed.interest <- metadata.analyzed.interest %>% filter(!is.na(date.bin))

metadata.analyzed.interest.bin <- metadata.analyzed.interest %>% group_by(date.bin.num,hap_Id) %>% summarize(count = n()) %>% ungroup()
metadata.analyzed.interest.bin.spread <- metadata.analyzed.interest.bin %>% spread(key=hap_Id,value = count)
metadata.analyzed.interest.bin.spread[is.na(metadata.analyzed.interest.bin.spread)] <- 0
metadata.analyzed.interest.bin.spread <- metadata.analyzed.interest.bin.spread %>% select(date.bin.num,all_of(hap_Id.interest.v))
rm(metadata.analyzed.merged)
rm(metadata.analyzed.interest.bin)

TS <- as.matrix(data.frame(X0 = 1, X1 = metadata.analyzed.interest.bin.spread$date.bin.num))

Y <- metadata.analyzed.interest.bin.spread %>% select(- date.bin.num)
# group.df <- data.frame(group_Id = 1:ncol(Y), hap_Id = colnames(Y))

# Y %>% apply(2, sum)

Y <- Y %>% as.matrix()
Y_sum.v <- apply(Y,1,sum)

dim(COEF)
dim(Y)

num.coef <- ncol(COEF) #+ ncol(clade.dummy.mat)

#data for stan
data.stan <- list(K = ncol(Y),
                  N = nrow(Y),
                  D = num.coef,
                  S_w = S_w.value,
                  TS = TS,
                  COEF = COEF,
                  Y = Y,
                  generation_time = generation_time,
                  bin_size = bin.size,
                  Y_sum = Y_sum.v)

multi_nomial_model <- cmdstan_model(stan_f.name)

#stan fitting
fit.stan <- multi_nomial_model$sample(
  data=data.stan,
  iter_sampling=1000,
  iter_warmup=500,
  seed=1234,
  parallel_chains = 4,
  #adapt_delta = 0.999,
  max_treedepth = tree.depth,
  chains=4)

fit.stan$summary("S_b")
stat.info.growth_gain <- fit.stan$summary("growth_gain") %>% as.data.frame()
stat.info.growth_rate <- fit.stan$summary("growth_rate") %>% as.data.frame()

stat.info.growth_gain.q <- fit.stan$summary("growth_gain", ~quantile(.x, probs = c(0.005,0.995))) %>% as.data.frame() %>% rename(q0.5 = `0.5%`, q99.5 = `99.5%`)
stat.info.growth_rate.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.005,0.995))) %>% as.data.frame() %>% rename(q0.5 = `0.5%`, q99.5 = `99.5%`)

stat.info.growth_gain <- stat.info.growth_gain %>% inner_join(stat.info.growth_gain.q,by="variable")
stat.info.growth_rate <- stat.info.growth_rate %>% inner_join(stat.info.growth_rate.q,by="variable")

stat.info.growth_rate$hap_Id <- colnames(Y)
stat.info.growth_gain$feature <- colnames(COEF)


#growth gain
stat.info.growth_gain <- stat.info.growth_gain %>% mutate(signif = ifelse(q0.5 > 1,'positive',
                                                                          ifelse(q99.5 < 1,'negative', 'not signif')))
write.table(stat.info.growth_gain,txt.growth.gain.name,col.names=T,row.names=F,sep="\t",quote=F)

#growth rate
stat.info.growth_rate <- stat.info.growth_rate
stat.info.growth_rate <- stat.info.growth_rate %>%
  mutate(signif = ifelse(q0.5 > 1,'positive',ifelse(q99.5 < 1,'negative', 'not signif')))

stat.info.growth_rate <- stat.info.growth_rate %>%
  arrange(mean) %>% mutate(hap_Id = factor(hap_Id, levels=hap_Id))
mut.diff.sum.v <- c()

for(hap_Id.interest in stat.info.growth_rate$hap_Id){
  mut.diff.v <- mut.info.filtered.long.mut_group %>%
    filter(hap_Id == hap_Id.interest, value > 0) %>% pull(mut_group) %>% unique() %>% sort()
  mut.diff <- paste(mut.diff.v,collapse=",")
  #mut.diff <- gsub("Spike_", "", mut.diff )
  mut.diff.sum.v <- c(mut.diff.sum.v,mut.diff)
}

stat.info.growth_rate$label <- mut.diff.sum.v
write.table(stat.info.growth_rate,txt.growth.rate.name,col.names=T,row.names=F,sep="\t",quote=F)
