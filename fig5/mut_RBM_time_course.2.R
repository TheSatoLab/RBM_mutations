library(tidyverse)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(scales)
library(circlize)
options(scipen=1000)

### dir="/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/"
###
# ${R_path} --vanilla --slave --args \
#     ${input_fasta_prefix} \
#     Spike_L452 Spike_T478 Spike_E484 Spike_N501 < \
#     script/mut_RBM_time_course.2.R
###

input_fasta_prefix <-  commandArgs()[5]
locus_list <- commandArgs()[6:length(commandArgs())]


metadata_mut_path <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/output/", input_fasta_prefix, "/metadata_", input_fasta_prefix, ".mut_RBM.tsv")
out_plot_mut_path <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/output/", input_fasta_prefix, "/histgram_mut_rate.metadata_", input_fasta_prefix, ".mut_RBM.modified.pdf")
out_plot_strain_path <- str_c("/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/output/", input_fasta_prefix,"/histgram_mut_rate.metadata_", input_fasta_prefix, ".strain.pdf")
# locus_list <- c("Spike_L452", "Spike_T478", "Spike_E484", "Spike_N501")


metadata_mut <- read_tsv(metadata_mut_path)
metadata_mut$date[nchar(metadata_mut$date)!=10] <- NA
metadata_mut$date <- as.Date(metadata_mut$date)
first_date  <-  min(na.omit(as.numeric(metadata_mut$date)))
metadata_mut  <- metadata_mut %>% mutate(date.num = as.numeric(date) - first_date)
metadata_mut <- filter(metadata_mut , !is.na(date.num), !is.na(pango_lineage))
#metadata_mut <- select(metadata_mut, -date, -region_exposure)
#metadata_mut

##each mut##
major_mut.L452.v <- metadata_mut %>% filter(!is.na(Spike_L452)) %>% group_by(Spike_L452) %>%
  summarise(count_L452=n()) %>%
  top_n(3,count_L452) %>%
  arrange(desc(count_L452)) %>% pull(Spike_L452)

major_mut.T478.v <- metadata_mut %>% filter(!is.na(Spike_T478)) %>% group_by(Spike_T478) %>%
  summarise(count_T478=n()) %>%
  top_n(3,count_T478) %>%
  arrange(desc(count_T478)) %>% pull(Spike_T478)

major_mut.E484.v <- metadata_mut %>% filter(!is.na(Spike_E484)) %>% group_by(Spike_E484) %>%
  summarise(count_E484=n()) %>%
  top_n(3,count_E484) %>%
  arrange(desc(count_E484)) %>% pull(Spike_E484)

major_mut.N501.v <- metadata_mut %>% filter(!is.na(Spike_N501)) %>% group_by(Spike_N501) %>%
  summarise(count_N501=n()) %>%
  top_n(3,count_N501) %>%
  arrange(desc(count_N501)) %>% pull(Spike_N501)

pal.base <- brewer.pal(4, "Dark2")
pal1 <- colorRamp2(c(0, 5), c('white', pal.base[1]))
pal2 <- colorRamp2(c(0, 5), c('white', pal.base[2]))
pal3 <- colorRamp2(c(0, 5), c('white', pal.base[3]))
pal4 <- colorRamp2(c(0, 5), c('white', pal.base[4]))

mut_color.df <- tibble(mut=major_mut.L452.v, color=pal1(c(5, 3, 1)))
mut_color.df <- mut_color.df %>% bind_rows(tibble(mut=major_mut.T478.v, color=pal2(c(5, 3, 1))))
mut_color.df <- mut_color.df %>% bind_rows(tibble(mut=major_mut.E484.v, color=pal3(c(5, 3, 1))))
mut_color.df <- mut_color.df %>% bind_rows(tibble(mut=major_mut.N501.v, color=pal4(c(5, 3, 1))))

###modify###
Spike_L452_p.v <- ifelse(metadata_mut$Spike_L452 %in% major_mut.L452.v,
                         metadata_mut$Spike_L452,
                         ifelse(is.na(metadata_mut$Spike_L452), NA, "others"))
Spike_L452_p.v <- factor(Spike_L452_p.v, levels = c(major_mut.L452.v, "others"))

Spike_T478_p.v <- ifelse(metadata_mut$Spike_T478 %in% major_mut.T478.v,
                         metadata_mut$Spike_T478,
                         ifelse(is.na(metadata_mut$Spike_T478), NA, "others"))
Spike_T478_p.v <- factor(Spike_T478_p.v, levels = c(major_mut.T478.v, "others"))

Spike_E484_p.v <- ifelse(metadata_mut$Spike_E484 %in% major_mut.E484.v,
                         metadata_mut$Spike_E484,
                         ifelse(is.na(metadata_mut$Spike_E484), NA, "others"))
Spike_E484_p.v <- factor(Spike_E484_p.v, levels = c(major_mut.E484.v, "others"))

Spike_N501_p.v <- ifelse(metadata_mut$Spike_N501 %in% major_mut.N501.v,
                         metadata_mut$Spike_N501,
                         ifelse(is.na(metadata_mut$Spike_N501), NA, "others"))
Spike_N501_p.v <- factor(Spike_N501_p.v, levels = c(major_mut.N501.v, "others"))


metadata_mut.modified <- tibble(gisaid_epi_isl=metadata_mut$gisaid_epi_isl,
                                date=metadata_mut$date,
                                date.num=metadata_mut$date.num,
                                Spike_L452=Spike_L452_p.v,
                                Spike_T478=Spike_T478_p.v,
                                Spike_E484=Spike_E484_p.v,
                                Spike_N501=Spike_N501_p.v)

my_pal <- mut_color.df %>% filter(str_detect(mut, locus_list[1])) %>% pull(color) %>% c("grey50")
g <- metadata_mut.modified[!is.na(metadata_mut.modified[,4]),] %>% ggplot() +
  geom_histogram(aes_string(x="date", fill=locus_list[1]), binwidth = 10)+
  theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))+
  theme(legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.title =  element_text(colour="black", size=12),
        axis.line.x = element_line(color="black", size=0.5, lineend = "square"),
        axis.line.y = element_line(color="black", size=0.5, lineend = "square"),
        axis.text.y = element_text(colour="black", size=12, angle = 45),
        axis.text.x = element_text(colour="black", size=12),
        axis.ticks = element_line(color="black", size=0.5, lineend = "square"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right")+
  ggtitle(locus_list[1])+
  scale_fill_manual(values = my_pal)+
  scale_y_continuous(labels=scientific)

for (i in 2:length(locus_list)) {
  my_pal <- mut_color.df %>% filter(str_detect(mut, locus_list[i])) %>% pull(color) %>% c("grey50")
  g <- g + metadata_mut.modified[!is.na(metadata_mut.modified[,(3+i)]),] %>% ggplot() +
    geom_histogram(aes_string(x="date", fill=locus_list[i]), binwidth = 10)+
    theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))+
    theme(legend.key.size = unit(0.3, 'cm'),
          legend.text = element_text(size=12),
          legend.title = element_blank(),
          axis.line.x = element_line(color="black", size=0.5, lineend = "square"),
          axis.line.y = element_line(color="black", size=0.5, lineend = "square"),
          axis.title = element_blank(),
          axis.text.y = element_text(colour="black", size=12, angle = 45),
          axis.text.x = element_text(colour="black", size=12),
          axis.ticks = element_line(color="black", size=0.5, lineend = "square"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "right")+
    ggtitle(locus_list[i])+
    scale_fill_manual(values = my_pal)+
    scale_y_continuous(labels=scientific)
}
g <- g + plot_layout(ncol = 2)
ggsave(filename=out_plot_mut_path, plot = g, width=12, height=8)


###stain count###
metadata_mut$pango_lineage <- gsub("(BA\\.[0-9])\\..+", "\\1", metadata_mut$pango_lineage)
metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="B.1.1.7"|str_detect(metadata_mut$pango_lineage, "^Q\\."), "Alpha", "others")
metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="B.1.617.2"|str_detect(metadata_mut$pango_lineage, "^AY\\."), "Delta",metadata_mut$pango_lineage2)
metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="B.1.1.529", "Omicron", metadata_mut$pango_lineage2)
metadata_mut$pango_lineage2 <- ifelse(str_detect(metadata_mut$pango_lineage, "^BA\\."), metadata_mut$pango_lineage, metadata_mut$pango_lineage2)


# metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="P.1"|str_detect(metadata_mut$pango_lineage, "P\\.1\\."), "Gamma", metadata_mut$pango_lineage2)
# metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="B.1.351"|str_detect(metadata_mut$pango_lineage, "B\\.1\\.351"), "Beta", metadata_mut$pango_lineage2)
# metadata_mut$pango_lineage2 <- ifelse(metadata_mut$pango_lineage=="B.1.621" | metadata_mut$pango_lineage=="BB.2" |str_detect(metadata_mut$pango_lineage, "B\\.1\\.621\\."), "Mu", metadata_mut$pango_lineage2)


#metadata_mut$pango_lineage2

metadata_mut %>% group_by(pango_lineage2) %>% summarise(count=n()) %>% arrange(desc(count))

strain_order.v <- metadata_mut$pango_lineage2 %>% unique()  %>% sort()
strain_order.v <- c("Alpha", "Delta", "Omicron", strain_order.v[-which(strain_order.v %in% c("Alpha", "Delta", "Omicron"))])
metadata_mut$pango_lineage2 <- factor(metadata_mut$pango_lineage2, levels = strain_order.v)

my_pal <- c(brewer.pal(length(strain_order.v)-1, "Paired"), "grey50")
g <- metadata_mut %>% ggplot() +
  geom_histogram(aes_string(x="date", fill="pango_lineage2"), binwidth = 10)+
  theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))+
  theme(legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.title =  element_text(colour="black", size=12),
        axis.line.x = element_line(color="black", size=0.5, lineend = "square"),
        axis.line.y = element_line(color="black", size=0.5, lineend = "square"),
        axis.text.y = element_text(colour="black", size=12, angle = 45),
        axis.text.x = element_text(colour="black", size=12),
        axis.ticks = element_line(color="black", size=0.5, lineend = "square"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right")+
  scale_fill_manual(values=my_pal) +
  scale_y_continuous(labels=scientific)
g
ggsave(filename=out_plot_strain_path, plot = g, width = 12, height = 8)








