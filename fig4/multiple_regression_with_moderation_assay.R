#! usr/bin/env R

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(stats)
library(coin)
library(exactRankTests)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(broom)
library(MuMIn)
#library(latex2exp)
setwd('/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/assay')

########################read assay data#######################
assay.data.name <- '210924_assay_data.csv'
data <- read_csv(assay.data.name)
data$L452R <- factor(data$L452R, levels = c("0", "1"))
data$T478K <- factor(data$T478K, levels = c("0", "1"))
data$E484K <- factor(data$E484K, levels = c("0", "1"))
data$N501Y <- factor(data$N501Y, levels = c("0", "1"))

data$assay <- ifelse(data$assay =="yeast-binding", "ACE2-binding",
                     ifelse(data$assay =="pseudovirus", "transmissibility",
                            ifelse(data$assay =="fusion", "fusion activity",data$assay)))
mut_list <- colnames(data)[1:4]
mut_set_df <- data %>% select(-assay, -value) %>% unique()

###moderna###
moderna.data.name <- '/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/assay/moderna_NT50_VOC_fold_change.csv'
moderna.data <- read_csv(moderna.data.name)
moderna.data  <- gather(data=moderna.data, -sample, key=mut_number, value = fold_change)
mutset_number <- moderna.data$mut_number %>% unique()
mut_set_df$number <- mutset_number
moderna.data <- moderna.data %>% left_join(mut_set_df, by=c("mut_number"="number")) %>% select(-mut_number)
moderna.data <- rename(moderna.data, value=fold_change)
moderna.data$assay <- "moderna"
moderna.data <- moderna.data %>% select(L452R, T478K, E484K, N501Y,assay, value)
data <- bind_rows(data, moderna.data)

###pfizer data###
pfizer.data.name <- '/Users/yuuka.m/Dropbox/Research/cov-phylo/RBM_project/data/assay/pfizer_NT50_VOC_fold_change.csv'
pfizer.data <- read_csv(pfizer.data.name)
pfizer.data  <- gather(data=pfizer.data, -sample, key=mut_number, value = fold_change)
pfizer.data <- pfizer.data %>% left_join(mut_set_df, by=c("mut_number"="number")) %>% select(-mut_number)
pfizer.data <- rename(pfizer.data, value=fold_change)
pfizer.data$assay <- "pfizer"
pfizer.data <- pfizer.data %>% select(L452R, T478K, E484K, N501Y,assay, value)
data <- bind_rows(data, pfizer.data)

data$value <- ifelse(data$assay=="ACE2-binding", data$value*(-1), data$value)
#assay_list <- data$assay %>% unique()

###normalization###
data <- data %>% group_by(assay) %>% mutate(value_nor=scale(value))
data <- data %>% mutate(assay=ifelse((assay == "moderna" | assay == "pfizer"),
                                     "neutralization", assay))
assay_list <- data$assay %>% unique()
############################################################


###multiple regression with moderation between two variants###
r2_array <- c()
BIC_array <- c()
res_lm_m_df <- tibble(term=character(), estimate=numeric(), std.error=numeric(), statistic=numeric(), p.value=numeric(), assay=character())
for (assay_name in assay_list){
  test_df <- filter(data, assay == assay_name) %>% as.data.frame()
  n_data <- nrow(test_df)
  #test_df$value <-  test_df$value - ave(test_df$value)
  res <- lm(value_nor~L452R+T478K+E484K+N501Y+L452R:T478K+L452R:E484K+L452R:N501Y+T478K:E484K+T478K:N501Y+E484K:N501Y,
            data = test_df)
  res_lm_m_df_tmp <- tidy(res)
  res_lm_m_df_tmp$assay <- assay_name
  res_lm_m_df <- bind_rows(res_lm_m_df, res_lm_m_df_tmp)
  print(summary(res))
  r.squared <- summary(res)$adj.r.squared
  BIC_array <-c(BIC_array,  extractAIC(res,k=log(n_data))[2])
  r2_array <- c(r2_array, r.squared)
  # f.stat <- summary(res)$fstatistic
  # p.value <- 1-pf(f.stat["value"],f.stat["numdf"],f.stat["dendf"])[[1]]
  # pval_array <- c(pval_array, p.value)
  print(step(res, k=log(n_data)))
}
BIC_array
res_lm_m_df
res_lm_tbl <- tibble(assay = assay_list, r_squares_adj = r2_array, BIC=BIC_array)
res_lm_tbl
res_lm_m_df <- left_join(res_lm_m_df, res_lm_tbl)
write_tsv(res_lm_m_df, file = 'res_multiple_regression_with_moderation_between_two_variables.txt')


heatmap.col.pal <-  brewer.pal(5, "RdBu")
heatmap.col.pal <-  brewer.pal(5, "PiYG")
###Heatmap###
plots <- list()
for (i in 1:4){
  assay_name <- assay_list[i]
  test_df <- filter(res_lm_m_df, assay == assay_name)
  test_df <- select(test_df, term, estimate, p.value) %>% tail(10)
  test_df$var1 <- substr(test_df$term, 1, 5)
  test_df$var2 <- str_sub(test_df$term, -6, -2)
  test_df  <- select(test_df, -term)
  test_df$var1 <- factor(test_df$var1, levels=c("L452R", "T478K", "E484K", "N501Y"))
  test_df$var2 <- factor(test_df$var2, levels=c("L452R", "T478K", "E484K", "N501Y"))
  test_df$estimate <- round(test_df$estimate, digits = 2)
  test_df$pval_code <- ifelse(test_df$p.value < 0.05, "*", "")
  test_df$pval_code <- ifelse(test_df$p.value < 0.01, "**", test_df$pval_code)
  test_df$pval_code <- ifelse(test_df$p.value < 0.001, "***", test_df$pval_code)
  plots[[i]] <- test_df %>% ggplot() +
    geom_tile(aes(x=var1, y=var2, fill=estimate))+
    geom_text(aes(x=var1, y=var2, label=estimate))+
    geom_text(aes(x=var1, y=var2, label=pval_code), position=position_nudge(y = -0.2))+
    theme_minimal()+
    #scale_fill_gradient2(low="#0068b7", mid="white", high="orange", midpoint=0)+
    scale_fill_gradient2(low=heatmap.col.pal[5], mid=heatmap.col.pal[3], high=heatmap.col.pal[1], midpoint=0)+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          # axis.line = element_blank(),
          axis.line = element_line(color="black", size=0.5, lineend = "square"),
          axis.ticks = element_line(color="black", size=0.5, lineend = "square"),
          axis.text.x.bottom = element_text(colour = "black", size = 12, angle=60, hjust = 0, vjust= 0),
          axis.text.y.left = element_text(colour = "black", size = 12, angle=60,hjust = 0.7, vjust= 0),
          axis.title = element_blank(),
          legend.title = element_text(colour = "black", size = 12),
          legend.text = element_text(colour = "black", size = 12))+
    labs(title = assay_name,
         subtitle = bquote("R"^2==.(round(res_lm_tbl[[i, 'r_squares_adj']],3))~~"BIC"==.(round(res_lm_tbl[[i, 'BIC']],3)) ) )
}
#plots

g <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]]
g <- g + plot_layout(ncol = 4)
g
# out.name <- "heatmap_estimate_coefficient_regression_each_assay.pdf"
out.name <- "heatmap_estimate_coefficient_regression_each_assay_normalized.pdf"
ggsave(filename = out.name, plot = g, width = 16, height = 4)


###compare###
compare_df <- tibble(assay=assay_list)
compare_df$all_BIC <- 0
compare_df$partial_BIC <- 0
compare_df$no_BIC <- 0
for (i in 1:length(assay_list)){
  assay_name <- assay_list[i]
  test_df <- filter(data, assay == assay_name) %>% as.data.frame()
  n_data <- nrow(test_df)
  test_df$value <-  test_df$value - ave(test_df$value)
  res_n <- lm(value~L452R+T478K+E484K+N501Y, data = test_df)
  res_p <- lm(value~L452R+T478K+E484K+N501Y+L452R:T478K+L452R:E484K+L452R:N501Y+T478K:E484K+T478K:N501Y+E484K:N501Y,
            data = test_df)
  res_a <- lm(value~L452R*T478K*E484K*N501Y, data = test_df)
  compare_df$partial_BIC[i] <- extractAIC(res_p, k=log(n_data))[2]
  compare_df$all_BIC[i] <- extractAIC(res_a, k=log(n_data))[2]
  compare_df$no_BIC[i] <- extractAIC(res_n, k=log(n_data))[2]
}
compare_df
write_tsv(compare_df, file = 'res_BIC_multiple_regression_compare.txt')


####with no moderation####
r2_array <- c()
BIC_array <- c()
res_lm_df <- tibble(term=character(), estimate=numeric(), std.error=numeric(), statistic=numeric(), p.value=numeric(), assay=character())
for (assay_name in assay_list){
  test_df <- filter(data, assay == assay_name) %>% as.data.frame()
  n_data <- nrow(test_df)
  #test_df$value <-  test_df$value - ave(test_df$value)
  res <- lm(value_nor~L452R+T478K+E484K+N501Y,
            data = test_df)
  res_lm_df_tmp <- tidy(res)
  res_lm_df_tmp$assay <- assay_name
  res_lm_df <- bind_rows(res_lm_df, res_lm_df_tmp)
  print(summary(res))
  r.squared <- summary(res)$adj.r.squared
  BIC_array <-c(BIC_array,  extractAIC(res,k=log(n_data))[2])
  r2_array <- c(r2_array, r.squared)
  # f.stat <- summary(res)$fstatistic
  # p.value <- 1-pf(f.stat["value"],f.stat["numdf"],f.stat["dendf"])[[1]]
  # pval_array <- c(pval_array, p.value)
  print(step(res, k=log(n_data)))
}
BIC_array
res_lm_df

res_lm_df <- filter(res_lm_df, term!="(Intercept)")
res_lm_df$estimate <- round(res_lm_df$estimate,digits = 2)
res_lm_df$term <- substr(res_lm_df$term, 1, 5)
res_lm_df$term <- factor(res_lm_df$term, levels=c("L452R", "T478K", "E484K", "N501Y"))
res_lm_df$pval_code <- ifelse(res_lm_df$p.value < 0.05, "*", "")
res_lm_df$pval_code <- ifelse(res_lm_df$p.value < 0.01, "**", res_lm_df$pval_code)
res_lm_df$pval_code <- ifelse(res_lm_df$p.value < 0.001, "***", res_lm_df$pval_code)

g <- ggplot(res_lm_df)+
  geom_tile(aes(x=term, y=assay, fill=estimate))+
  geom_text(aes(x=term, y=assay, label=estimate))+
  geom_text(aes(x=term, y=assay, label=pval_code), position=position_nudge(y = -0.2))+
  theme_minimal()+
  #scale_fill_gradient2(low="#0068b7", mid="white", high="orange", midpoint=0)+
  scale_fill_gradient2(low=heatmap.col.pal[5], mid=heatmap.col.pal[3], high=heatmap.col.pal[1], midpoint=0)+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(colour = "black", size = 12),
        legend.text = element_text(colour = "black", size = 12),
        axis.line = element_line(color="black", size=0.5, lineend = "square"),
        axis.ticks = element_line(color="black", size=0.5, lineend = "square"),
        axis.title = element_text(colour = "black", size = 12),
        axis.text.x.bottom = element_text(colour = "black", size = 12, angle=60,hjust = 0, vjust= 0),
        axis.text.y.left = element_text(colour = "black", size = 12, angle=60,hjust = 0.7, vjust= 0))+
  labs(x="Mut", y="Assay")
g

# out.name <- "heatmap_estimate_coefficient_regression_each_assay.with_no_morderation.pdf"
out.name <- "heatmap_estimate_coefficient_regression_each_assay.with_no_morderation_normalized.pdf"
ggsave(filename = out.name, plot = g, width = 5, height = 6)









# ###anova###
# res_anova_df <- tibble(term=character(),
#                  sumsq=numeric(),
#                  meansq=numeric(),
#                  statistic=numeric(),
#                  p.value=numeric(),
#                  assay=character())
# for (assay_name in assay_list){
#   print(assay_name)
#   res_anova_df <- filter(data, assay == assay_name) %>% as.data.frame()
#   res_anova_df$value <-  res_anova_df$value - ave(res_anova_df$value)
#   res <- aov(value~L452R+T478K+E484K+N501Y+L452R:T478K+L452R:E484K+L452R:N501Y+T478K:E484K+T478K:N501Y+E484K:N501Y, data = res_anova_df)
#   print(summary(res))
#   res_anova_df_tmp <- tidy(res)
#   res_anova_df_tmp$assay <- assay_name
#   res_anova_df_tmp <- select(res_anova_df_tmp, -df)
#   res_anova_df <- bind_rows(res_anova_df, res_anova_df_tmp)
# }
# res_anova_df  <- res_anova_df %>% select(assay, everything())
# write_tsv(res_anova_df, file = 'res_anova_with_moderation_between_two_variables.txt')
#
#
# res_anova_all_df <- tibble(term=character(),
#                  sumsq=numeric(),
#                  meansq=numeric(),
#                  statistic=numeric(),
#                  p.value=numeric(),
#                  assay=character())
# for (assay_name in assay_list){
#   print(assay_name)
#   test_df <- filter(data, assay == assay_name) %>% as.data.frame()
#   test_df$value <-  test_df$value - ave(test_df$value)
#   res <- aov(value~L452R*T478K*E484K*N501Y, data = test_df)
#   print(summary(res))
#   res_anova_all_df_tmp <- tidy(res)
#   res_anova_all_df_tmp$assay <- assay_name
#   res_anova_all_df_tmp <- select(res_anova_all_df_tmp, -df)
#   res_anova_all_df <- bind_rows(res_anova_all_df, res_anova_all_df_tmp)
# }
# res_anova_all_df  <- res_anova_all_df %>% select(assay, everything())
# res_anova_all_df
#
# ######simple main effect######
# ####yeast-binding###
# ##T478K+N501Y##
# assay_name <- assay_list[1]
# test_df_T478K <- filter(data, assay == assay_name, T478K==1) %>% as.data.frame()
# t.test(test_df_T478K$value~test_df_T478K$N501Y, var.equal=T)
# test_df_T478K <- filter(data, assay == assay_name, T478K==0) %>% as.data.frame()
# t.test(test_df_T478K$value~test_df_T478K$N501Y, var.equal=T)
#
# test_df_N501Y <- filter(data, assay == assay_name, N501Y==1) %>% as.data.frame()
# t.test(test_df_N501Y$value~test_df_N501Y$T478K, var.equal=T)
# test_df_N501Y <- filter(data, assay == assay_name, N501Y==0) %>% as.data.frame()
# t.test(test_df_N501Y$value~test_df_N501Y$T478K, var.equal=T)
#
# ##E484K+N501Y##
# test_df_E484K <- filter(data, assay == assay_name, E484K==1) %>% as.data.frame()
# t.test(test_df_E484K$value~test_df_E484K$N501Y, var.equal=T)
# test_df_E484K <- filter(data, assay == assay_name, E484K==0) %>% as.data.frame()
# t.test(test_df_E484K$value~test_df_E484K$N501Y, var.equal=T)
#
# test_df_N501Y <- filter(data, assay == assay_name, N501Y==1) %>% as.data.frame()
# t.test(test_df_N501Y$value~test_df_N501Y$E484K, var.equal=T)
# test_df_N501Y <- filter(data, assay == assay_name, N501Y==0) %>% as.data.frame()
# t.test(test_df_N501Y$value~test_df_N501Y$E484K, var.equal=T)
#
# ####pseudovirus###
# ##T478K+E484K##
# assay_name <- assay_list[2]
# test_df_T478K <- filter(data, assay == assay_name, T478K==1) %>% as.data.frame()
# t.test(test_df_T478K$value~test_df_T478K$E484K, var.equal=T)
# test_df_T478K <- filter(data, assay == assay_name, T478K==0) %>% as.data.frame()
# t.test(test_df_T478K$value~test_df_T478K$E484K, var.equal=T)
#
# test_df_E484K <- filter(data, assay == assay_name, E484K==1) %>% as.data.frame()
# t.test(test_df_E484K$value~test_df_E484K$T478K, var.equal=T)
# test_df_E484K <- filter(data, assay == assay_name, E484K==0) %>% as.data.frame()
# t.test(test_df_E484K$value~test_df_E484K$T478K, var.equal=T)
#
#
#
#
# ###between two factors###
# assay_name <- assay_list[1]
# test_df <- filter(data, assay == assay_name) %>% as.data.frame()
# test_df$value <-  test_df$value - ave(test_df$value)
# res <- aov(value~L452R*T478K, data = test_df)
# summary(res)
# interaction.plot(test_df$L452R, test_df$T478K, test_df$value)
#
# res <- aov(value~L452R*E484K, data = test_df)
# summary(res)
# interaction.plot(test_df$L452R, test_df$E484K, test_df$value)
#
# summary(aov(value~T478K*E484K, data = test_df))
# interaction.plot(test_df$T478K, test_df$E484K, test_df$value)
#
# summary(aov(value~T478K*N501Y, data = test_df))
# interaction.plot(test_df$T478K, test_df$N501Y, test_df$value)
#
# summary(aov(value~E484K*N501Y, data = test_df))
# interaction.plot(test_df$E484K, test_df$N501Y, test_df$value)
