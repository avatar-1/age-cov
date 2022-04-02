# This script was used to determine transcriptional noise in within each age group for the COVID-19 paper #
# Code developed by Alexander Capraro October 2021

library(tidyverse)
library(magrittr)
library(rstatix)
library(data.table)
library(ggpmisc)

tpm <- read.table(file="tpm.txt",sep = "\t", header = T)

# Mock ####

# 1. Select mock samples
tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B.")) %>%
  select(-contains("V2."))

# 2. Divide genes into 10 equally sized bins based on mean expression and exclude bin 1 and 10
# 3. Calculate co-efficent of variation for all genes

tpm_retained <- tpm %>%
  data.frame(rowMeans(.)) %>%
  set_rownames(row.names(tpm)) %>%
  arrange(rowMeans...) %>% #sort by ascending value
  mutate(bins = ntile(.$rowMeans...,10)) %>% #creates 10 equally sized bins
  filter(bins != 10) %>% #remove top bin
  filter(bins != 1) %>% #remove bottom bin
  mutate(sd = apply(.[,1:26],1,sd)) %>% #calculate standard deviation of genes
  mutate(cv = (sd/rowMeans...)* 100) %>% #calculate co-efficient of variation
  arrange(cv) %>%
  mutate(cv_bins = ntile(.$cv,10)) #split into 10 groups based on co-efficient of variation

low_cv <- tpm_mean %>%
  filter(cv_bins == "1") #only keep 10% of genes with lowest CV

# 4. Calculate euclidean distance between each sample & the sample mean within each age group

tpm_OA <- low_cv %>%
  select(contains("M.OA.")) %>%
  mutate(mean = rowMeans(.))
ED_OA <- as.matrix(dist(t(tpm_OA))) %>%
  data.frame() %>%
  select(mean) %>%
  mutate(AgeGroup = "Older Adult") %>%
  mutate(Condition = "Mock")

tpm_YA <- low_cv %>%
  select(contains("M.YA.")) %>%
  mutate(mean = rowMeans(.))
ED_YA <- as.matrix(dist(t(tpm_YA))) %>%
  data.frame() %>%
  select(mean) %>%
  mutate(AgeGroup = "Young Adult") %>%
  mutate(Condition = "Mock")

tpm_CH <- low_cv %>%
  select(contains("M.CH.")) %>%
  mutate(mean = rowMeans(.))
ED_CH <- as.matrix(dist(t(tpm_CH))) %>%
  data.frame() %>%
  select(mean) %>%
  mutate(AgeGroup = "Child") %>%
  mutate(Condition = "Mock")

ED_mock <- rbind(ED_OA,ED_YA,ED_CH) %>%
  filter(mean > 0) #removes "mean" vs. "mean" distance

write.table(ED_mock, file="ed_mock.txt", sep="\t",quote =F)

# Virus 0.2 MOI ####

# 1. Select mock samples
tpm <- tpm %>%
  select(-contains(".B.")) %>%
  select(contains("V2"))

# 2. Divide genes into 10 equally sized bins based on mean expression and exclude bin 1 and 10
# 3. Calculate co-efficent of variation for all genes

tpm_retained <- tpm %>%
  data.frame(rowMeans(.)) %>%
  set_rownames(row.names(tpm)) %>%
  arrange(rowMeans...) %>% #sort by ascending value
  mutate(bins = ntile(.$rowMeans...,10)) %>% #creates 10 equally sized bins
  filter(bins != 10) %>% #remove top bin
  filter(bins != 1) %>% #remove bottom bin
  mutate(sd = apply(.[,1:29],1,sd)) %>% #calculate standard deviation of genes
  mutate(cv = (sd/rowMeans...)* 100) %>% #calculate co-efficient of variation
  arrange(cv) %>%
  mutate(cv_bins = ntile(.$cv,10)) #split into 10 groups based on co-efficient of variation

low_cv <- tpm_retained %>%
  filter(cv_bins == "1") #only keep 10% of genes with lowest CV

# 4. Calculate euclidean distance between each sample & the sample mean within each age group

tpm_OA <- low_cv %>%
  select(contains("V2.OA.")) %>%
  mutate(mean = rowMeans(.))
ED_OA <- as.matrix(dist(t(tpm_OA))) %>%
  data.frame() %>%
  select(mean) %>%
  mutate(AgeGroup = "Older Adult") %>%
  mutate(Condition = "Virus 0.2")

tpm_YA <- low_cv %>%
  select(contains("V2.YA.")) %>%
  mutate(mean = rowMeans(.))
ED_YA <- as.matrix(dist(t(tpm_YA))) %>%
  data.frame() %>%
  select(mean) %>%
  mutate(AgeGroup = "Young Adult") %>%
  mutate(Condition = "Virus 0.2")

tpm_CH <- low_cv %>%
  select(contains("V2.CH.")) %>%
  mutate(mean = rowMeans(.))
ED_CH <- as.matrix(dist(t(tpm_CH))) %>%
  data.frame() %>%
  select(mean) %>%
  mutate(AgeGroup = "Child") %>%
  mutate(Condition = "Virus 0.2")

ED_v2 <- rbind(ED_OA,ED_YA,ED_CH) %>%
  filter(mean > 0) #removes "mean" vs. "mean" distance

write.table(ED_v2, file="ed_virus0.2.txt", sep="\t",quote =F)

ED_all <- rbind(ED_mock, ED_v2)

write.table(ED_all, file="ed_all.txt", sep="\t",quote =F)


# Plotting ####

ED_all <- read.table(file = "ed_all.txt", sep = "\t", header = T)

stat_test <- ED_all %>%
  group_by(Condition) %>%
  wilcox_test(mean ~ AgeGroup, p.adjust.method = "none")
stat_test

ED_all %>%
  mutate(AgeGroup = factor(AgeGroup, levels = c("Older Adult","Young Adult","Child"))) %>%
  mutate(Condition = factor(Condition, levels = c("Mock","Virus 0.2"))) %>%
  ggplot(aes(x=AgeGroup,y=mean,fill=AgeGroup)) +
  geom_boxplot() +
  geom_point(aes(size=3,color=AgeGroup)) +
  scale_color_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) + 
  scale_fill_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) + 
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Transcriptional noise (Euclidean distance to mean)") +
  xlab("") +
  ylim(5,20) +
  facet_wrap(~Condition)
  
# Delta Euclidean distance (change in transcriptional noise) ####

delta_ed <- read.table(file = "delta_ed.txt", sep = "\t", header = T)

stat_test <- delta_ed %>%
  wilcox_test(Delta_ED ~ AgeGroup, p.adjust.method = "none")
stat_test

delta_ed %>%
  mutate(AgeGroup = factor(AgeGroup, levels = c("Older Adult","Young Adult","Child"))) %>%
  ggplot(aes(x=AgeGroup,y=Delta_ED,fill=AgeGroup)) +
  geom_boxplot() +
  geom_point(aes(size=3,color=AgeGroup)) +
  scale_color_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) + 
  scale_fill_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) + 
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Delta Transcriptional Noise (Euclidean distance to mean)") +
  xlab("")

# Transcriptional Noise Correlation ####

ED_all <- read.table(file = "ed_all.txt", sep = "\t", header = T)

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B."))


TN_Corr <- function(genes){
  input_tpm <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
  input_tpm$gene <- row.names(input_tpm)
  
  input_tpm$gene <- row.names(input_tpm)
  input_tpm <- reshape2::melt(input_tpm)
  input_tpm <- input_tpm[order(input_tpm$gene),] #order by gene
  
  input_tpm$Condition <- ifelse(grepl("V2.", input_tpm$variable), "Virus", "Mock") 
  input_tpm$Age <- ifelse(grepl(".CH", input_tpm$variable), "Child", 
                          ifelse(grepl("YA", input_tpm$variable), "Young Adult","Older Adult"))
  input_tpm$Age <- factor(x = input_tpm$Age, levels = c("Older Adult","Young Adult", "Child")) #Keeps plot in order of child > YA > OA
  input_tpm$gene <- factor(x = input_tpm$gene, levels = genes) #Keeps order of genes
  
  input_tpm <- input_tpm %>%
    filter(Condition == "Virus")
  
  merged <- merge(input_tpm, ED_all, by.x = "variable", by.y = 0)
  merged$Group <- paste0(merged$Condition.x,"_",merged$Age)
  
  merged_dt <- data.table(merged)
  corr_test_all <- merged_dt[ ,cor.test(x=value, y=mean, method = 'spearman')[-2]]
  corr_test <- merged_dt[ ,cor.test(x=value, y=mean, method = 'spearman')[-2], by="gene"]
  
  
  plot(merged %>%
         ggplot(aes(x = mean, y = value)) +
         geom_point(aes(size=3, shape=gene, color=Age, fill=Age)) +
         scale_fill_manual(values =  rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2)) +
         scale_color_manual(values =  rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2)) +
         scale_shape_manual(values=c(15,16,17,23,25,9))+
         theme_bw() +
         theme(legend.position = "none") +
         scale_linetype_discrete(name="Age") +
         geom_smooth(method=lm, formula = y ~ x,fullrange=T) +
         xlab("Transcriptional Noise")+
         ylab(" Gene Expression (log2 TPM)") +
         ylim(-6,10) +
         xlim(5,20) +
         stat_poly_eq(formula = y ~ x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                      parse=TRUE,label.x.npc = "right")) #use to get R^2
  
  print(corr_test_all)
  print(corr_test)
}


#Chemokines
TN_Corr(genes = c("CXCL1", "CXCL2", "CXCL3","CXCL16", "CXCL17","CCL22"))


#Interleukin 1

TN_Corr(genes = c("IL1B"))

#IFN

TN_Corr(genes = c("IFNB1","IFNL1","IFNL2","IFNL3"))

TN_Corr(genes = C5)
