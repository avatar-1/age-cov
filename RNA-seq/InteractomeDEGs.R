# This script was to analyse the expression of SARS-CoV-2 interacting genes for the COVID-19 paper #
# Code developed by Alexander Capraro, January 2021

#### SARS-CoV-2 interactome gene expression ####

####MOCK ####
tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm[,c(1:16,18:40,42:44)] #Remove BC_B1
tpm_plots <- tpm[,grepl('Mock', colnames(tpm))] 
#tpm_plots <- tpm[,grepl('Virus', colnames(tpm))] 

interactome_genes <- read.table(file = "interactome_genes.txt", sep = "\t", header = T)
agegroups <- read.table(file = "boxplot_zscore_groups_All.txt", sep = "\t", header = T)

interactome_tpm <- tpm_plots[interactome_genes$Protein, ]
interactome_tpm <- interactome_tpm[complete.cases(interactome_tpm),] #Remove NAs
zscore <- data.frame(scale(t(interactome_tpm)))
zscore <- data.frame(t(zscore))
zscore$geneID <- row.names(zscore)

zscore <- zscore %>% gather(ID, zscore, Mock_CH_B_376AC_B2_t1:Mock_YA_N_SB_B2_t1)

zscore_long <- merge(zscore, agegroups, by = "ID")
zscore_long <- merge(zscore_long, interactome_genes,  by.x = "geneID", by.y= "Protein")
zscore_long$Condition <- factor(zscore_long$Condition, levels = c("OA_Mock","YA_Mock","CH_Mock","B_Mock"))

stat_test <- zscore_long %>% 
  group_by(Type) %>%
  wilcox_test(zscore ~ Condition, p.adjust.method = "bonferroni")
stat_test

ggplot(zscore_long, aes(x = Condition, y = zscore, fill = Condition)) + 
  geom_boxplot() +
  theme_minimal() +
  ylim(-4,4) +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D")) +
  facet_wrap(~Type, ncol = 5,
             scales = "fixed")

####VIRUS #### 

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm[,c(1:16,18:40,42:44)] #Remove BC_B1
tpm_plots <- tpm[,grepl('Virus', colnames(tpm))] 

interactome_genes <- read.table(file = "interactome_genes.txt", sep = "\t", header = T)
agegroups <- read.table(file = "boxplot_zscore_groups_All.txt", sep = "\t", header = T)

interactome_tpm <- tpm_plots[interactome_genes$Protein, ]
interactome_tpm <- interactome_tpm[complete.cases(interactome_tpm),] #Remove NAs
zscore <- data.frame(scale(t(interactome_tpm)))
zscore <- data.frame(t(zscore))
zscore$geneID <- row.names(zscore)

zscore <- zscore %>% gather(ID, zscore, Virus_CH_B_376AC_B2_t1:Virus_YA_N_SB_B2_t1)

zscore_long <- merge(zscore, agegroups, by = "ID")
zscore_long <- merge(zscore_long, interactome_genes,  by.x = "geneID", by.y= "Protein")
zscore_long$Condition <- factor(zscore_long$Condition, levels = c("OA_Virus","YA_Virus","CH_Virus","B_Virus"))
my_comparisons <- list(c("CH", "OA"), c("CH", "YA"), c("YA", "OA"))

stat_test <- zscore_long %>% 
  group_by(Type) %>%
  wilcox_test(zscore ~ Condition, p.adjust.method = "bonferroni")
stat_test

ggplot(zscore_long, aes(x = Condition, y = zscore, fill = Condition)) + 
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D")) +
  ylim(-4,4) +
  facet_wrap(~Type, ncol = 5,
             scales = "fixed")

#### Combined ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm_plots <- tpm[,c(1:16,18:40,42:44)] #Remove BC_B1

interactome_genes <- read.table(file = "interactome_genes.txt", sep = "\t", header = T)
agegroups <- read.table(file = "boxplot_zscore_groups_All.txt", sep = "\t", header = T)

interactome_tpm <- tpm_plots[interactome_genes$Protein, ]
interactome_tpm <- interactome_tpm[complete.cases(interactome_tpm),] #Remove NAs
zscore <- data.frame(scale(t(interactome_tpm)))
zscore <- data.frame(t(zscore))
zscore$geneID <- row.names(zscore)

zscore <- zscore %>% gather(ID, zscore, Mock_CH_B_376AC_B2_t1:Virus_YA_N_SB_B2_t1)

zscore_long <- merge(zscore, agegroups, by = "ID")
zscore_long <- merge(zscore_long, interactome_genes,  by.x = "geneID", by.y= "Protein")
zscore_long$Condition <- factor(zscore_long$Condition, levels = c("OA_Mock","YA_Mock","CH_Mock","B_Mock",
                                                                  "OA_Virus","YA_Virus","CH_Virus","B_Virus"))

stat_test <- zscore_long %>% 
  group_by(Type) %>%
  wilcox_test(zscore ~ Condition, p.adjust.method = "bonferroni")
stat_test

ggplot(zscore_long, aes(x = Condition, y = zscore, fill = Condition)) + 
  geom_boxplot() +
  theme_minimal() +
  ylim(-4,4) +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D","#BEBADA","#FFFFB3","#8DD3C7","#F7941D")) +
  facet_wrap(~CoV2, ncol = 5,
             scales = "fixed")

rep <- subset(zscore_long, Type %in% c("Replication complex"))

stat_test <- rep %>%
  group_by(CoV2) %>%
  wilcox_test(zscore ~ Condition, p.adjust.method = "bonferroni")
stat_test

ggplot(rep, aes(x = Condition, y = zscore, fill = Condition)) + 
  geom_violin() +
  theme_minimal() +
  ylim(-4,4) +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D","#BEBADA","#FFFFB3","#8DD3C7","#F7941D"))  +
  facet_wrap(~CoV2, ncol = 3,
             scales = "fixed") +
  stat_summary(fun.y=mean, geom="point", shape=16, size=2, position = position_dodge(0.9))

cap <- subset(zscore_long, Type %in% c("Capping enzymes"))

stat_test <- cap %>%
  group_by(CoV2) %>%
  wilcox_test(zscore ~ Condition, p.adjust.method = "bonferroni")
stat_test

ggplot(cap, aes(x = Condition, y = zscore, fill = Condition)) + 
  geom_violin() +
  theme_minimal() +
  ylim(-4,4) +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D","#BEBADA","#FFFFB3","#8DD3C7","#F7941D"))  +
  facet_wrap(~CoV2, ncol = 5,
             scales = "fixed") +
  stat_summary(fun.y=mean, geom="point", shape=16, size=2, position = position_dodge(0.9))

nsp1 <- subset(zscore_long, CoV2 %in% c("nsp1"))

stat_test <- nsp1 %>%
  wilcox_test(zscore ~ Condition, p.adjust.method = "bonferroni")
stat_test

ggplot(nsp1, aes(x = Condition, y = zscore, fill = Condition)) + 
  geom_violin() +
  theme_minimal() +
  ylim(-4,4) +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D","#BEBADA","#FFFFB3","#8DD3C7","#F7941D"))  +
  facet_wrap(~CoV2, ncol = 3,
             scales = "fixed") +
  stat_summary(fun.y=mean, geom="point", shape=16, size=2, position = position_dodge(0.9))

#### Interactome replication complex reads ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm_plots <- tpm[,c(1:16,18:40,42:44)] #Remove BC_B1

interactome_genes <- read.table(file = "interactome_genes.txt", sep = "\t", header = T)
genes <- interactome_genes[interactome_genes$CoV2 == "nsp12",]
nsp12 <- data.frame((subset(tpm_plots, rownames(tpm_plots) %in% genes$Protein)))
write.table(nsp12, file = "/Users/z3416833/Downloads/nsp12.txt", sep = "\t", quote = F)

#### Interactome - viral load correlation ####

corr <- read.table(file = "rep_vl_correlation_VIRUS.txt", sep = "\t", header = T)

corr_test <- cor.test(x=corr$nsp8, y=corr$virus, method = 'spearman')
corr_test

ggplot(corr, aes(x = virus, y = nsp8)) +
  geom_point(color = "#00AFBB") +
  theme_minimal() +
  geom_smooth(method=lm, color="#00AFBB", fill = "#00AFBB") +
  ylim(4,6) + xlim(5,15) #4.5x4

corr_test <- cor.test(x=corr$nsp12, y=corr$virus, method = 'spearman')
corr_test 

ggplot(corr, aes(x = virus, y = nsp12)) +
  geom_point(color = "#E7B800") +
  theme_minimal() +
  geom_smooth(method=lm, color="#E7B800", fill = "#E7B800") +
  ylim(4,6) + xlim(5,15) #4.5x4

#### NSP RNA interactome ####

library(tidyverse)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(viridis)
library(rmulti)
library(rstatix)
library(dplyr)
library(ggpubr)

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm[,c(1:16,18:40,42:44)] #Remove BC_B1
tpm_plots <- tpm #keep all

boxplot_tpm <- function(genes, scales, cols){
  
  input_tpm <- data.frame((subset(tpm_plots, rownames(tpm_plots) %in% genes))) #Subset only genes in list
  input_tpm$gene <- row.names(input_tpm)
  input_tpm <- reshape2::melt(input_tpm)
  input_tpm <- input_tpm[order(input_tpm$gene),] #order by gene
  
  input_tpm$Condition <- ifelse(grepl("Mock", input_tpm$variable), "Mock", "Virus") 
  input_tpm$Age <- ifelse(grepl("CH_N", input_tpm$variable), "Child", 
                          ifelse(grepl("YA", input_tpm$variable), "Young Adult", 
                                 ifelse(grepl("OA", input_tpm$variable), "Older Adult", "Bronchial")))
  input_tpm$Age <- factor(x = input_tpm$Age, levels = c("Bronchial","Child", "Young Adult", "Older Adult")) #Keeps plot in order of child > YA > OA
  input_tpm$gene <- factor(x = input_tpm$gene, levels = genes) #Keeps order of genes
  
  #Plot individual samples
  ggplot(input_tpm, aes (x = Age, y = value, fill = Condition)) + 
    geom_boxplot(alpha = 1) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3")+ 
    theme(legend.position = "none") +
    facet_wrap(~gene, ncol = 3,
               scales = "free") +
    ylab("log2(TPM)")
}
heatmap_tpm <- function(genes){
  heatmap_plot <- data.frame((subset(tpm_plots, rownames(tpm_plots) %in% genes))) #Subset only genes in list
  heatmap.2(as.matrix(heatmap_plot), scale = "row", trace = 'none', density.info = 'none', col = viridis)
  
} #Plots heatmap version of boxplot_tpm
heatmap_average <- function(genes){
  heatmap_plot <- data.frame((subset(tpm_plots, rownames(tpm_plots) %in% genes))) #Subset only genes in list
  average_plot <- data.frame(row.names = row.names(heatmap_plot))
  average_plot$Mock_OA <- rowMeans(heatmap_plot[,11:15])
  average_plot$Mock_YA <- rowMeans(heatmap_plot[,16:19])
  average_plot$Mock_CH <- rowMeans(heatmap_plot[,6:10])
  average_plot$Mock_Bronchial <- rowMeans(heatmap_plot[,1:5])
  average_plot$Virus_OA <- rowMeans(heatmap_plot[,33:38])
  average_plot$Virus_YA <- rowMeans(heatmap_plot[,39:42])
  average_plot$Virus_CH <- rowMeans(heatmap_plot[,25:32])
  average_plot$Virus_Bronchial <- rowMeans(heatmap_plot[,20:24])
  Mean <- summarise_all(average_plot, mean)
  average_plot <- rbind(average_plot, Mean)
  rownames(average_plot)[rownames(average_plot) == "1"] <- "Mean"
  heatmap.2(as.matrix(average_plot), dendrogram='none', Rowv=row.names(average_plot), breaks=seq(-1.5,1.5,0.2), Colv=FALSE, scale = "row", trace = 'none', density.info = 'none', col = viridis)
} #Heatmap of average for each condition + age group 

boxplot_zscore <- function(genes, scales, cols){
  groups <- read.table(file = "boxplot_zscore_groups_All.txt", sep = "\t", header = T)
  input_tpm <- data.frame((subset(tpm_plots, rownames(tpm_plots) %in% genes))) #Subset only genes in list
  zscore <- data.frame(scale(t(input_tpm)))
  zscore <- data.frame(t(zscore))
  zscore$geneID <- row.names(zscore)
  zscore_long <- zscore %>% gather(ID, zscore, Mock_CH_B_376AC_B2_t1:Virus_YA_N_SB_B2_t1)
  zscore_long <- merge(zscore_long, groups, by = "ID")
  zscore_long$Condition <- factor(zscore_long$Condition, levels = c("OA_Mock","YA_Mock","CH_Mock","B_Mock",
                                                                    "OA_Virus","YA_Virus","CH_Virus","B_Virus"))
  
  .GlobalEnv$stat_test <- zscore_long %>% 
    wilcox_test(zscore ~ Condition, p.adjust.method = "bonferroni")
  print(stat_test)
  
  plot<-  ggplot(zscore_long, aes(x = Condition, y = zscore, fill = Condition)) + 
    geom_violin() +
    theme_minimal() +
    ylim(-4,4) +
    scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D","#BEBADA","#FFFFB3","#8DD3C7","#F7941D"))  +
    stat_summary(fun.y=mean, geom="point", shape=16, size=2, position = position_dodge(0.9))
  #export 4.5x6
  return(list(plot=plot, stat_test = stat_test))
}

# NSP1 - rRNA 18s 

boxplot_tpm(genes = c("RN7SL1"))
boxplot_tpm(genes = c("RNA18SN1", "RNA18SN2", "RNA18SN3", "RNA18SN4", "RNA18SN5"))
boxplot_zscore(genes = c("RNA18SN1", "RNA18SN2", "RNA18SN3", "RNA18SN4", "RNA18SN5"))

heatmap_average((genes = c("RNA18SN1", "RNA18SN2", "RNA18SN3", "RNA18SN4", "RNA18SN5")))
