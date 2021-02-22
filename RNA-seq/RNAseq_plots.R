#Code is used for analysis of bronchial and nasal RNA-seq data for COVID-19 2020 project
#Code developed by Alexander Capraro, December 2020
library(tidyverse)
library(edgeR)
library(plotly)
library(gprofiler2)
library(gplots)
library(RUVSeq)
library(Hmisc)
library(pvclust)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(viridis)
library(EnhancedVolcano)
library(rmulti)
library(rstatix)
library(dplyr)
library(ggpubr)

#### Gene Plots ####

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
  
  ggplot(zscore_long, aes(x = Condition, y = zscore, fill = Condition)) + 
    geom_boxplot() +
    theme_minimal() +
    scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D","#BEBADA","#FFFFB3","#8DD3C7","#F7941D")) +
    ylim(-3,4)
  #export 4.5x6
}

#Viral genes

heatmap_average(genes = c("GU280_gp01", "GU280_gp02", "GU280_gp03", "GU280_gp04", "GU280_gp05", "GU280_gp06", 
                          "GU280_gp07", "GU280_gp08", "GU280_gp09", "GU280_gp10", "GU280_gp11"))

#Viral read percentage

viral_reads <- read.table(file = "viral_reads.txt", sep = "\t", header = T)
viral_reads$Condition <- factor(viral_reads$Condition, levels = c("OA","YA","CH","HBE"))

stat_test <- viral_reads %>% 
  wilcox_test(Percentage ~ Condition, p.adjust.method = "bonferroni")
print(stat_test)

ggplot(viral_reads, aes(x = Condition, y = Percentage, fill = Condition)) +
        geom_boxplot() +
        theme_minimal() +
        scale_fill_brewer(palette = "Set3")

#Viral entry factors

heatmap_average(genes = c("ACE","ACE2", "TMPRSS2","FURIN","CTSL", "ANPEP","DPP4","TMPRSS4","NRP1","ATP1A1","ATP1B1"))
boxplot_zscore(genes = c("ACE","ACE2", "TMPRSS2","FURIN","CTSL", "ANPEP","DPP4","TMPRSS4","NRP1","ATP1A1","ATP1B1"))

#IFN + IRF + IFNR

heatmap_average(genes = c("IFNB1", "IFNL1","IFNL2",
                          "IRF1","IRF2","IRF3","IRF5","IRF6","IRF7","IRF9",
                          "IFNAR1","IFNAR2","IFNGR1","IFNGR2","IFNLR1"))

#ISGs 
interferon_genes <- read.table(file = "interferon_genes.txt", sep = "\t", header = T)
ISG_rna_processing <- read.table(file = "ISG_rna_processing.txt", sep = "\t", header = T)
ISG_DEGs <- read.table(file = "ISG_DEGs_All.txt", sep = "\t", header = T)
ISG_Down <- read.table(file = "ISG_Down.txt", sep = "\t", header = T)

heatmap_average(genes = interferon_genes$Gene[interferon_genes$Cluster == "2"]) #Change number for different clusters
boxplot_zscore(genes = interferon_genes$Gene[interferon_genes$Cluster == "5"])

heatmap_average(genes = ISG_rna_processing$Gene) #Combined C1 and C2
boxplot_zscore(genes = ISG_rna_processing$Gene)

heatmap_average(genes = ISG_DEGs$Gene[ISG_DEGs$Cluster == "3"])

heatmap_average(genes = ISG_Down$Gene)

#Chemokines

heatmap_average(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                          "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                          "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                          "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"))

boxplot_zscore(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                         "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                         "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                         "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"))

#Interleukins

heatmap_average(genes = c("IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                          "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"))

boxplot_zscore(genes = c("IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                         "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"))

#Cytokines
heatmap_average(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                          "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                          "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                          "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17",
                          "IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                          "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"))
boxplot_zscore(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                         "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                         "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                         "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17",
                         "IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                         "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"))

#TLRs
boxplot_tpm(genes = c("TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6"))
heatmap_average(genes = c("TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6","DDX58","IFIH1","DHX58"))
boxplot_zscore(genes = c("TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6","DDX58","IFIH1","DHX58"))

#TLR4 associated genes

boxplot_tpm(genes = c("IRAK1","IRAK2","IRAK4","TIRAP","TRAM1","TRAF3","TRAF6","NFKBIA","NFKBIB","NFKBIE","TBK1"))
heatmap_average(genes = c("IRAK1","IRAK2","IRAK4","TIRAP","TRAM1","TRAF3","TRAF6","NFKBIA","NFKBIB","NFKBIE","TBK1"))
boxplot_zscore(genes = c("IRAK1","IRAK2","IRAK4","TIRAP","TRAM1","TRAF3","TRAF6","NFKBIA","NFKBIB","NFKBIE","TBK1"))

#RIG-I-like receptors

boxplot_tpm(genes = c("DDX58","IFIH1","DHX58"))
heatmap_average(genes = c("DDX58","IFIH1","DHX58"))
boxplot_zscore(genes = c("DDX58","IFIH1","DHX58"))


boxplot_tpm(genes = c("OAS1","OAS2","IFIH1","RNASEL","EIF2AK2"))

#### IFN Correlation ####

#### Viral load ####
tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm[,c(1:16,18:40,42:44)] #Remove BC_B1
tpm_plots <- tpm #keep all

genes = c("IFNB1","IFNL1","IFNL2")

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

viral_reads <- read.table(file = "viral_tpm.txt", sep = "\t", header = T)

mean_ratio <- merge(input_tpm, viral_reads, by.x = "variable", by.y = "ID")
mean_ratio$Condition <- paste(mean_ratio$Condition, mean_ratio$Age, sep="_")


corr_test <- cor.test(x=mean_ratio$value, y=mean_ratio$virus, method = 'spearman')
corr_test

ggplot(mean_ratio, aes(x = virus, y = value, color = Condition)) +
  geom_point() +
  scale_color_manual(values = c("#F7941D","#8DD3C7","#BEBADA","#FFFFB3","#F7941D","#8DD3C7","#BEBADA","#FFFFB3")) +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_smooth(method=lm, color="red", fill = "red", formula = y~x) +
  ylim(-4,2.5)#4.5x4
   
tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm[,c(1:16,18:40,42:44)] #Remove BC_B1
tpm_plots <- tpm #keep all

#Virus age correlation

mean_ratio <- merge(input_tpm, viral_reads, by.x = "variable", by.y = "ID")
mean_ratio$Group <- paste(mean_ratio$Condition, mean_ratio$Age, sep="_")

mean_ratio <- mean_ratio[mean_ratio$Condition == "Virus",]

corr_test <- cor.test(x=mean_ratio$value, y=mean_ratio$age, method = 'spearman')
corr_test

ggplot(mean_ratio, aes(x = age, y = value, color = Group, shape = gene)) +
  geom_point() +
  scale_color_manual(values = c("#F7941D","#8DD3C7","#BEBADA","#FFFFB3","#F7941D","#8DD3C7","#BEBADA","#FFFFB3")) +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_smooth(method=lm, color="red", fill = "red", formula = y~x) +
  ylim(-4,2.5)#4.5x4

