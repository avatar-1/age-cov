#Code is used for plotting figures of RNA-seq analysis for COVID-19 2020 project
#Code developed by Alexander Capraro, August 2021
library(tidyverse)
library(gplots)
library(ggplot2)
library(gridExtra)
library(viridis)
library(ComplexHeatmap)
library(RColorBrewer)
library(rstatix)
library(ggpmisc)
library(data)
library(magrittr)
library(data.table)
library(edgeR)
library(vegan)
library(circlize)
library(EnhancedVolcano)


tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)

# Gene Plots Functions ####

boxplot_tpm <- function(genes, scales, cols){
  
  input_tpm <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
  input_tpm$gene <- row.names(input_tpm)
  input_tpm <- input_tpm
  
  input_tpm <- reshape2::melt(input_tpm) #wide-form to long-form
  input_tpm <- input_tpm[order(input_tpm$gene),] #order by gene
  
  input_tpm$Condition <- ifelse(grepl("V2.", input_tpm$variable), "Virus 0.2", 
                                ifelse(grepl("Nafa.", input_tpm$variable), "Nafa",
                                       ifelse(grepl("V6.", input_tpm$variable), "Virus 0.6","Mock")))
  input_tpm$Tissue <- ifelse(grepl(".N.", input_tpm$variable), "HNE", "HBE")
  input_tpm$Age <- ifelse(grepl(".CH", input_tpm$variable), "Child", 
                          ifelse(grepl("YA", input_tpm$variable), "Young Adult","Older Adult"))
  input_tpm$Age <- factor(x = input_tpm$Age, levels = c("Older Adult","Young Adult", "Child")) #Keeps plot in order of child > YA > OA
  input_tpm$gene <- factor(x = input_tpm$gene, levels = genes) #Keeps order of genes
  
  input_tpm <- input_tpm %>%
    filter(Age == "Child")  %>%
    filter(Condition == "Mock" | Condition == "Virus 0.2")
  
  #Plot individual samples
  ggplot(input_tpm, aes (x = Condition, y = value, fill = Tissue)) + 
    geom_boxplot(alpha = 1) +
    geom_point(aes(fill = Age), size = 2, shape = 21, position = position_jitterdodge()) +
    theme_bw() +
    scale_fill_manual(values=rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2))+ 
    theme(legend.position = "none") +
    facet_wrap(~gene, ncol = 4,
               scales = "free") +
    ylab("log2(TPM)")
}
heatmap_tpm <- function(genes){
  heatmap_plot <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
  heatmap_plot <- heatmap_plot %>%
    select(-contains("V6.")) %>% 
    select(-contains("Nafa.")) %>%
    select(contains(".CH."))
  heatmap.2(as.matrix(heatmap_plot), scale = "row", trace = 'none', density.info = 'none', col = viridis)
  
} #Plots heatmap version of boxplot_tpm
heatmap_average <- function(genes,mean_only){
  heatmap_plot <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
  heatmap_plot <- heatmap_plot %>%
    select(-contains("V6.")) %>% 
    select(-contains("Nafa.")) %>%
    select(contains(".CH."))
  names(heatmap_plot) <- sub(".M.",".",names(heatmap_plot))
  names(heatmap_plot) <- sub(".F.",".",names(heatmap_plot)) # Remove sex to perform rowmeans easily
  
  average_plot <- data.frame(row.names = row.names(heatmap_plot))
  average_plot <- heatmap_plot %>%
    filter(row.names(heatmap_plot) %in% genes) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.N")))) %>%
    mutate(Mock_HBE = rowMeans(select(.,contains("M.CH.B")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.N")))) %>%
    mutate(Virus_HBE = rowMeans(select(.,contains("V2.CH.B")))) %>%
    select(Mock_CH:Virus_HBE)
  
  Mean <- summarise_all(average_plot, mean)
  average_plot <- rbind(average_plot, Mean)
  rownames(average_plot)[rownames(average_plot) == "1"] <- "Mean"
  
  zscore <- data.frame(t(scale(t(average_plot))))
  col_fun = colorRamp2(c(-1.5,-0.75,0,0.75,1.5), 
                       c("#440154FF","#3B528BFF","#21908CFF","#5DC863FF","#FDE725FF"))
  
  if (mean_only == TRUE){
    Heatmap(as.matrix(zscore["Mean",]),
            cluster_rows = F, 
            cluster_columns = F, 
            col = col_fun, 
            top_annotation = HeatmapAnnotation("Tissue" = c("N","B","N","B"), 
                                               annotation_name_side = "left", 
                                               col = list("Tissue" = rep(c("N"="#8DD3C7","B"="#F7941D"),2)),
                                               gp = gpar(col = "black"),
                                               show_legend = F),
            heatmap_legend_param = list(
              title = "Gene Expression Z-score",
              direction = "horizontal",
              title_position = "topleft",
              at = c(-2,-1,0,1,2),
              legend_width = unit(4, "cm")))
  } else {
    Heatmap(as.matrix(zscore),
            cluster_rows = F, 
            cluster_columns = F, 
            col = col_fun, 
            top_annotation = HeatmapAnnotation("Tissue" = c("N","B","N","B"), 
                                               annotation_name_side = "left", 
                                               col = list("Tissue" = rep(c("N"="#8DD3C7","B"="#F7941D"),2)),
                                               gp = gpar(col = "black"),
                                               show_legend = F),
            heatmap_legend_param = list(
              title = "Gene Expression Z-score",
              direction = "horizontal",
              title_position = "topleft",
              at = c(-2,-1,0,1,2),
              legend_width = unit(4, "cm")))
  }
  
} #Heatmap of average for each condition + age group
wilcoxon_test <- function(genes){
  input_tpm <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
  zscore <- data.frame(t(scale(t(input_tpm))))
  zscore$gene <- row.names(zscore)
  zscore <- zscore %>%
    select(-contains("V6.")) %>% #Remove Virus 0.6
    select(-contains("Nafa.")) %>% # remove bronchial
    select(contains("CH.")) %>%
    mutate(gene = row.names(.))
  
  zscore_long <- reshape2::melt(zscore) #wide-form to long-form
  zscore_long <- zscore_long[order(zscore_long$gene),] #order by gene
  
  zscore_long$Condition <- ifelse(grepl("V2.", zscore_long$variable), "Virus", "Mock") 
  zscore_long$Age <- ifelse(grepl("CH", zscore_long$variable), "CH", 
                            ifelse(grepl("YA", zscore_long$variable), "YA","OA"))
  zscore_long$Tissue <- ifelse(grepl(".N.", zscore_long$variable),"HNE","HBE")
  zscore_long$Group <- paste0(zscore_long$Condition,"_",zscore_long$Tissue)
  
  .GlobalEnv$stat_test <- zscore_long %>% 
    wilcox_test(value ~ Group, p.adjust.method = "bonferroni")
  print(stat_test)
}

# Viral Boxplot Figure S2B ####

viral_tpm <- read.table(file = "viral_tpm.txt", sep = "\t", header = T, row.names = 1)

input_tpm <- viral_tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".YA")) %>%
  select(-contains(".OA")) %>%
  select(contains(".NT.") | contains(".TJ.") | contains(".IB.") |
           contains(".376AC.") | contains(".AM."))

input_tpm$Gene <- row.names(input_tpm)
input_tpm <- reshape2::melt(input_tpm)
input_tpm$Condition <- ifelse(grepl("V2.", input_tpm$variable), "Virus 0.2", 
                              ifelse(grepl("V6.", input_tpm$variable), "Virus 0.6","Nafa")) 
input_tpm$Age <- ifelse(grepl(".CH.", input_tpm$variable), "Child", 
                        ifelse(grepl("YA", input_tpm$variable), "Young Adult","Older Adult"))
input_tpm$Tissue <- ifelse(grepl("\\.N\\.", input_tpm$variable), "HNE", "HBE") 

input_tpm$Tissue <- factor(x = input_tpm$Tissue, levels = c("HNE","HBE"))

stat_test <- input_tpm %>%
  wilcox_test(value ~ Tissue, p.adjust.method = "bonferroni")
print(stat_test)

mean_virus <- input_tpm %>%
  group_by(Condition) %>%
  group_by(variable, add =T) %>%
  group_by(Age, add=T) %>%
  group_by(Tissue, add=T) %>%
  summarize(Mean=mean(value))

ggplot(mean_virus, aes(x=Tissue, y=Mean, fill=Tissue)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#8ED2C6","#F7941D")) +
  scale_color_manual(values = c("#8ED2C6","#F7941D")) +
  geom_point(aes(color=Tissue),position=position_dodge(width=0.1), size=3) +
  geom_text(aes(label=variable),position=position_dodge(width=0.3)) +
  theme_bw() +
  theme(legend.position = "none") # 3x5 portrait

#Average Virus for Figure 3E

tpm <- read.table(file = "viral_tpm_mockinc.txt", sep = "\t", header = T, row.names = 1)
tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B."))


virusgenes <- paste("GU280_gp", str_pad(1:11,2,pad=0),sep = "")

heatmap_average(genes = virusgenes, mean_only = T)

# DGE volcano plot - Figure S2C ####

dx <- read.table(file = "diffex_glm.CH_HBE.txt",sep="\t",header=T)

{
  dx$PValue <- as.numeric(as.character(dx$PValue))
  
  keyvals <- ifelse(
    dx$logFC <= -0.585 & dx$FDR < 0.05, '#E8495C',
    ifelse(dx$logFC >= 0.585 & dx$FDR < 0.05, '#38ACE2',
           'grey'))
  names(keyvals)[keyvals == '#38ACE2'] <- 'Upregulated'
  names(keyvals)[keyvals == 'grey'] <- 'Not significant'
  names(keyvals)[keyvals == '#E8495C'] <- 'Downregulated'
  
  EnhancedVolcano(dx,
                  lab = row.names(dx),
                  selectLab = c('SPIKE'),
                  x='logFC',
                  y='FDR',
                  pCutoff = 0.05,
                  FCcutoff = 0.585,
                  colCustom = keyvals,
                  colAlpha = 1, 
                  legendPosition = 'none',
                  xlim=c(-10,10), ylim=c(0,25)) #export pdf as 5.5x6
}

# IFNs Figure S2F ####

heatmap_average(genes= c("IFNB1","IFNE","IFNL1","IFNL2","IFNL3",
                          "IRF1","IRF2","IRF3","IRF5","IRF6","IRF7","IRF9",
                          "IFNAR1","IFNAR2","IFNGR1","IFNGR2","IFNLR1"),mean_only = F)
heatmap_average(genes=c("IFNAR1","IFNAR2","IFNGR1","IFNGR2","IFNLR1"), mean_only = T)

wilcoxon_test(genes = c("IFNB1","IFNE","IFNL1","IFNL2","IFNL3"))
wilcoxon_test(genes = c("IRF1","IRF2","IRF3","IRF5","IRF6","IRF7","IRF9"))
wilcoxon_test(genes = c("IFNAR1","IFNAR2","IFNGR1","IFNGR2","IFNLR1"))

# IFN Cluster plots Figure S2G ####

isgs <- read.table(file = "interferon_genes.txt", sep = "\t", header = T)

#Create heatmap of mean cluster expression across age groups
{
  C12 <- isgs %>%
    filter(Cluster == "1" | Cluster == "2") %>%
    pull(Gene)
  C3 <- isgs %>%
    filter(Cluster == "3") %>%
    pull(Gene)
  C4 <- isgs %>%
    filter(Cluster == "4") %>%
    pull(Gene)
  C5 <- isgs %>%
    filter(Cluster == "5") %>%
    pull(Gene)
  
  isg_tpm <- tpm
  names(isg_tpm) <- sub(".M.",".",names(isg_tpm))
  names(isg_tpm) <- sub(".F.",".",names(isg_tpm)) # Remove sex to perform rowmeans easily

  C12_mean <- isg_tpm %>%
    filter(row.names(isg_tpm) %in% C12) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.N")))) %>%
    mutate(Mock_HBE = rowMeans(select(.,contains("M.CH.B")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.N")))) %>%
    mutate(Virus_HBE = rowMeans(select(.,contains("V2.CH.B")))) %>%
    select(Mock_CH:Virus_HBE) %>%
    summarise_all(.,mean) %>%
    set_rownames("C1/2 Mean")
  
  C3_mean <- isg_tpm %>%
    filter(row.names(isg_tpm) %in% C3) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.N")))) %>%
    mutate(Mock_HBE = rowMeans(select(.,contains("M.CH.B")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.N")))) %>%
    mutate(Virus_HBE = rowMeans(select(.,contains("V2.CH.B")))) %>%
    select(Mock_CH:Virus_HBE) %>%
    summarise_all(.,mean) %>%
    set_rownames("C3 Mean")
  
  C4_mean <- isg_tpm %>%
    filter(row.names(isg_tpm) %in% C4) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.N")))) %>%
    mutate(Mock_HBE = rowMeans(select(.,contains("M.CH.B")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.N")))) %>%
    mutate(Virus_HBE = rowMeans(select(.,contains("V2.CH.B")))) %>%
    select(Mock_CH:Virus_HBE) %>%
    summarise_all(.,mean) %>%
    set_rownames("C4 Mean")
  
  C5_mean <- isg_tpm %>%
    filter(row.names(isg_tpm) %in% C5) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.N")))) %>%
    mutate(Mock_HBE = rowMeans(select(.,contains("M.CH.B")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.N")))) %>%
    mutate(Virus_HBE = rowMeans(select(.,contains("V2.CH.B")))) %>%
    select(Mock_CH:Virus_HBE) %>%
    summarise_all(.,mean) %>%
    set_rownames("C5 Mean")
  
  cluster_means <- bind_rows(C12_mean,C3_mean,C4_mean,C5_mean)
  
  zscore <- data.frame(t(scale(t(cluster_means))))
  col_fun = colorRamp2(c(-1.5,-0.75,0,0.75,1.5), 
                       c("#440154FF","#3B528BFF","#21908CFF","#5DC863FF","#FDE725FF"))
  
  Heatmap(as.matrix(zscore),
          cluster_rows = F, 
          cluster_columns = F, 
          col = col_fun, 
          top_annotation = HeatmapAnnotation("Tissue" = c("N","B","N","B"), 
                                             annotation_name_side = "left", 
                                             col = list("Tissue" = rep(c("N"="#8DD3C7","B"="#F7941D"),2)),
                                             gp = gpar(col = "black"),
                                             show_legend = F),
          heatmap_legend_param = list(
            title = "Gene Expression Z-score",
            direction = "horizontal",
            title_position = "topleft",
            legend_width = unit(4, "cm")))
}

heatmap_average(genes = C5,mean_only = F)

#Stats test 
#C1
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "1" | isgs$Cluster == "2"])
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "3"])
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "4"])
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "5"])

# Cytokine Plots - Figure S2H ####

heatmap_average(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                          "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                          "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                          "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"), mean_only = F)

heatmap_average(genes = c("IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                          "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"), mean_only = F)


wilcoxon_test(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                        "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                        "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                        "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"))

wilcoxon_test(genes = c("IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                        "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"))
# Viral Receptors and Entry Factors - Figure S2I ####

heatmap_average(genes = c("DDX58","DHX58","IFIH1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6"),mean_only = F)
wilcoxon_test(gen = c("DDX58","DHX58","IFIH1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6"))

heatmap_average(genes = c("ACE","ACE2","ANPEP","DPP4","ATP1A1","ATP1B1","NRP1",
                          "CTSL","FURIN","TMPRSS2","TMPRSS4"),mean_only = F)
wilcoxon_test(genes = c("ACE","ACE2","ANPEP","DPP4","ATP1A1","ATP1B1","NRP1",
                        "CTSL","FURIN","TMPRSS2","TMPRSS4"))


# NSP Interactome Z-Score - Figure S2J ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(contains("CH."))

interactome_genes <- read.table(file = "interactome_genes.txt", sep = "\t", header = T)

interactome_tpm <- data.frame((subset(tpm, rownames(tpm) %in% interactome_genes$Protein))) #Subset only genes in list
interactome_tpm <- interactome_tpm[complete.cases(interactome_tpm),] #Remove NAs

zscore <- data.frame(t(scale(t(interactome_tpm))))
zscore$geneID <- row.names(zscore)

zscore_long <- reshape2::melt(zscore) #wide-form to long-form
zscore_long <- zscore_long[order(zscore_long$geneID),] #order by gene

zscore_long$Condition <- ifelse(grepl("V2.", zscore_long$variable), "Virus", "Mock") 
zscore_long$Age <- ifelse(grepl("CH", zscore_long$variable), "CH", 
                          ifelse(grepl("YA", zscore_long$variable), "YA","OA"))
zscore_long$Age <- factor(x = zscore_long$Age, levels = c("OA","YA", "CH")) #Keeps plot in order of child > YA > OA
zscore_long$Tissue <- ifelse(grepl(".N.", zscore_long$variable),"HNE","HBE")
zscore_long$Tissue <- factor(x=zscore_long$Tissue, levels = c("HNE","HBE"))

zscore_long$Group <- paste0(zscore_long$Condition,"_",zscore_long$Tissue)


zscore_long <- merge(zscore_long, interactome_genes,  by.x = "geneID", by.y= "Protein")

stat_test <- zscore_long %>%
  group_by(Type) %>%
  wilcox_test(value ~ Group, p.adjust.method = "bonferroni")
stat_test

#All Types

ggplot(zscore_long, aes(x = Condition, y = value, fill = Tissue)) + 
  geom_violin() +
  theme_bw() +
  ylim(-4,4) +
  scale_fill_manual(values = c("#8DD3C7","#F7941D")) +
  theme(legend.position = "none") +
  facet_wrap(~Type, ncol = 5,
             scales = "fixed") +
  ylab("Z-score of gene expression") + 
  stat_summary(fun.y=median, geom="point", size=2, color="black", position=position_dodge(1))

#Figure 3G

stat_test <- subset(zscore_long, CoV2 %in% c("nsp7","nsp8","nsp12","nsp13")) %>%
  group_by(CoV2) %>%
  wilcox_test(value ~ Group, p.adjust.method = "bonferroni")
stat_test

ggplot(subset(zscore_long,CoV2 %in% c("nsp7","nsp8","nsp12","nsp13")), aes(x = Condition, y = value, fill = Tissue)) + 
  geom_violin() +
  theme_bw() +
  ylim(-4,4) +
  scale_fill_manual(values = c("#8DD3C7","#F7941D")) +
  theme(legend.position = "none") +
  facet_wrap(~CoV2, ncol = 5,
             scales = "fixed") +
  ylab("Z-score of gene expression") + 
  stat_summary(fun.y=mean, geom="point", size=2, color="black", position=position_dodge(1))

#Replication Complex

rep <- subset(zscore_long, Type %in% c("Replication complex"))

stat_test <- rep %>%
  group_by(CoV2) %>%
  group_by(Tissue, .add=T) %>%
  wilcox_test(value ~ Condition, p.adjust.method = "bonferroni")
stat_test

ggplot(rep, aes(x = Condition, y = value, fill = Age)) + 
  geom_violin() +
  theme_bw() +
  ylim(-4,4) +
  scale_fill_manual(values = rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2)) +
  theme(legend.position = "none") +
  facet_wrap(~CoV2, ncol = 5,
             scales = "fixed") +
  ylab("Z-score of gene expression") + 
  stat_summary(fun.y=median, geom="point", size=2, color="black", position=position_dodge(1))

#Capping Enzymes

cap <- subset(zscore_long, Type %in% c("Capping enzymes"))

stat_test <- cap %>%
  group_by(CoV2) %>%
  wilcox_test(value ~ Condition, p.adjust.method = "bonferroni")
stat_test

ggplot(cap, aes(x = Condition, y = value, fill = Age)) + 
  geom_violin() +
  theme_bw() +
  ylim(-4,4) +
  scale_fill_manual(values = rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2)) +
  theme(legend.position = "none") +
  facet_wrap(~CoV2, ncol = 5,
             scales = "fixed") +
  ylab("Z-score of gene expression") + 
  stat_summary(fun.y=median, geom="point", size=2, color="black", position=position_dodge(1))

#NSP1 

nsp1 <- subset(zscore_long, CoV2 %in% c("nsp1"))

stat_test <- nsp1 %>%
  wilcox_test(zscore ~ Condition, p.adjust.method = "bonferroni")
stat_test

ggplot(nsp1, aes(x = Condition, y = value, fill = Age)) + 
  geom_violin() +
  theme_bw() +
  ylim(-4,4) +
  scale_fill_manual(values = rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2)) +
  theme(legend.position = "none") +
  facet_wrap(~CoV2, ncol = 5,
             scales = "fixed") +
  ylab("Z-score of gene expression") + 
  stat_summary(fun.y=median, geom="point", size=2, color="black", position=position_dodge(1))
