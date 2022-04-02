#Code is used for plotting figures of proteomics analysis for COVID-19 2020 project
#Code developed by Alexander Capraro, September 2021
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

log2lfq <- read.table(file="log2_lfq.txt", sep = "\t", header=T)
sinfo <- read.table(file= "sampleinfo.txt", sep="\t", header=T)

# Protein Plots Functions - Virus 0.2 vs. Mock ####

heatmap_average <- function(genes,mean_only){
  heatmap_plot <- data.frame((subset(log2lfq, rownames(log2lfq) %in% genes))) #Subset only genes in list
  average_plot <- data.frame(row.names = row.names(heatmap_plot))
  
  average_plot <- log2lfq %>%
    filter(row.names(log2lfq) %in% genes) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("Mock.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("Mock.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("Mock.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    select(Mock_OA:Virus_CH)
  
  Mean <- summarise_all(average_plot, mean)
  average_plot <- rbind(average_plot, Mean)
  average_plot[average_plot == 0] <- NA 
  
  rownames(average_plot)[rownames(average_plot) == "1"] <- "Mean"
  
  zscore <- data.frame(t(scale(t(average_plot))))
  col_fun = colorRamp2(c(-1.5,-0.75,0,0.75,1.5), 
                       c("#440154FF","#3B528BFF","#21908CFF","#5DC863FF","#FDE725FF"))
  
  if (mean_only == TRUE){
    Heatmap(as.matrix(zscore["Mean",]),
            cluster_rows = F, 
            cluster_columns = F, 
            col = col_fun, 
            top_annotation = HeatmapAnnotation("Age Group" = rep(c("OA","YA","CH"),2), 
                                               annotation_name_side = "left", 
                                               col = list("Age Group" = rep(c("OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7"),2)),
                                               gp = gpar(col = "black"),
                                               show_legend = F),
            heatmap_legend_param = list(
              title = "Protein Expression Z-score",
              direction = "horizontal",
              title_position = "topleft",
              at = c(-2,-1,0,1,2),
              legend_width = unit(4, "cm")))
  } else {
    Heatmap(as.matrix(zscore),
            cluster_rows = F, 
            cluster_columns = F, 
            col = col_fun, 
            top_annotation = HeatmapAnnotation("Age Group" = rep(c("OA","YA","CH"),2), 
                                               annotation_name_side = "left", 
                                               col = list("Age Group" = rep(c("OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7"),2)),
                                               gp = gpar(col = "black"),
                                               show_legend = F),
            heatmap_legend_param = list(
              title = "Protein Expression Z-score",
              direction = "horizontal",
              title_position = "topleft",
              at = c(-2,-1,0,1,2),
              legend_width = unit(4, "cm")))
  }
  
} #Heatmap of average for each condition + age group
wilcoxon_test <- function(genes){
  input_tpm <- data.frame((subset(log2lfq, rownames(log2lfq) %in% genes))) #Subset only genes in list
  zscore <- data.frame(t(scale(t(input_tpm))))
  zscore$gene <- row.names(zscore)
  zscore <- zscore %>%
    select(-contains("V6.")) %>% #Remove Virus 0.6
    select(-contains("Nafa.")) %>%
    select(-contains(".B.")) # remove bronchial
  
  zscore_long <- reshape2::melt(zscore) #wide-form to long-form
  zscore_long <- zscore_long[order(zscore_long$gene),] #order by gene
  
  zscore_long$Condition <- ifelse(grepl("V2.", zscore_long$variable), "Virus", "Mock") 
  zscore_long$Age <- ifelse(grepl("CH", zscore_long$variable), "CH", 
                            ifelse(grepl("YA", zscore_long$variable), "YA","OA"))
  zscore_long$Group <- paste0(zscore_long$Condition,"_",zscore_long$Age)
  
  .GlobalEnv$stat_test <- zscore_long %>% 
    wilcox_test(value ~ Group, p.adjust.method = "bonferroni")
  print(stat_test)
}

# Chao et al Age Associated Protein Expression - Figure S1B ####

{log2lfq <- read.table(file="log2_lfq.txt", sep = "\t", header=T)

log2lfq <- log2lfq %>%
  select(contains("Mock."))

sinfo <- read.table(file= "sampleinfo.txt", sep="\t", header=T)
sinfo <- sinfo %>%
  filter(treatment == "Mock" & tissue == "N") %>%
  select(sname,individualid,agegroup,treatment,technicalrep) %>%
  mutate_at("sname",str_replace_all,"-",".")

lung_age <- read.table(file = "aging-genes-lung.revised.txt", sep = "\t", header = T)

lfq_up <- log2lfq %>%
  filter(row.names(.) %in% lung_age$genes)
zscore_up_long <- lfq_up %>%
  t() %>% scale() %>% t() %>% data.frame() %>%
  mutate("geneID" = row.names(lfq_up)) %>%
  reshape2::melt() %>%
  merge(., sinfo, by.x = "variable", by.y="sname") %>%
  merge(.,lung_age, by.x = "geneID", by.y="genes")

zscore_up_long$agegroup <- factor(zscore_up_long$agegroup, levels = c("OA","YA","CH"))
zscore_up_long$cluster <- factor(zscore_up_long$cluster, levels = c("AgeUp","AgeDown"))
my_comparisons <- list(c("CH", "OA"), c("CH", "YA"), c("YA", "OA"))

stat_test <- zscore_up_long %>%
  group_by(cluster) %>%
  wilcox_test(value ~ agegroup, p.adjust.method = "bonferroni")
stat_test

ggplot(zscore_up_long, aes(x = agegroup, y = value, fill = agegroup)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2))+
  ylim(-4,4) +
  facet_wrap(~cluster) #6.5x6
}

# Viral Protein Expression - Figure S1E ####

virus <- log2lfq %>%
  mutate(Gene.names = row.names(.)) %>%
  filter(grepl('_SARS2', Gene.names))

virus_input <- reshape2::melt(virus)
virus_input <- virus_input[order(virus_input$Gene.names),] #order by gene
virus_input$Age <- ifelse(grepl(".CH.", virus_input$variable), "Child", 
                          ifelse(grepl(".YA.", virus_input$variable), "Young Adult","Older Adult"))
virus_input$Tissue <- ifelse(grepl(".N.",virus_input$variable), "HNE", "HBE")
virus_input$Condition <- ifelse(grepl("Mock.", virus_input$variable), "Mock", 
                                ifelse(grepl("V2", virus_input$variable), "Virus 0.2",
                                       ifelse(grepl("V6", virus_input$variable), "Virus 0.6",
                                              ifelse(grepl("NafaM", virus_input$variable), "Nafa Mock",
                                                     ifelse(grepl("CuM", virus_input$variable), "Copper Mock",
                                                            ifelse(grepl("Nafa2", virus_input$variable), "Nafa 0.2",
                                                                   ifelse(grepl("Cu2", virus_input$variable), "Copper 0.2", "Copper 0.6")))))))
virus_input$Age <- factor(x = virus_input$Age, levels = c("Older Adult","Young Adult","Child"))
virus_input$Condition <- factor(x = virus_input$Condition, levels = c("Mock","Copper Mock","Nafa Mock","Virus 0.2","Nafa 0.2","Copper 0.2","Virus 0.6","Copper 0.6"))
virus_input <- sinfo %>%
  select(sname,individualid,agegroup,treatment,technicalrep,tissue,batch) %>%
  merge(virus_input,.,by.x="variable",by.y="sname")

plot_input <- virus_input %>%
  filter(Tissue == "HNE") %>% #filter out HBE
  filter(Condition == "Virus 0.2")

stat_test <- plot_input %>%
  wilcox_test(value ~ Age, p.adjust.method = "bonferroni")
stat_test

#Horizontal
ggplot(plot_input, aes (x = Age, y = value, fill = Age)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  ylab("log2 LFQ intensity") +
  ylim(15,25) +
  theme(legend.position = "none") +
  facet_wrap(~Gene.names)

# Nafa Virus Fold Change - Figure 2E ####

nafa_prot_fc <- read.table(file = "NafaVirus_FoldChange.txt", sep = "\t", header = T)

nafa_prot_fc %>%
  mutate(Age = factor(Age, levels = c("Older Adult","Young Adult","Child"))) %>%
  ggplot(aes(x=Age, y=FoldChange, fill=Age)) +
  geom_boxplot(alpha = 1) +
  geom_point(aes(size=0.3,fill=Age,color=Age)) +
  theme_bw() +
  ylim(-15,0) +
  scale_color_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  scale_fill_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  theme(legend.position = "none") #4x6

# Mean Viral Protein Expression HNE & HBE - Figure 1F ####

virus <- log2lfq %>%
  mutate(Gene.names = row.names(.)) %>%
  filter(grepl('_SARS2', Gene.names))

virus_input <- reshape2::melt(virus)
virus_input <- merge(virus_input, sinfo, by.x =2, by.y=11)
virus_input <- virus_input %>%
  filter(treatment == "Virus_0.2")
virus_input$agegroup <- paste(virus_input$agegroup,virus_input$tissue)
virus_input$agegroup <- factor(virus_input$agegroup, levels=c("OA N","YA N", "CH N","CH B"))

mean_virus_prot <- virus_input %>%
  filter(value > 0) %>%
  group_by(treatment) %>%
  group_by(variable, add =T) %>%
  group_by(tissue, add=T) %>%
  group_by(agegroup, add=T) %>%
  group_by(individualid, add=T) %>%
  summarize(Mean=mean(value)) %>%
  mutate(tissue = factor(tissue, levels=c("N","B")))

stat_test <- virus_input %>%
  wilcox_test(value ~ agegroup, p.adjust.method = "bonferroni")
stat_test

# Mean Age Comparison
ggplot(mean_virus_prot[which(mean_virus_prot$Mean > 0), ], aes (x = agegroup, y = Mean, fill = agegroup)) + 
  geom_boxplot(alpha = 1) +
  geom_point(aes(size=1,fill=agegroup,color=agegroup))+
  theme_bw() +
  scale_color_manual(values= c("#BEBADA","#FFFFB3","#8DD3C7","#F89521")) +
  scale_fill_manual(values= c("#BEBADA","#FFFFB3","#8DD3C7","#F89521")) +
  ylab("log2 LFQ intensity") +
  ylim(20,25) +
  theme(legend.position = "none")

# Nafa & 0.2 MOI - Viral Protein Expression - Figure S1F ####

virus <- log2lfq %>%
  mutate(Gene.names = row.names(.)) %>%
  filter(grepl('_SARS2', Gene.names))

virus_input <- reshape2::melt(virus)
virus_input <- virus_input[order(virus_input$Gene.names),] #order by gene
virus_input$Age <- ifelse(grepl(".CH.", virus_input$variable), "Child", 
                          ifelse(grepl(".YA.", virus_input$variable), "Young Adult","Older Adult"))
virus_input$Tissue <- ifelse(grepl(".N.",virus_input$variable), "HNE", "HBE")
virus_input$Condition <- ifelse(grepl("Mock.", virus_input$variable), "Mock", 
                                ifelse(grepl("V2", virus_input$variable), "Virus 0.2",
                                       ifelse(grepl("V6", virus_input$variable), "Virus 0.6",
                                              ifelse(grepl("NafaM", virus_input$variable), "Nafa Mock",
                                                     ifelse(grepl("CuM", virus_input$variable), "Copper Mock",
                                                            ifelse(grepl("Nafa2", virus_input$variable), "Nafa 0.2",
                                                                   ifelse(grepl("Cu2", virus_input$variable), "Copper 0.2", "Copper 0.6")))))))
virus_input$Age <- factor(x = virus_input$Age, levels = c("Older Adult","Young Adult","Child"))
virus_input$Condition <- factor(x = virus_input$Condition, levels = c("Mock","Copper Mock","Nafa Mock","Virus 0.2","Nafa 0.2","Copper 0.2","Virus 0.6","Copper 0.6"))

virus_input <- virus_input %>%
  filter(Tissue == "HNE") %>% #filter out HBE
  filter(Condition == "Virus 0.2")

stat_test <- virus_input %>%
  group_by(Condition) %>%
  wilcox_test(value ~ Age, p.adjust.method = "bonferroni")
stat_test

stat_test <- wilcox.test(virus_input$value ~ virus_input$Condition)
stat_test

#Age Comparison

mean_virus_prot <- virus_input %>%
  filter(value > 0) %>%
  group_by(Condition) %>%
  group_by(variable, add =T) %>%
  group_by(Age, add=T) %>%
  summarize(Mean=mean(value))

mean_virus_prot %>%
  filter(Mean > 0) %>%
  ggplot(aes (x = Age, y = Mean, fill = Age)) + 
  geom_boxplot(alpha = 1) +
  geom_point(aes(size=3,color=Age))+
  theme_bw() +
  scale_color_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  scale_fill_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  ylab("log2 LFQ intensity")  +
  theme(legend.position = "none") +
  ylim(15,25) +
  facet_wrap(~Condition, ncol = 4,
             scales = "fixed")




# Interferon, Cytokines, Viral receptors, Entry Factors - Figure S# ####

#IRFs 

heatmap_average(genes = c("IRF1","IRF2","IRF3","IRF5","IRF6","IRF7","IRF9"),mean_only = F)
wilcoxon_test(genes = c("IRF1","IRF2","IRF3","IRF5","IRF6","IRF7","IRF9"))

#Cytokines
heatmap_average(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                          "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                          "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                          "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"), mean_only = F)
wilcoxon_test(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                        "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                        "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                        "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"))
heatmap_average(genes = c("IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                          "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"), mean_only = F)
wilcoxon_test(genes = c("IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                        "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"))
#Receptors
heatmap_average(genes = c("DDX58","DHX58","IFIH1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6"),mean_only = F)
wilcoxon_test(genes = c("DDX58","DHX58","IFIH1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6"))
#Entry factors
heatmap_average(genes = c("ACE","ACE2","ANPEP","DPP4","ATP1A1","ATP1B1","NRP1",
                          "CTSL","FURIN","TMPRSS2","TMPRSS4"),mean_only = F)
wilcoxon_test(genes = c("ACE","ACE2","ANPEP","DPP4","ATP1A1","ATP1B1","NRP1",
                        "CTSL","FURIN","TMPRSS2","TMPRSS4"))
#Nsp8 & Nsp12 

interactome_genes <- read.table(file = "interactome_genes.txt", sep = "\t", header = T)
nsp8 <- interactome_genes[interactome_genes$CoV2 == "nsp8",]
nsp12 <- interactome_genes[interactome_genes$CoV2 == "nsp12",]

heatmap_average(genes = nsp8$Protein,mean_only = F)
heatmap_average(genes = nsp12$Protein,mean_only = F)
wilcoxon_test()
wilcoxon_test()

# ISG Clusters - Figure S2F ####

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
  C35 <- isgs %>%
    filter(Cluster == "3"| Cluster == "4" | Cluster == "5") %>%
    pull(Gene)
  
  C12_mean <- log2lfq %>%
    filter(row.names(log2lfq) %in% C12) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    select(Mock_OA:Virus_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C1/2 Mean")
  
  C3_mean <- log2lfq %>%
    filter(row.names(log2lfq) %in% C3) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    select(Mock_OA:Virus_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C3 Mean")
  
  C4_mean <- log2lfq %>%
    filter(row.names(log2lfq) %in% C4) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    select(Mock_OA:Virus_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C4 Mean")
  
  C5_mean <- log2lfq %>%
    filter(row.names(log2lfq) %in% C5) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    select(Mock_OA:Virus_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C5 Mean")
  
  C35_mean <- log2lfq %>%
    filter(row.names(log2lfq) %in% C35) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    select(Mock_OA:Virus_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C3-5 Mean")
  
  cluster_means <- bind_rows(C12_mean,C3_mean,C4_mean,C5_mean,C35_mean)
  
  zscore <- data.frame(t(scale(t(cluster_means))))
  
  Heatmap(as.matrix(zscore),
          cluster_rows = F, 
          cluster_columns = F, 
          col = viridis(10), 
          top_annotation = HeatmapAnnotation("Age Group" = rep(c("OA","YA","CH"),2), 
                                             annotation_name_side = "left", 
                                             col = list("Age Group" = rep(c("OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7"),2)),
                                             gp = gpar(col = "black"),
                                             show_legend = F),
          heatmap_legend_param = list(
            title = "Gene Expression Z-score",
            direction = "horizontal",
            title_position = "topleft",
            at = c(-2,0,2),
            legend_width = unit(4, "cm")))
}
#Stats test
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "1" | isgs$Cluster == "2"])
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "3"])
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "4"])
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "5"])
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "3" | isgs$Cluster == "4" | isgs$Cluster == "5"])

# Cell-type Markers - Figure 3C ####

cell_markers <- read.table(file = "cellmarkers.txt", sep = "\t", header = T)

prot <- log2lfq %>%
  mutate(Gene.names = row.names(.))

prot_input <- data.frame((subset(prot, rownames(prot) %in% cell_markers$Protein))) #Subset only genes in list
prot_input <- reshape2::melt(prot)
prot_input <- prot_input[order(prot_input$Gene.names),] #order by gene
prot_input$Age <- ifelse(grepl(".CH.", prot_input$variable), "Child", 
                         ifelse(grepl(".YA.", prot_input$variable), "Young Adult","Older Adult"))
prot_input$Tissue <- ifelse(grepl(".N.",prot_input$variable), "HNE", "HBE")
prot_input$Condition <- ifelse(grepl("Mock.", prot_input$variable), "Mock", 
                               ifelse(grepl("V2", prot_input$variable), "Virus 0.2",
                                      ifelse(grepl("V6", prot_input$variable), "Virus 0.6",
                                             ifelse(grepl("NafaM", prot_input$variable), "Nafa Mock",
                                                    ifelse(grepl("CuM", prot_input$variable), "Copper Mock",
                                                           ifelse(grepl("Nafa2", prot_input$variable), "Nafa 0.2",
                                                                  ifelse(grepl("Cu2", prot_input$variable), "Copper 0.2", "Copper 0.6")))))))
prot_input$Age <- factor(x = prot_input$Age, levels = c("Older Adult","Young Adult","Child"))
prot_input$Condition <- factor(x = prot_input$Condition, levels = c("Mock","Copper Mock","Nafa Mock","Virus 0.2","Nafa 0.2","Copper 0.2","Virus 0.6","Copper 0.6"))
prot_input <- merge(prot_input, cell_markers, by.x = "Gene.names", by.y = "Protein")

prot_input <- prot_input %>%
  filter(Tissue == "HNE") %>% #filter out HBE
  filter(Condition == "Mock" | Condition == "Virus 0.2")

ggplot(prot_input[prot_input$value > 0, ], aes (x = Celltype, y = value, fill = Age)) + 
  geom_violin(alpha = 1) +
  theme_bw() +
  scale_fill_manual(values=rep(c("#BEBADA","#FFFFB3","#8DD3C7"),3)) +
  theme(legend.position = "none") +
  ylab("log2 LFQ intensity")+
  stat_summary(fun=mean, geom="point", shape=16, size=2, position = position_dodge(0.9)) +
  facet_wrap(~Condition) #5.5x7

stat_test <- prot_input %>%
  group_by(Celltype) %>%
  group_by(Age, add=T) %>%
  wilcox_test(value ~ Condition, p.adjust.method = "bonferroni")
stat_test


# PCoA Analysis - Figure 3J ####

log2lfq <- read.table(file="log2_lfq.txt", sep = "\t", header=T)
sinfo <- read.table(file= "sampleinfo.txt", sep="\t", header=T)

group <- data.frame(names(log2lfq))

group <- sinfo %>%
  select(sname,individualid,agegroup,treatment,technicalrep,tissue) %>%
  merge(group,.,by.x="names.log2lfq.",by.y="sname")

mds <- plotMDS(log2lfq,plot = T,gene.selection = "common", top = 20000)

mds_plot <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  set_colnames(c("Dim1", "Dim2")) %>%
  rownames_to_column("SampleID")

gg <- merge(mds_plot,group,by.x="SampleID",by.y="names.log2lfq.")
gg$agegroup <- factor(gg$agegroup, levels = c("OA","YA","CH"))
gg <- gg %>%
  filter(treatment == "Mock" | treatment == "Virus_0.2") %>%
  filter(tissue == "N")

ggplot(gg, aes(x = Dim1, y = Dim2, colour = agegroup, fill = agegroup)) +
  geom_point(size = 3) + 
  scale_colour_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  scale_fill_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_text(aes(label=paste(individualid,treatment,tissue,sep=" ")),hjust=0, vjust=0)

#Vegan ADONIS

prot <- log2lfq %>%
  select(contains("Mock.") | contains("V2")) %>%
  select(contains(".N.")) %>%
  select(gg$SampleID)

adon.results<-adonis(t(prot) ~ gg$treatment * gg$agegroup,perm=999)
print(adon.results)

# RNA-Protein Correlation - Figure S3D ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
log2lfq <- read.table(file="log2_lfq.txt", sep = "\t", header=T)

{rna <- tpm %>%
  select(contains("M.CH") | contains("M.YA") | contains("M.OA"))
prot <- log2lfq %>%
  select(contains("Mock."))

rna$meanRNA <- rowMeans(rna)
prot$meanProt <- rowMeans(prot)

RNA_prot <- merge(rna,prot,by = 0) %>%
  set_rownames(.$Row.names) %>%
  select(meanRNA,meanProt) %>%
  filter(meanProt > 0) #remove all missing data

corr <- cor.test(x=RNA_prot$meanRNA, y=RNA_prot$meanProt, method = 'spearman')
print(corr)
ggplot(RNA_prot, aes(x = meanProt, y = meanRNA)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method=lm, color='red') + 
  ylim(0,15)
} # Mock

{rna <- tpm %>%
    select(contains("V2."))
  prot <- log2lfq %>%
    select(contains("V2."))
  
  rna$meanRNA <- rowMeans(rna)
  prot$meanProt <- rowMeans(prot)
  
  RNA_prot <- merge(rna,prot,by = 0) %>%
    set_rownames(.$Row.names) %>%
    select(meanRNA,meanProt) %>%
    filter(meanProt > 0) #remove all missing data
  
  corr <- cor.test(x=RNA_prot$meanRNA, y=RNA_prot$meanProt, method = 'spearman')
  print(corr)
  ggplot(RNA_prot, aes(x = meanProt, y = meanRNA)) +
    geom_point() +
    theme_bw() +
    geom_smooth(method=lm, color='red') + 
    ylim(0,15)
} # Virus 0.2

# Copper - Viral Protein Expression - Figure 4E ####

virus <- log2lfq %>%
  mutate(Gene.names = row.names(.)) %>%
  filter(grepl('_SARS2', Gene.names))

virus_input <- reshape2::melt(virus)
virus_input <- virus_input[order(virus_input$Gene.names),] #order by gene
virus_input$Age <- ifelse(grepl(".CH.", virus_input$variable), "Child", 
                          ifelse(grepl(".YA.", virus_input$variable), "Young Adult","Older Adult"))
virus_input$Tissue <- ifelse(grepl(".N.",virus_input$variable), "HNE", "HBE")
virus_input$Condition <- ifelse(grepl("Mock.", virus_input$variable), "Mock", 
                                ifelse(grepl("V2", virus_input$variable), "Virus 0.2",
                                       ifelse(grepl("V6", virus_input$variable), "Virus 0.6",
                                              ifelse(grepl("NafaM", virus_input$variable), "Nafa Mock",
                                                     ifelse(grepl("CuM", virus_input$variable), "Copper Mock",
                                                            ifelse(grepl("Nafa2", virus_input$variable), "Nafa 0.2",
                                                                   ifelse(grepl("Cu2", virus_input$variable), "Copper 0.2", "Copper 0.6")))))))
virus_input$Age <- factor(x = virus_input$Age, levels = c("Older Adult","Young Adult","Child"))
virus_input$Condition <- factor(x = virus_input$Condition, levels = c("Mock","Copper Mock","Nafa Mock","Virus 0.2","Nafa 0.2","Copper 0.2","Virus 0.6","Copper 0.6"))

virus_input <- virus_input %>%
  filter(Tissue == "HNE") %>% #filter out HBE
  filter(!grepl("Mock", Condition)) %>% # filter out mocks
  filter(!grepl("Nafa", Condition)) %>%
  filter(grepl("466BM", variable) | grepl("475HL",variable) | grepl("489ED",variable) | grepl("412JA",variable)
         | grepl("380HK",variable) | grepl("476LT",variable) | grepl("024AT",variable) | grepl("022RG",variable)
         | grepl("023RB",variable) | grepl("025PT",variable) | grepl("010SLW",variable) | grepl("019MB",variable)
         | grepl("021CsH",variable))

mean_virus_prot <- virus_input %>%
  group_by(Condition) %>%
  group_by(variable, add =T) %>%
  group_by(Age, add=T) %>%
  summarize(Mean=mean(value))

stat_test <- virus_input %>%
  wilcox_test(value ~ Condition, p.adjust.method = "bonferroni")
stat_test

#Age Comparison
ggplot(mean_virus_prot[which(mean_virus_prot$Mean > 0), ], aes (x = Condition, y = Mean, fill = Condition)) + 
  geom_boxplot(alpha = 1) +
  geom_point(aes(size=1,fill=Condition,color=Condition))+
  theme_bw() +
  scale_color_manual(values=c("#A6CEE2","#1F78B4","#B4D88A","#34A048")) +
  scale_fill_manual(values=c("#A6CEE2","#1F78B4","#B4D88A","#34A048")) +
  ylab("log2 LFQ intensity")  +
  facet_wrap(~Age, ncol = 4,
             scales = "fixed") +
  theme(legend.position = "none")

#All ages
ggplot(mean_virus_prot[which(mean_virus_prot$Mean > 0), ], aes (x = Condition, y = Mean, fill = Condition)) + 
  geom_boxplot(alpha = 1) +
  geom_point(aes(size=1,fill=Condition,color=Condition))+
  theme_bw() +
  scale_color_manual(values=c("#A6CEE2","#1F78B4","#B4D88A","#34A048")) +
  scale_fill_manual(values=c("#A6CEE2","#1F78B4","#B4D88A","#34A048")) +
  ylab("log2 LFQ intensity") +
  theme(legend.position = "none")

virus_input %>% filter(Condition == "Virus 0.2" | Condition == "Virus 0.6") %>%
  filter(value > 0) %>%
  ggplot(aes (x = Condition, y = value, fill = Condition)) + 
  geom_boxplot(alpha = 1) +
  theme_bw() +
  ylab("log2 LFQ intensity") +
  theme(legend.position = "none")


# Copper PCoA - Figure # ####

log2lfq <- read.table(file="log2_lfq.txt", sep = "\t", header=T)
sinfo <- read.table(file= "sampleinfo.txt", sep="\t", header=T)

group <- data.frame(names(log2lfq))

group <- sinfo %>%
  select(sname,individualid,agegroup,treatment,technicalrep,tissue) %>%
  merge(group,.,by.x="names.log2lfq.",by.y="sname")

mds <- plotMDS(log2lfq,plot = T,gene.selection = "common", top = 20000)

mds_plot <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  set_colnames(c("Dim1", "Dim2")) %>%
  rownames_to_column("SampleID")

gg <- merge(mds_plot,group,by.x="SampleID",by.y="names.log2lfq.")
gg$agegroup <- factor(gg$agegroup, levels = c("OA","YA","CH"))
gg <- gg %>%
  filter(treatment == "Virus_0.2" | treatment == "CuCl_0.2" | treatment == "Mock") %>%
  filter(tissue == "N")%>%
  filter(grepl("466BM", individualid) | grepl("475HL",individualid) | grepl("489ED",individualid) | grepl("412JA",individualid)
         | grepl("380HK",individualid) | grepl("476LT",individualid) | grepl("024AT",individualid) | grepl("022RG",individualid)
         | grepl("023RB",individualid) | grepl("025PT",individualid) | grepl("010SLW",individualid) | grepl("019MB",individualid)
         | grepl("021CsH",individualid))

ggplot(gg, aes(x = Dim1, y = Dim2, colour = agegroup, fill = agegroup)) +
  geom_point(size = 5,aes(shape=treatment)) + 
  scale_colour_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  scale_fill_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  theme_bw() +
  geom_text(aes(label=paste(individualid,treatment,tissue,sep=" ")),hjust=0, vjust=0)

#Vegan ADONIS

prot <- log2lfq %>% select(contains("V6.") | contains("Cu6.")) %>%
  select(contains(".466BM.") | contains(".475HL.") | contains(".489ED.") | contains(".412JA.") | contains(".380HK.") | contains(".476LT.") |
           contains("010SLW") | contains("019MB") | contains("021CsH") |
           contains("024AT") | contains("022RG") | contains("023RB") | contains("025PT"))

adon.results<-adonis(t(prot) ~ gg$treatment * gg$agegroup,perm=999)
print(adon.results)

# Protein Plots Functions - Virus 0.2 vs. Mock ####

log2lfq <- read.table(file="log2_lfq.txt", sep = "\t", header=T)
log2lfq <- log2lfq %>% select(contains(".N.")) %>%
  select(contains("466BM") | contains("475HL") | contains("489ED") | contains("412JA") | contains("380HK") | contains("476LT") |
           contains("010SLW") | contains("019MB") | contains("021CsH") |
           contains("024AT") | contains("022RG") | contains("023RB") | contains("025PT"))

heatmap_average2 <- function(genes,mean_only){
  heatmap_plot <- data.frame((subset(log2lfq, rownames(log2lfq) %in% genes))) #Subset only genes in list
  average_plot <- data.frame(row.names = row.names(heatmap_plot))
  
  average_plot <- log2lfq %>%
    filter(row.names(log2lfq) %in% genes) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("Mock.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("Mock.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("Mock.CH.")))) %>%
    
    mutate(CuMock_OA = rowMeans(select(.,contains("CuM.OA.")))) %>%
    mutate(CuMock_YA = rowMeans(select(.,contains("CuM.YA.")))) %>%
    mutate(CuMock_CH = rowMeans(select(.,contains("CuM.CH.")))) %>%
    
    mutate(NafaMock_OA = rowMeans(select(.,contains("NafaM.OA.")))) %>%
    mutate(NafaMock_YA = rowMeans(select(.,contains("NafaM.YA.")))) %>%
    mutate(NafaMock_CH = rowMeans(select(.,contains("NafaM.CH.")))) %>%
    
    mutate(Virus0.2_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus0.2_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus0.2_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    
    mutate(Cu0.2_OA = rowMeans(select(.,contains("Cu2.OA.")))) %>%
    mutate(Cu0.2_YA = rowMeans(select(.,contains("Cu2.YA.")))) %>%
    mutate(Cu0.2_CH = rowMeans(select(.,contains("Cu2.CH.")))) %>%
    
    mutate(Nafa0.2_OA = rowMeans(select(.,contains("Nafa2.OA.")))) %>%
    mutate(Nafa0.2_YA = rowMeans(select(.,contains("Nafa2.YA.")))) %>%
    mutate(Nafa0.2_CH = rowMeans(select(.,contains("Nafa2.CH.")))) %>%
    
    mutate(Virus0.6_OA = rowMeans(select(.,contains("V6.OA.")))) %>%
    mutate(Virus0.6_YA = rowMeans(select(.,contains("V6.YA.")))) %>%
    mutate(Virus0.6_CH = rowMeans(select(.,contains("V6.CH.")))) %>%
    
    mutate(Cu0.6_OA = rowMeans(select(.,contains("Cu6.OA.")))) %>%  
    mutate(Cu0.6_YA = rowMeans(select(.,contains("Cu6.YA.")))) %>%  
    mutate(Cu0.6_CH = rowMeans(select(.,contains("Cu6.CH.")))) %>%
    
    select(Mock_OA:Cu0.6_CH)
  
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
            top_annotation = HeatmapAnnotation("Age Group" = c("OA","YA","CH", "OA","YA","CH", "OA","YA","CH", "OA","YA","CH",
                                                               "OA","YA","CH", "OA","YA","CH", "OA","YA","CH", "OA","YA","CH"), 
                                               annotation_name_side = "left", 
                                               col = list("Age Group" = c("OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7", 
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7")),
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
            top_annotation = HeatmapAnnotation("Age Group" = c("OA","YA","CH", "OA","YA","CH", "OA","YA","CH", "OA","YA","CH",
                                                               "OA","YA","CH", "OA","YA","CH", "OA","YA","CH", "OA","YA","CH"), 
                                               annotation_name_side = "left", 
                                               col = list("Age Group" = c("OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7", 
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7",
                                                                          "OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7")),
                                               gp = gpar(col = "black"),
                                               show_legend = F),
            heatmap_legend_param = list(
              title = "Gene Expression Z-score",
              direction = "horizontal",
              title_position = "topleft",
              at = c(-2,-1,0,1,2),
              legend_width = unit(4, "cm")))
  }
}
wilcoxon_test2 <- function(genes){
  
  input_lfq <- data.frame((subset(log2lfq, rownames(log2lfq) %in% genes))) #Subset only genes in list
  zscore <- data.frame(t(scale(t(input_lfq))))
  zscore$gene <- row.names(zscore)
  zscore <- zscore %>%
    select(-contains(".B.")) # remove bronchial
  
  zscore_long <- reshape2::melt(zscore) #wide-form to long-form
  zscore_long <- zscore_long[order(zscore_long$gene),] #order by gene
  
  zscore_long$Condition <- ifelse(grepl("V2.", zscore_long$variable), "Virus 0.2", 
                                  ifelse(grepl("V6.", zscore_long$variable), "Virus 0.6", 
                                         ifelse(grepl("Nafa.", zscore_long$variable), "Nafa", "Mock")))
  zscore_long$Age <- ifelse(grepl("CH", zscore_long$variable), "CH", 
                            ifelse(grepl("YA", zscore_long$variable), "YA","OA"))
  zscore_long$Group <- paste0(zscore_long$Condition,"_",zscore_long$Age)
  
  .GlobalEnv$stat_test <- zscore_long %>% 
    wilcox_test(value ~ Group, p.adjust.method = "bonferroni")
  print(stat_test)
}

#Virus

heatmap_average2(genes = c("SPIKE_SARS2","AP3A_SARS2","VME1_SARS2","NS7A_SARS2","NS8_SARS2","NCAP_SARS2","R1AB_SARS2","ORF9B_SARS2"), mean_only=T)

#IRFs 

heatmap_average2(genes = c("IRF1","IRF2","IRF3","IRF5","IRF6","IRF7","IRF9"),mean_only = T)
wilcoxon_test(genes = c("IRF1","IRF2","IRF3","IRF5","IRF6","IRF7","IRF9"))

#Cytokines
heatmap_average2(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                          "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                          "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                          "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"), mean_only = T)
wilcoxon_test2(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                        "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                        "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                        "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"))
heatmap_average2(genes = c("IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                          "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"), mean_only = T)
wilcoxon_test(genes = c("IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                        "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"))
#Receptors
heatmap_average2(genes = c("DDX58","DHX58","IFIH1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6"),mean_only = T)
wilcoxon_test(genes = c("DDX58","DHX58","IFIH1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6"))
#Entry factors
heatmap_average2(genes = c("ACE","ACE2","ANPEP","DPP4","ATP1A1","ATP1B1","NRP1",
                          "CTSL","FURIN","TMPRSS2","TMPRSS4"),mean_only = T)
wilcoxon_test(genes = c("ACE","ACE2","ANPEP","DPP4","ATP1A1","ATP1B1","NRP1",
                        "CTSL","FURIN","TMPRSS2","TMPRSS4"))
#Nsp8 & Nsp12 

interactome_genes <- read.table(file = "interactome_genes.txt", sep = "\t", header = T)
nsp8.12 <- interactome_genes[interactome_genes$CoV2 == "nsp8" | interactome_genes$CoV2 == "nsp12",]

heatmap_average2(genes = nsp8.12$Protein,mean_only = T)
heatmap_average2(genes = nsp12$Protein,mean_only = T)
wilcoxon_test()
wilcoxon_test()
# Nafa & 0.6 MOI - Viral Protein Expression - Figure S4C ####

virus <- log2lfq %>%
  mutate(Gene.names = row.names(.)) %>%
  filter(grepl('_SARS2', Gene.names))

{virus_input <- reshape2::melt(virus)
virus_input <- virus_input[order(virus_input$Gene.names),] #order by gene
virus_input$Age <- ifelse(grepl(".CH.", virus_input$variable), "Child", 
                          ifelse(grepl(".YA.", virus_input$variable), "Young Adult","Older Adult"))
virus_input$Tissue <- ifelse(grepl(".N.",virus_input$variable), "HNE", "HBE")
virus_input$Condition <- ifelse(grepl("Mock.", virus_input$variable), "Mock", 
                                ifelse(grepl("V2", virus_input$variable), "Virus 0.2",
                                       ifelse(grepl("V6", virus_input$variable), "Virus 0.6",
                                              ifelse(grepl("NafaM", virus_input$variable), "Nafa Mock",
                                                     ifelse(grepl("CuM", virus_input$variable), "Copper Mock",
                                                            ifelse(grepl("Nafa2", virus_input$variable), "Nafa 0.2",
                                                                   ifelse(grepl("Cu2", virus_input$variable), "Copper 0.2", "Copper 0.6")))))))
virus_input$Age <- factor(x = virus_input$Age, levels = c("Older Adult","Young Adult","Child"))
virus_input$Condition <- factor(x = virus_input$Condition, levels = c("Mock","Copper Mock","Nafa Mock","Virus 0.2","Nafa 0.2","Copper 0.2","Virus 0.6","Copper 0.6"))
  }

virus_input <- virus_input %>%
  filter(Tissue == "HNE") %>% #filter out HBE
  filter(!grepl("Mock", Condition)) %>% # filter out mocks
  filter(Condition != "Virus 0.6") %>%
  filter(Condition != "Copper 0.6")
mean_virus_prot <- virus_input %>%
  filter(value > 0) %>%
  group_by(Condition) %>%
  group_by(variable, add =T) %>%
  group_by(Age, add=T) %>%
  summarize(Mean=mean(value))

mean_prot_all <- mean_virus_prot %>%
  mutate(Age = "All") %>%
  rbind(., mean_virus_prot) %>%
  mutate(Age = factor(Age, levels = c("All","Older Adult","Young Adult","Child"))) %>%
  mutate(Condition = factor(Condition, levels=c("Virus 0.2","Copper 0.2", "Nafa 0.2")))

test <- virus_input %>%
  filter(Condition == "Virus 0.2" | Condition == "Nafa 0.2")
wilcox.test(value ~ Condition,data = test)

test <- virus_input %>%
  filter(Condition == "Copper 0.2" | Condition == "Nafa 0.2")
wilcox.test(value ~ Condition,data = test)

test <- virus_input %>%
  filter(Condition == "Virus 0.2" | Condition == "Copper 0.2")
wilcox.test(value ~ Condition,data = test)

# Mean Age Comparison
ggplot(mean_virus_prot[which(mean_virus_prot$Mean > 0), ], aes (x = Age, y = Mean, fill = Age)) + 
  geom_boxplot(alpha = 1) +
  geom_point(aes(size=1,fill=Age,color=Age))+
  theme_bw() +
  scale_color_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  scale_fill_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  ylab("log2 LFQ intensity")  +
  facet_wrap(~Condition, ncol = 4,
             scales = "fixed") +
  theme(legend.position = "none") +
  ylim(0,30)

# Mean All ages
ggplot(mean_prot_all[mean_prot_all$Mean > 0, ], aes(x=Condition, y=Mean, fill=Condition)) +
  geom_boxplot(alpha = 1) +
  geom_point(aes(size=0.3,fill=Condition,color=Condition))+
  theme_bw() +
  scale_color_manual(values=c("#A7A9AC","#A47F5D","#90C4C4")) +
  scale_fill_manual(values=c("#A7A9AC","#A47F5D","#90C4C4")) +
  theme(legend.position = "none") +
  facet_wrap(~Age, ncol = 4)

stat_test <- virus_input %>%
  wilcox_test(value ~ Condition, p.adjust.method = "bonferroni")
stat_test

# Nafa & 0.6 MOI - PCoA - Figure S4J ####

log2lfq <- read.table(file="log2_lfq.txt", sep = "\t", header=T)
sinfo <- read.table(file= "sampleinfo.txt", sep="\t", header=T)

group <- data.frame(names(log2lfq))

group <- sinfo %>%
  select(sname,individualid,agegroup,treatment,technicalrep,tissue) %>%
  merge(group,.,by.x="names.log2lfq.",by.y="sname")

mds <- plotMDS(log2lfq,plot = T,gene.selection = "common", top = 20000)

mds_plot <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  set_colnames(c("Dim1", "Dim2")) %>%
  rownames_to_column("SampleID")

gg <- merge(mds_plot,group,by.x="SampleID",by.y="names.log2lfq.")
gg$agegroup <- factor(gg$agegroup, levels = c("OA","YA","CH"))
gg <- gg %>%
  filter(tissue == "N") %>% #filter out HBE
  filter(!grepl("Mock", treatment)) %>% # filter out mocks
  filter(!grepl("Cu", treatment)) %>%
  filter(grepl("466BM", individualid) | grepl("475HL",individualid) | grepl("489ED",individualid) | grepl("412JA",individualid)
         | grepl("380HK",individualid) | grepl("476LT",individualid) | grepl("024AT",individualid) | grepl("022RG",individualid)
         | grepl("023RB",individualid) | grepl("025PT",individualid) | grepl("010SLW",individualid) | grepl("019MB",individualid)
         | grepl("021CsH",individualid))

ggplot(gg, aes(x = Dim1, y = Dim2, colour = agegroup, fill = agegroup)) +
  geom_point(size = 5,aes(shape=treatment)) + 
  scale_colour_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  scale_fill_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_text(aes(label=paste(individualid,treatment,tissue,sep=" ")),hjust=0, vjust=0)

#Vegan ADONIS

prot <- log2lfq %>% select(!contains("Mock.") | !contains("Cu")) %>%
  select(contains(".466BM.") | contains(".475HL.") | contains(".489ED.") | contains(".412JA.") | contains(".380HK.") | contains(".476LT.") |
           contains("010SLW") | contains("019MB") | contains("021CsH") |
           contains("024AT") | contains("022RG") | contains("023RB") | contains("025PT"))

adon.results<-adonis(t(prot) ~ gg$treatment * gg$agegroup,perm=999)
print(adon.results)

