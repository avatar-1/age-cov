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
library(reshape2)
library(usedist)

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)

tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B."))

# Gene Plots Functions ####

boxplot_tpm <- function(genes, scales, cols){
  
  input_tpm <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
  input_tpm$gene <- row.names(input_tpm)
  input_tpm <- input_tpm %>%
    select(-contains(".B.")) # remove bronchial
  
  input_tpm <- reshape2::melt(input_tpm) #wide-form to long-form
  input_tpm <- input_tpm[order(input_tpm$gene),] #order by gene
  
  input_tpm$Condition <- ifelse(grepl("V2.", input_tpm$variable), "Virus 0.2", 
                                ifelse(grepl("Nafa.", input_tpm$variable), "Nafa",
                                       ifelse(grepl("V6.", input_tpm$variable), "Virus 0.6","Mock")) )
  
  input_tpm$Age <- ifelse(grepl(".CH", input_tpm$variable), "Child", 
                          ifelse(grepl("YA", input_tpm$variable), "Young Adult","Older Adult"))
  input_tpm$Age <- factor(x = input_tpm$Age, levels = c("Older Adult","Young Adult", "Child")) #Keeps plot in order of child > YA > OA
  input_tpm$gene <- factor(x = input_tpm$gene, levels = genes) #Keeps order of genes
  
  #Plot individual samples
  ggplot(input_tpm, aes (x = Condition, y = value, fill = Age)) + 
    geom_boxplot(alpha = 1) +
    geom_point(aes(fill = Age), size = 2, shape = 21, position = position_jitterdodge()) +
    theme_bw() +
    scale_fill_manual(values=rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2))+ 
    theme(legend.position = "none") +
    facet_wrap(~gene, ncol = cols,
               scales = scales) +
    ylab("log2(TPM)")
}
heatmap_tpm <- function(genes){
  heatmap_plot <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
  heatmap_plot <- heatmap_plot %>%
    select(-contains("V6.")) %>% #Remove Virus 0.6
    select(-contains(".B.")) # remove bronchial
  heatmap.2(as.matrix(heatmap_plot), scale = "row", trace = 'none', density.info = 'none', col = viridis)
  
} #Plots heatmap version of boxplot_tpm
heatmap_average <- function(genes,mean_only){
  heatmap_plot <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
  average_plot <- data.frame(row.names = row.names(heatmap_plot))
  
  average_plot <- tpm %>%
    filter(row.names(tpm) %in% genes) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    select(Mock_OA:Virus_CH)
  
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
            top_annotation = HeatmapAnnotation("Age Group" = rep(c("OA","YA","CH"),2), 
                                               annotation_name_side = "left", 
                                               col = list("Age Group" = rep(c("OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7"),2)),
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
            top_annotation = HeatmapAnnotation("Age Group" = rep(c("OA","YA","CH"),2), 
                                               annotation_name_side = "left", 
                                               col = list("Age Group" = rep(c("OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7"),2)),
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
boxplot_zscore <- function(genes, scales, cols){
  input_tpm <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
  zscore <- data.frame(t(scale(t(input_tpm))))
  zscore$gene <- row.names(zscore)
  zscore <- zscore %>%
    select(-contains("V6.")) %>% #Remove Virus 0.6
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
  
  zscore_long$Age <- factor(x = zscore_long$Age, levels = c("OA","YA", "CH")) #Keeps plot in order of child > YA > OA
  
  ggplot(zscore_long, aes(x = Condition, y = value, fill = Age)) + 
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = rep(c("OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7"),2)) +
    ylim(-3,4) +
    facet_wrap(~gene,
               scales = scales,
               ncol = cols) +
    ylab("Z-score of Gene Expression")
  #export 4.5x6
}
wilcoxon_test <- function(genes){
  input_tpm <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
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

# Viral Heatmap - Figure 1E ####

viral_tpm <- read.table(file = "viral_tpm.txt", sep = "\t", header = T, row.names = 1)

input_tpm <- viral_tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B."))

input_tpm$Gene <- row.names(input_tpm)
input_tpm <- reshape2::melt(input_tpm)
input_tpm$Condition <- ifelse(grepl("V2.", input_tpm$variable), "Virus 0.2", 
                              ifelse(grepl("V6.", input_tpm$variable), "Virus 0.6","Nafa")) 
input_tpm$Age <- ifelse(grepl(".CH.", input_tpm$variable), "Child", 
                        ifelse(grepl("YA", input_tpm$variable), "Young Adult","Older Adult"))
input_tpm$Tissue <- ifelse(grepl("\\.N\\.", input_tpm$variable), "HNE", "HBE") 
input_tpm$Age <- factor(x = input_tpm$Age, levels = c("Older Adult", "Young Adult","Child")) #Keeps plot in order of child > YA > OA


stat_test <- input_tpm %>%
  wilcox_test(value ~ Age, p.adjust.method = "bonferroni")
print(stat_test)

ggplot(input_tpm, aes (x = Age, y = value, fill = Age)) + 
  geom_boxplot(alpha = 1) +
  theme_minimal() +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7")) +
  theme(legend.position = "none") +
  ylab("log2(TPM)")

zscore <- data.frame(t(scale(t(viral_tpm)))) %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B.")) %>%
  select(-contains(".AC.t2")) %>% #Remove technical replicates
  select(-contains(".NT.t2")) %>%
  select(-contains(".NiT.t2")) %>%
  select(-contains(".WB.t2")) %>%
  select(-contains(".JW.t2"))

agegroup <- ifelse(grepl(".CH", colnames(zscore)), "Child", 
                   ifelse(grepl("YA", colnames(zscore)), "Young Adult","Older Adult"))


Heatmap(as.matrix(zscore),
        cluster_rows = F, 
        cluster_columns = T, 
        col = viridis(10), 
        top_annotation = HeatmapAnnotation(type = agegroup, annotation_name_side = "left",
                                           gp = gpar(col = "black"),
                                           show_legend = T),
        heatmap_legend_param = list(
          title = "Gene Expression Z-score",
          direction = "horizontal",
          title_position = "topleft",
          at = c(-2,0,2),
          legend_width = unit(4, "cm")))

#Average Virus for Figure 3E

tpm <- read.table(file = "viral_tpm_mockinc.txt", sep = "\t", header = T, row.names = 1)
tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B."))


virusgenes <- paste("GU280_gp", str_pad(1:11,2,pad=0),sep = "")

heatmap_average(genes = virusgenes, mean_only = T)
wilcoxon_test(genes = virusgenes)

#Batch effect

sinfo <- read.table(file = "sinfo.txt", sep = "\t", header = T)

viral_tpm$Gene <- row.names(viral_tpm)
batch_tpm <- reshape2::melt(viral_tpm)

batch_tpm <- sinfo %>%
  select(sname,individualid,agegroup,treatment,technicalrepid,tissue,experiment) %>%
  mutate_at("sname",str_replace_all,"-",".") %>%
  merge(batch_tpm,.,by.x="variable",by.y="sname")

batch_tpm <- batch_tpm %>%
  filter(tissue == "N") %>% #filter out HBE
  filter(treatment == "Virus" | treatment == "Virus6") %>%
  mutate(experiment = factor(experiment)) %>%
  mutate(agegroup = factor(agegroup, levels = c("OA","YA","CH")))

ggplot(batch_tpm, aes(x=treatment, y=value)) +
  geom_boxplot(aes(fill=batch)) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~agegroup)

stat_test <- batch_tpm %>%
  group_by(agegroup) %>%
  wilcox_test(value ~ treatment, p.adjust.method = "bonferroni")
print(stat_test)

# Virus Boxplot V2 + Nafa - Figure 1F ####

viral_tpm <- read.table(file = "viral_tpm.txt", sep = "\t", header = T, row.names = 1)

input_tpm <- viral_tpm %>%
  select(-contains("V6")) %>%
  select(-contains(".B."))

input_tpm$Gene <- row.names(input_tpm)
input_tpm <- reshape2::melt(input_tpm)
input_tpm$Condition <- ifelse(grepl("V2.", input_tpm$variable), "Virus 0.2", 
                              ifelse(grepl("V6.", input_tpm$variable), "Virus 0.6","Nafa")) 
input_tpm$Age <- ifelse(grepl(".CH.", input_tpm$variable), "Child", 
                        ifelse(grepl("YA", input_tpm$variable), "Young Adult","Older Adult"))
input_tpm$Tissue <- ifelse(grepl("\\.N\\.", input_tpm$variable), "HNE", "HBE") 
input_tpm$AgeGroup <- paste(input_tpm$Age,input_tpm$Tissue)
input_tpm$AgeGroup <- factor(x = input_tpm$AgeGroup, levels = c("Older Adult HNE", "Young Adult HNE","Child HNE","Child HBE")) #Keeps plot in order of child > YA > OA

stat_test <- input_tpm %>%
  group_by(Condition) %>%
  wilcox_test(value ~ AgeGroup, p.adjust.method = "bonferroni")
print(stat_test)

stat_test <- input_tpm %>%
  wilcox_test(value ~ Condition, p.adjust.method = "bonferroni")
print(stat_test)

mean_virus <- input_tpm %>%
  group_by(Condition, add=T) %>%
  group_by(variable, add =T) %>%
  group_by(AgeGroup, add=T) %>%
  summarize(Mean=mean(value)) %>%
  mutate(AgeGroup = factor(AgeGroup, levels = c("Older Adult HNE","Young Adult HNE","Child HNE", "Child HBE"))) %>%
  mutate(Condition = factor(Condition, levels = c("Virus 0.2","Nafa")))

ggplot(mean_virus, aes (x = AgeGroup, y = Mean, fill = AgeGroup)) + 
  geom_boxplot(alpha = 1) +
  geom_point(aes(size=3,color=AgeGroup)) +
  theme_bw() +
  scale_color_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F89521")) +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F89521")) +
  theme(legend.position = "none") +
  ylab("log2(TPM)") + 
  facet_wrap(~Condition)+
  ylim(0,18)#pdf 6x6


# Chao et al - Age Associated Genes - Figure S1B ####

{tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)

tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains("V2.")) %>%
  select(-contains(".B."))

sinfo <- read.table(file = "sinfo.txt", sep = "\t", header = T)
sinfo <- sinfo %>%
  filter(treatment == "Mock" & tissue == "N") %>%
  select(sname,individualid,agegroup,treatment,technicalrepid) %>%
  mutate_at("sname",str_replace_all,"-",".")

lung_age <- read.table(file = "Chow et al 2020/aging-genes-lung.revised.txt", sep = "\t", header = T)

tpm_up <- tpm %>%
  filter(row.names(.) %in% lung_age$genes)
zscore_up_long <- tpm_up %>%
  t() %>% scale() %>% t() %>% data.frame() %>%
  mutate("geneID" = row.names(tpm_up)) %>%
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

# Technical Replicate Viral Reads Figure ####
tpm <- read.table(file = "viral_tpm.txt", sep = "\t", header = T)
sinfo <- read.table(file = "sinfo.txt", sep = "\t", header = T)
tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B."))
tpm$gene <- row.names(tpm)
  
{input_tpm <- reshape2::melt(tpm) #wide-form to long-form
  input_tpm <- input_tpm[order(input_tpm$gene),] #order by gene
  
  input_tpm$Condition <- ifelse(grepl("V2.", input_tpm$variable), "Virus", "Mock") 
  
  input_tpm$Age <- ifelse(grepl(".CH", input_tpm$variable), "Child", 
                          ifelse(grepl("YA", input_tpm$variable), "Young Adult","Older Adult"))
  input_tpm$Age <- factor(x = input_tpm$Age, levels = c("Older Adult","Young Adult", "Child")) #Keeps plot in order of child > YA > OA
  
  input_tpm$RepID <- ifelse(grepl(".t1", input_tpm$variable), "1", "2") 
  input_tpm <- sinfo %>%
    select(sname,individualid) %>%
    mutate_at("sname",str_replace_all,"-",".") %>%
    merge(input_tpm,.,by.x="variable",by.y="sname") #Get age from sinfo and merge with mean_ratio_virus
  
  rep_tpm <- filter(input_tpm, individualid == "AC" | 
                      individualid == "JW" |
                      individualid == "Nit" |
                      individualid == "NT" |
                      individualid == "WB")
  rep_tpm$individualid <- factor(x=rep_tpm$individualid, levels = c("NT","AC","JW","WB"))
  
  
  #Plot individual samples
  plot <- ggplot(rep_tpm, aes (x = Age, y = value, fill=RepID)) + 
    geom_boxplot(alpha = 1) +
    geom_point(aes(fill = RepID), size = 2, shape = 21, position = position_dodge(0.75)) +
    theme_bw() +
    ylab(label = "Viral Expression (log2(TPM))") +
    xlab(label = "Individual ID") +
    guides(fill=guide_legend(title="Technical Replicate"))
  plot(plot)
  
  stat_test <- rep_tpm %>%
    group_by(Age) %>%
    wilcox_test(value ~ RepID, p.adjust.method = "bonferroni")
  print(stat_test)
}
  

# Interferon Plots - Figure 2A####

#IFNs, IRFs & IFNRs 

heatmap_average(genes = c("IFNB1","IFNE","IFNL1","IFNL2","IFNL3"),mean_only = F)
heatmap_average(genes = c("IRF1","IRF2","IRF3","IRF5","IRF6","IRF7","IRF9"),mean_only = F)
heatmap_average(genes = c("IFNAR1","IFNAR2","IFNGR1","IFNGR2","IFNLR1"),mean_only = F)

wilcoxon_test(genes = c("IFNB1","IFNE","IFNL1","IFNL2","IFNL3"))
wilcoxon_test(genes = c("IRF1","IRF2","IRF3","IRF5","IRF6","IRF7","IRF9"))
wilcoxon_test(genes = c("IFNAR1","IFNAR2","IFNGR1","IFNGR2","IFNLR1"))

# IFN Correlation - Figure 2B ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains(".B."))

#Use for Mock vs. Virus 0.2

sinfo <- read.table("sinfo.txt", sep="\t",header=T)
genes = c("IFNB1","IFNL1","IFNL2","IFNL3")

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

viral_tpm <- read.table(file = "viral_tpm_mockinc.txt", sep = "\t", header = T, row.names = 1)
viral_tpm <- viral_tpm %>%
  select(-contains("V6")) %>%
  select(-contains(".B."))

viral_tpm <- viral_tpm %>%
  colSums %>%
  as.data.frame()

viral_tpm$Condition <- ifelse(grepl("V2.", row.names(viral_tpm)), "Virus", 
                              ifelse(grepl("Nafa.", row.names(viral_tpm)), "Nafa", "Mock"))
viral_tpm$Age <- ifelse(grepl("CH", row.names(viral_tpm)), "Child", 
                        ifelse(grepl("YA", row.names(viral_tpm)), "Young Adult","Older Adult"))
viral_tpm$Age <- factor(x = viral_tpm$Age, levels = c("Older Adult","Young Adult", "Child")) #Keeps plot in order of child > YA > OA

mean_ratio <- merge(input_tpm, viral_tpm, by.x = "variable", by.y = 0)
mean_ratio$Group <- paste0(mean_ratio$Condition.x,"_",mean_ratio$Age.x)

mean_ratio_dt <- data.table(mean_ratio)
corr_test <- mean_ratio_dt[ ,cor.test(x=value, y=., method = 'spearman')[-2], by=Age.x]
corr_test

mean_ratio %>%
  ggplot(aes(x = ., y = value, color=Age.x)) +
  geom_point(aes(size=1, shape=gene)) +
  scale_color_manual(values =  rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_linetype_discrete(name="Age.x") +
  geom_smooth(method=lm, formula = y ~ x,fullrange=T) +
  xlab("SARS-CoV-2 Gene Expression (log2 TPM)")+
  ylab("Interferon Gene Expression (log2 TPM)") #+
  stat_poly_eq(formula = y ~ x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse=TRUE,label.x.npc = "right") #use to get R^2


#Virus age correlation Figure 

mean_ratio_virus <- mean_ratio %>%
  filter(Condition.x != "Mock")

mean_ratio_virus <- sinfo %>%
  select(sname, age) %>%
  mutate_at("sname",str_replace_all,"-",".") %>%
  merge(mean_ratio_virus,.,by.x="variable",by.y="sname") #Get age from sinfo and merge with mean_ratio_virus

corr_test <- cor.test(x=mean_ratio_virus$value, y=mean_ratio_virus$age, method = 'spearman')
corr_test

ggplot(mean_ratio_virus, aes(x = age, y = value)) +
  geom_point(aes(size=1,shape=gene, color = Group)) +
  scale_color_manual(values = c("#8DD3C7","#BEBADA","#FFFFB3")) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_smooth(method=lm, color="red", fill = "red", formula = y ~ x)+
  xlab("Age (years)")+
  ylab("Interferon Gene Expression (log2 TPM)") +
  xlim(0,70) +
  stat_poly_eq(formula = y ~ x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse=TRUE,label.x.npc = "right")


# IFN Cluster Plots - Figure 2C ####

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
  
  C12_mean <- tpm %>%
    filter(row.names(tpm) %in% C12) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
    mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
    mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
    select(Mock_OA:Nafa_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C1/2 Mean")
  
  C3_mean <- tpm %>%
    filter(row.names(tpm) %in% C3) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
    mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
    mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
    select(Mock_OA:Nafa_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C3 Mean")
  
  C4_mean <- tpm %>%
    filter(row.names(tpm) %in% C4) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
    mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
    mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
    select(Mock_OA:Nafa_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C4 Mean")
  
  C5_mean <- tpm %>%
    filter(row.names(tpm) %in% C5) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
    mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
    mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
    select(Mock_OA:Nafa_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C5 Mean")
  
  C35_mean <- tpm %>%
    filter(row.names(tpm) %in% C35) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
    mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
    mutate(Virus_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Virus_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Virus_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
    mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
    mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
    select(Mock_OA:Nafa_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C3-5 Mean")
  
  cluster_means <- bind_rows(C12_mean,C3_mean,C4_mean,C5_mean,C35_mean)
  
  zscore <- data.frame(t(scale(t(cluster_means))))
  
  Heatmap(as.matrix(zscore),
          cluster_rows = F, 
          cluster_columns = F, 
          col = viridis(10), 
          top_annotation = HeatmapAnnotation("Age Group" = rep(c("OA","YA","CH"),3), 
                                             annotation_name_side = "left", 
                                             col = list("Age Group" = rep(c("OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7"),3)),
                                             gp = gpar(col = "black"),
                                             show_legend = F),
          heatmap_legend_param = list(
            title = "Gene Expression Z-score",
            direction = "horizontal",
            title_position = "topleft",
            at = c(-2,0,2),
            legend_width = unit(4, "cm")))
}

heatmap_average(genes = C5,mean_only = F)

#Stats test 
#C1
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "1" | isgs$Cluster == "2"])
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "3"])
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "4"])
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "5"])
wilcoxon_test(genes = isgs$Gene[isgs$Cluster == "3" | isgs$Cluster == "4" | isgs$Cluster == "5"])




# Chemokines and Interleukins - Figure 2D ####

heatmap_average(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                          "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                          "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                          "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"), mean_only = T)

heatmap_average(genes = c("IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                          "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"), mean_only = T)


wilcoxon_test(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                         "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                         "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                         "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"))
wilcoxon_test(genes = c("IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                        "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"))





# CIBERSORTx - Figure 4B ####

celltype_prop <- read.table(file = "CIBERSORTx_Output_Normalised_NoTransitional.txt", sep = "\t", header = T)

long_prop <- celltype_prop %>%
  reshape2::melt() %>%
  mutate(Condition = ifelse(grepl("V2.", .$Mixture), "Virus 0.2", 
                            ifelse(grepl("V6.", .$Mixture), "Virus 0.6",
                                   ifelse(grepl("Nafa.", .$Mixture), "Nafa","Mock")))) %>%
  mutate(Condition = factor(.$Condition, levels= c("Mock","Virus 0.2", "Nafa","Virus 0.6"))) %>%
  mutate(Age = ifelse(grepl(".CH.", .$Mixture), "Child", 
                      ifelse(grepl("YA", .$Mixture), "Young Adult","Older Adult"))) %>%
  mutate(Age = factor(.$Age, levels= c("Older Adult","Young Adult", "Child"))) %>%
  mutate(Tissue = ifelse(grepl(".N.", .$Mixture), "HNE", "HBE")) %>%
  mutate(value = value*100) %>% #converts decimal into percentage
  filter(Tissue == "HNE")

ggplot(long_prop, aes(x = Age, y = value, fill = variable)) + 
  geom_bar(position = "stack", stat = "summary", alpha = 0.9, fun = "mean") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  facet_wrap(~Condition, ncol = 4)#6x7.5

ggplot(long_prop, aes(x = Condition, y = value, fill = variable)) + 
  geom_bar(position = "stack", stat = "summary", alpha = 0.9, fun = "mean") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()

stat_test <- long_prop %>%
  filter(Condition == "Mock") %>%
  group_by(variable) %>%
  wilcox_test(value ~ Age, p.adjust.method = "bonferroni", paired = F)
stat_test

stat_test <- long_prop %>%
  filter(Condition == "Virus 0.2") %>%
  group_by(variable) %>%
  wilcox_test(value ~ Age, p.adjust.method = "bonferroni", paired = F)
stat_test

stat_test <- long_prop %>%
  group_by(variable, add = TRUE) %>%
  wilcox_test(value ~ Condition, p.adjust.method = "bonferroni", paired = F)
stat_test

#Ciliated heatmap

t_celltypes <- celltype_prop %>%
  set_rownames(.$Mixture) %>% 
  select(Ciliated,Basal,Secretory) %>%
  t() %>% data.frame

celltype_heatmap <- t_celltypes %>%
  mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
  mutate(Mock_YA = rowMeans(select(.,contains("M.YA.")))) %>%
  mutate(Mock_CH = rowMeans(select(.,contains("M.CH.")))) %>%
  mutate(Virus0.2_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
  mutate(Virus0.2_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
  mutate(Virus0.2_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
  mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
  mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
  mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
  mutate(Virus0.6_OA = rowMeans(select(.,contains("V6.OA.")))) %>%
  mutate(Virus0.6_YA = rowMeans(select(.,contains("V6.YA.")))) %>%
  mutate(Virus0.6_CH = rowMeans(select(.,contains("V6.CH.")))) %>%
  select(Mock_OA:Virus0.6_CH)

celltype_zscore <- data.frame(t(scale(t(celltype_heatmap))))
col_fun = colorRamp2(c(-1.5,-0.75,0,0.75,1.5), 
                     c("#440154FF","#3B528BFF","#21908CFF","#5DC863FF","#FDE725FF"))

Heatmap(as.matrix(celltype_zscore),
        cluster_rows = F, 
        cluster_columns = F, 
        col = col_fun, 
        top_annotation = HeatmapAnnotation("Age Group" = rep(c("OA","YA","CH"),4), 
                                           annotation_name_side = "left", 
                                           col = list("Age Group" = rep(c("OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7"),4)),
                                           gp = gpar(col = "black"),
                                           show_legend = F),
        heatmap_legend_param = list(
          title = "Gene Expression Z-score",
          direction = "horizontal",
          title_position = "topleft",
          at = c(-2,-1,0,1,2),
          legend_width = unit(4, "cm")))

# Viral Receptors and Entry Factors - Figure 4C####

heatmap_average(genes = c("DDX58","DHX58","IFIH1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6"),mean_only = F)
wilcoxon_test(gen = c("DDX58","DHX58","IFIH1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6"))

heatmap_average(genes = c("ACE","ACE2","ANPEP","DPP4","ATP1A1","ATP1B1","NRP1",
                          "CTSL","FURIN","TMPRSS2","TMPRSS4"),mean_only = F)
wilcoxon_test(genes = c("ACE","ACE2","ANPEP","DPP4","ATP1A1","ATP1B1","NRP1",
                        "CTSL","FURIN","TMPRSS2","TMPRSS4"))

# Entry Factor Expression - SARS-CoV-2 Expression Correlation - Figure 4D ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B."))

genes = c("ACE2","TMPRSS2")

{input_tpm <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
input_tpm$gene <- row.names(input_tpm)

input_tpm$gene <- row.names(input_tpm)
input_tpm <- reshape2::melt(input_tpm)
input_tpm <- input_tpm[order(input_tpm$gene),] #order by gene

input_tpm$Condition <- ifelse(grepl("V2.", input_tpm$variable), "Virus", "Mock") 
input_tpm$Age <- ifelse(grepl(".CH", input_tpm$variable), "Child", 
                        ifelse(grepl("YA", input_tpm$variable), "Young Adult","Older Adult"))
input_tpm$Age <- factor(x = input_tpm$Age, levels = c("Older Adult","Young Adult", "Child")) #Keeps plot in order of child > YA > OA
input_tpm$gene <- factor(x = input_tpm$gene, levels = genes) #Keeps order of genes

viral_tpm <- read.table(file = "viral_tpm.txt", sep = "\t", header = T, row.names = 1)
viral_tpm <- viral_tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B."))

viral_tpm <- viral_tpm %>%
  colSums %>%
  as.data.frame()

viral_tpm$Condition <- ifelse(grepl("V2.", row.names(viral_tpm)), "Virus", "Mock") 
viral_tpm$Age <- ifelse(grepl("CH", row.names(viral_tpm)), "Child", 
                        ifelse(grepl("YA", row.names(viral_tpm)), "Young Adult","Older Adult"))
viral_tpm$Age <- factor(x = viral_tpm$Age, levels = c("Older Adult","Young Adult", "Child")) #Keeps plot in order of child > YA > OA

mean_ratio <- merge(input_tpm, viral_tpm, by.x = "variable", by.y = 0)
mean_ratio$Group <- paste0(mean_ratio$Condition.x,"_",mean_ratio$Age.x)

mean_ratio_dt <- data.table(mean_ratio)
}

corr_test <- mean_ratio_dt[ ,cor.test(x=value, y=., method = 'spearman')[-2], by=gene]
corr_test

ggplot(mean_ratio, aes(x = ., y = value,fill=gene)) +
  geom_point(aes(shape=gene,color=Age.x),size = 3) +
  scale_fill_manual(values =  c("#03AFBB","#E6B821")) +
  scale_color_manual(values =  c("#BEBADA","#FFFFB3","#8DD3C7")) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_smooth(method=lm, formula = y ~ x,fullrange=T) +
  xlab("SARS-CoV-2 Gene Expression (log2 TPM)")+
  ylab("Entry Factor Expression (log2 TPM)") +
  stat_poly_eq(formula = y ~ x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
             parse=TRUE,label.x.npc = "right") #use to get R^2


# NSP Interactome Z-Score - Figure 4E ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B."))

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
zscore_long$Group <- paste0(zscore_long$Condition,"_",zscore_long$Age)
zscore_long <- merge(zscore_long, interactome_genes,  by.x = "geneID", by.y= "Protein")

stat_test <- zscore_long %>%
  group_by(Type) %>%
  wilcox_test(value ~ Group, p.adjust.method = "bonferroni")
stat_test

#All Types

ggplot(zscore_long, aes(x = Condition, y = value, fill = Age)) + 
  geom_violin() +
  theme_bw() +
  ylim(-4,4) +
  scale_fill_manual(values = rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2)) +
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

ggplot(subset(zscore_long,CoV2 %in% c("nsp7","nsp8","nsp12","nsp13")), aes(x = Condition, y = value, fill = Age)) + 
  geom_violin() +
  theme_bw() +
  ylim(-4,4) +
  scale_fill_manual(values = rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2)) +
  theme(legend.position = "none") +
  facet_wrap(~CoV2, ncol = 5,
             scales = "fixed") +
  ylab("Z-score of gene expression") + 
  stat_summary(fun.y=mean, geom="point", size=2, color="black", position=position_dodge(1))

zscore_long %>%
  filter(CoV2 == c("nsp7","nsp8","nsp12","nsp13")) %>%
  mutate(CoV2 = factor(CoV2, levels = c("nsp7","nsp8","nsp12","nsp13"))) %>%
  ggplot(aes(x = Condition, y = value, fill = Age)) + 
  geom_boxplot() +
  theme_bw() +
  ylim(-3,3) +
  scale_fill_manual(values = rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2)) +
  theme(legend.position = "none") +
  facet_wrap(~CoV2, ncol = 5,
             scales = "fixed") +
  ylab("Z-score of gene expression") + 
  stat_summary(fun.y=mean, geom="point", size=1.5, color="black", position=position_dodge(0.75))


#Replication Complex

rep <- subset(zscore_long, Type %in% c("Replication complex"))

stat_test <- rep %>%
  group_by(Condition) %>%
  wilcox_test(value ~ Age, p.adjust.method = "bonferroni")
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
# Nsp Virus Correlation - Figure 4H ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm %>%
  select(-contains("V6")) %>%
  select(-contains("Nafa")) %>%
  select(-contains(".B."))

interactome_genes <- read.table(file = "interactome_genes.txt", sep = "\t", header = T)

#Nsp8
genes <- interactome_genes %>%
  subset(CoV2 %in% c("nsp8"))

#Nsp12
genes <- interactome_genes %>%
  subset(CoV2 %in% c("nsp12"))

genes <- interactome_genes %>%
  subset(CoV2 %in% c("nsp8","nsp12","nsp13")) %>%
  select(Protein)

heatmap_average(genes = nsps$Protein,mean_only = T)
wilcoxon_test(genes = nsps$Protein)


{
  input_tpm <- data.frame((subset(tpm, rownames(tpm) %in% genes$Protein))) #Subset only genes in list
  input_tpm$gene <- row.names(input_tpm)
  
  input_tpm$gene <- row.names(input_tpm)
  input_tpm <- reshape2::melt(input_tpm)
  input_tpm <- input_tpm[order(input_tpm$gene),] #order by gene
  
  input_tpm$Condition <- ifelse(grepl("V2.", input_tpm$variable), "Virus", "Mock") 
  input_tpm$Age <- ifelse(grepl(".CH", input_tpm$variable), "Child", 
                          ifelse(grepl("YA", input_tpm$variable), "Young Adult","Older Adult"))
  input_tpm$Age <- factor(x = input_tpm$Age, levels = c("Older Adult","Young Adult", "Child")) #Keeps plot in order of child > YA > OA
  input_tpm$gene <- factor(x = input_tpm$gene, levels = genes$Protein) #Keeps order of genes
  
  viral_tpm <- read.table(file = "viral_tpm.txt", sep = "\t", header = T, row.names = 1)
  viral_tpm <- viral_tpm %>%
    select(-contains("V6")) %>%
    select(-contains("Nafa")) %>%
    select(-contains(".B."))
  
  viral_tpm <- viral_tpm %>%
    colSums %>%
    as.data.frame()
  
  viral_tpm$Condition <- ifelse(grepl("V2.", row.names(viral_tpm)), "Virus", "Mock") 
  viral_tpm$Age <- ifelse(grepl("CH", row.names(viral_tpm)), "Child", 
                          ifelse(grepl("YA", row.names(viral_tpm)), "Young Adult","Older Adult"))
  viral_tpm$Age <- factor(x = viral_tpm$Age, levels = c("Older Adult","Young Adult", "Child")) #Keeps plot in order of child > YA > OA
  
  mean_ratio <- merge(input_tpm, viral_tpm, by.x = "variable", by.y = 0)
  mean_ratio$Group <- paste0(mean_ratio$Condition.x,"_",mean_ratio$Age.x)
  
  corr_test <- cor.test(x=mean_ratio$value, y=mean_ratio$., method = 'spearman')
  
  
  plot <- ggplot(mean_ratio, aes(x = ., y = value)) +
    geom_point(color = "#00AFBB") +
    theme_bw() +
    geom_smooth(method=lm, color="#00AFBB", fill = "#00AFBB")+
    xlab("SARS-CoV-2 Gene Expression (log2 TPM)")+
    ylab("Interacting Protein Expression (log2 TPM)") +
    stat_poly_eq(formula = y ~ x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse=TRUE,label.x.npc = "right") #use to get R^2
  
  plot(plot)
  print(corr_test)
  }

# PCoA Analysis - Figure 4I ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm %>%
  select(-contains("Nafa")) %>%
  select(-contains("V6")) %>%
  select(contains(".N."))

sinfo <- read.table(file = "sinfo.txt", sep = "\t", header = T)
viral_percent <- read.table(file = "viral_percent_incMock.txt", sep = "\t", header = T)

group <- data.frame(names(tpm))

group <- sinfo %>%
  select(sname,individualid,agegroup,treatment,technicalrepid) %>%
  mutate_at("sname",str_replace_all,"-",".") %>%
  merge(group,.,by.x="names.tpm.",by.y="sname") %>%
  mutate(group,Condition = paste0(agegroup,"_",treatment))

group <- merge(group,viral_percent,by.x="names.tpm.",by.y="sname")

group <- group[ order(group$names.tpm.), ]
mds <- plotMDS(tpm,plot = T,gene.selection = "common", top = 20000)

mds_plot <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  set_colnames(c("Dim1", "Dim2")) %>%
  rownames_to_column("SampleID")

gg <- merge(mds_plot,group,by.x="SampleID",by.y="names.tpm.")

ggplot(gg, aes(x = Dim1, y = Dim2, colour = Condition, fill = Condition)) +
  geom_point(size = 3) + 
  scale_colour_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_text(aes(label=paste(individualid,Condition)),hjust=0, vjust=0)

#Vegan ADONIS

counts <- tpm %>%
  select(gg$SampleID)

adon.results<-adonis(t(tpm) ~ gg$agegroup + gg$treatment + gg$percentage,perm=999)
print(adon.results)

#betadispr
groups <- sub( "(^[^.]+[.][^.]+)(.+$)", "\\1", names(tpm)) #Get Condition + Age Group
dst <- vegdist(t(tpm))
disp <- betadisper(dst, groups)
disp
plot(disp, hull = FALSE, ellipse = F, label = F)

anova(disp)
permutest(disp, pairwise = T)

df <- data.frame(dist=disp$distances,group=disp$group)

summ <- df %>% 
  group_by(group) %>% 
  summarise(mean = mean(dist), median = median(dist), sd = sd(dist))
summ

ggplot(df,aes(x=group,y=dist,fill=group)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_bw()

# Distances between centroids

CH_Virus <- grep("V2.CH", names(tpm), value = T)
YA_Virus <- grep("V2.YA", names(tpm), value = T)
OA_Virus <- grep("V2.OA", names(tpm), value = T)

dist_CH <- list()
dist_YA <- list()
dist_OA <- list()

for (i in CH_Virus){
  dist <- dist_between_centroids(dst, grep("M.CH",names(tpm)), i)
  dist_CH[[i]] <- dist
}
for (i in YA_Virus){
  dist <- dist_between_centroids(dst, grep("M.YA",names(tpm)), i)
  dist_YA[[i]] <- dist
}
for (i in OA_Virus){
  dist <- dist_between_centroids(dst, grep("M.OA",names(tpm)), i)
  dist_OA[[i]] <- dist
}

dist_CH <- data.frame(matrix(unlist(dist_CH), nrow=length(dist_CH), byrow=TRUE)) %>%
  set_rownames(CH_Virus) %>%
  rename("Dist_to_Mock" = 1) %>%
  mutate(AgeGroup = "CH")

dist_YA <- data.frame(matrix(unlist(dist_YA), nrow=length(dist_YA), byrow=TRUE)) %>%
  set_rownames(YA_Virus) %>%
  rename("Dist_to_Mock" = 1) %>%
  mutate(AgeGroup = "YA")

dist_OA <- data.frame(matrix(unlist(dist_OA), nrow=length(dist_OA), byrow=TRUE)) %>%
  set_rownames(OA_Virus) %>%
  rename("Dist_to_Mock" = 1) %>%
  mutate(AgeGroup = "OA")

dist_to_mock_centroid <- rbind(dist_CH,dist_YA,dist_OA)

dist_to_mock_centroid %>%
  mutate(AgeGroup = factor(AgeGroup, levels=c("OA","YA","CH"))) %>%
  ggplot(aes(x=AgeGroup, y=Dist_to_Mock, fill=AgeGroup)) +
  geom_boxplot() +
  geom_point(aes(size=3,color=AgeGroup)) +
  stat_summary(fun=mean, geom="point", shape=20, size=3) +
  scale_color_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  scale_fill_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  theme_bw()

dist_to_mock_centroid %>%
  wilcox_test(Dist_to_Mock ~ AgeGroup)

tapply(dist_to_mock_centroid$Dist_to_Mock,
       dist_to_mock_centroid$AgeGroup,
       summary)




dist_between_centroids(dst, grep("M.CH",names(tpm)), grep("V2.CH", names(tpm))) #Child


dist_between_centroids(dst, grep("M.YA",names(tpm)), grep("V2.YA",names(tpm))) #Young adult
dist_between_centroids(dst, grep("M.OA",names(tpm)), grep("V2.OA",names(tpm))) #Older adult

# Nafa & 0.6 Virus Heatmap - Figure S3B ####

viral_tpm <- read.table(file = "viral_tpm.txt", sep = "\t", header = T, row.names = 1)
viral_tpm <- viral_tpm %>% select(-contains(".B.")) %>%
  select(contains(".BM.") | contains(".HK.") | contains(".LT.") | contains("HL") |
           contains("021CsH") | contains(".MB.") | contains(".SLW.") |
           contains("RB") | contains("PT") | contains("RG") | contains("AT"))
input_tpm <- viral_tpm

input_tpm$Gene <- row.names(input_tpm)
input_tpm <- reshape2::melt(input_tpm)
input_tpm$Condition <- ifelse(grepl("V2.", input_tpm$variable), "Virus 0.2", 
                              ifelse(grepl("V6.", input_tpm$variable), "Virus 0.6","Nafa"))
input_tpm$Condition <- factor(input_tpm$Condition, levels = c("Virus 0.2", "Nafa", "Virus 0.6"))
input_tpm$Age <- ifelse(grepl("CH", input_tpm$variable), "Child", 
                        ifelse(grepl("YA", input_tpm$variable), "Young Adult","Older Adult"))
input_tpm$Tissue <- ifelse(grepl("\\.N\\.", input_tpm$variable), "HNE", "HBE") 

input_tpm$Age <- factor(x = input_tpm$Age, levels = c("Older Adult", "Young Adult","Child")) #Keeps plot in order of child > YA > OA

stat_test <- input_tpm %>%
  group_by(Condition) %>%
  wilcox_test(value ~ Condition, p.adjust.method = "bonferroni")
print(stat_test)

ggplot(input_tpm, aes(x=Condition,y=value,fill=Age)) +
  geom_boxplot() 

input_tpm %>% filter(Condition != "Nafa") %>%
  ggplot(aes (x = Condition, y = value, fill = Condition)) + 
  geom_boxplot(alpha = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("log2(TPM)")

zscore <- data.frame(t(scale(t(viral_tpm))))

agegroup <- ifelse(grepl(".CH", colnames(zscore)), "Child", 
                   ifelse(grepl("YA", colnames(zscore)), "Young Adult","Older Adult"))
condition <- ifelse(grepl("V2.",  colnames(zscore)), "Virus 0.2", 
                    ifelse(grepl("V6.",  colnames(zscore)), "Virus 0.6","Nafa"))


Heatmap(as.matrix(zscore),
        cluster_rows = F, 
        cluster_columns = T, 
        col = viridis(10), 
        top_annotation = HeatmapAnnotation(AgeGroup = agegroup, annotation_name_side = "left",
                                           gp = gpar(col = "black"),
                                           show_legend = T),
        bottom_annotation = HeatmapAnnotation(Condition = condition, annotation_name_side = "left",
                                              gp = gpar(col = "black"),
                                              show_legend = T),
        heatmap_legend_param = list(
          title = "Gene Expression Z-score",
          direction = "horizontal",
          title_position = "topleft",
          at = c(-2,0,2),
          legend_width = unit(4, "cm")))



# Heatmap Plots - Nafa & 0.6 MOI Figure S3E ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm %>% select(contains("V6") | contains("Nafa") | contains("V2.")) %>%
  select(contains(".BM.") | contains(".HK.") | contains(".LT.") | contains("HL") |
           contains("021CsH") | contains(".MB.") | contains(".SLW.") |
           contains("RB") | contains("PT") | contains("RG") | contains("AT"))

heatmap_average2 <- function(genes,mean_only){
  heatmap_plot <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
  average_plot <- data.frame(row.names = row.names(heatmap_plot))
  
  average_plot <- tpm %>%
    filter(row.names(tpm) %in% genes) %>%
    mutate(Virus0.2_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
    mutate(Virus0.6_OA = rowMeans(select(.,contains("V6.OA.")))) %>%
    mutate(Virus0.2_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
    mutate(Virus0.6_YA = rowMeans(select(.,contains("V6.YA.")))) %>%
    mutate(Virus0.2_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
    mutate(Virus0.6_CH = rowMeans(select(.,contains("V6.CH.")))) %>%
    select(Virus0.2_OA:Virus0.6_CH)
  
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
            top_annotation = HeatmapAnnotation("Age Group" = c("OA","OA","OA",
                                                               "YA","YA","YA",
                                                               "CH","CH","CH"), 
                                               annotation_name_side = "left", 
                                               col = list("Age Group" = c("OA"="#BEBADA","OA"="#BEBADA","OA"="#BEBADA",
                                                                         "YA"="#FFFFB3", "YA"="#FFFFB3", "YA"="#FFFFB3",
                                                                         "CH"="#8DD3C7", "CH"="#8DD3C7", "CH"="#8DD3C7")),
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
            top_annotation = HeatmapAnnotation("Age Group" = c("OA","OA","OA",
                                                               "YA","YA","YA",
                                                               "CH","CH","CH"), 
                                               annotation_name_side = "left", 
                                               col = list("Age Group" = c("OA"="#BEBADA","OA"="#BEBADA","OA"="#BEBADA",
                                                                          "YA"="#FFFFB3", "YA"="#FFFFB3", "YA"="#FFFFB3",
                                                                          "CH"="#8DD3C7", "CH"="#8DD3C7", "CH"="#8DD3C7")),
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

  input_tpm <- data.frame((subset(tpm, rownames(tpm) %in% genes))) #Subset only genes in list
  zscore <- data.frame(t(scale(t(input_tpm))))
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

# IFN, cytokine heatmaps
heatmap_average2(genes = c("IFNB1","IFNE","IFNL1","IFNL2","IFNL3",
                          "IRF1","IRF2","IRF3","IRF5","IRF6","IRF7","IRF9",
                          "IFNAR1","IFNAR2","IFNGR1","IFNGR2","IFNLR1"),mean_only = F)
wilcoxon_test2(genes=c("IFNAR1","IFNAR2","IFNGR1","IFNGR2","IFNLR1"))

heatmap_average2(genes = c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                          "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                          "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                          "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"), mean_only = F)
wilcoxon_test2(genes=c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                       "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                       "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                       "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"))
wilcoxon_test2(genes=c("CCL1","CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8","CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
                       "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", 
                       "CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "CXCL6", "CXCL9", "CXCL10", "CXCL11", 
                       "CXCL12", "CXCL13", "CXCL14", "CX3CL1", "CXCL16", "CXCL17"))

heatmap_average2(genes = c("IL1A","IL1B","IL6","IL7","CXCL8","IL11","IL12A","IL15","IL16","IL17C","IL18",
                          "IL18","IL19","IL23A","IL32","IL33","IL34","IL36G","IL37"), mean_only = F)

heatmap_average2(genes = c("DDX58","DHX58","IFIH1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6"),mean_only = F)
wilcoxon_test2(gen = c("DDX58","DHX58","IFIH1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6"))

heatmap_average2(genes = c("ACE","ACE2","ANPEP","DPP4","ATP1A1","ATP1B1","NRP1",
                          "CTSL","FURIN","TMPRSS2","TMPRSS4"),mean_only = F)
wilcoxon_test2(genes = c("ACE","ACE2","ANPEP","DPP4","ATP1A1","ATP1B1","NRP1",
                        "CTSL","FURIN","TMPRSS2","TMPRSS4"))

heatmap_average2(genes = c("FAP","H19"), mean_only = F)

# Nafa & 0.6 Virus - Replication complex ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm %>% select(contains("Nafa") | contains("V2.")) %>%
  select(contains(".BM.") | contains(".HK.") | contains(".LT.") | contains("HL") |
           contains("021CsH") | contains(".MB.") | contains(".SLW.") |
           contains("RB") | contains("PT") | contains("RG") | contains("AT"))

interactome_genes <- read.table(file = "interactome_genes.txt", sep = "\t", header = T)

interactome_tpm <- data.frame((subset(tpm, rownames(tpm) %in% interactome_genes$Protein))) #Subset only genes in list
interactome_tpm <- interactome_tpm[complete.cases(interactome_tpm),] #Remove NAs

zscore <- data.frame(t(scale(t(interactome_tpm))))
zscore$geneID <- row.names(zscore)

zscore_long <- reshape2::melt(zscore) #wide-form to long-form
zscore_long <- zscore_long[order(zscore_long$gene),] #order by gene

zscore_long$Condition <- ifelse(grepl("V2.", zscore_long$variable), "Virus 0.2", "Nafa")
zscore_long$Age <- ifelse(grepl("CH", zscore_long$variable), "CH", 
                          ifelse(grepl("YA", zscore_long$variable), "YA","OA"))
zscore_long$Age <- factor(x = zscore_long$Age, levels = c("OA","YA", "CH")) #Keeps plot in order of child > YA > OA
zscore_long$Group <- paste0(zscore_long$Condition,"_",zscore_long$Age)
zscore_long$Condition <- factor(x=zscore_long$Condition, levels=c("Virus 0.2","Nafa"))
zscore_long <- merge(zscore_long, interactome_genes,  by.x = "geneID", by.y= "Protein")

#Nsp7, 8, 12, 13 (Replication complex)

stat_test <- subset(zscore_long, CoV2 %in% c("nsp7","nsp8","nsp12","nsp13")) %>%
  group_by(CoV2) %>%
  wilcox_test(value ~ Group, p.adjust.method = "bonferroni")
stat_test

ggplot(subset(zscore_long,CoV2 %in% c("nsp7","nsp8","nsp12","nsp13")), aes(x = Condition, y = value, fill = Age)) + 
  geom_violin() +
  theme_bw() +
  ylim(-4,4) +
  scale_fill_manual(values = rep(c("#BEBADA","#FFFFB3","#8DD3C7"),2)) +
  theme(legend.position = "none") +
  facet_wrap(~CoV2, ncol = 5,
             scales = "fixed") +
  ylab("Z-score of gene expression") + 
  stat_summary(fun.y=mean, geom="point", size=2, color="black", position=position_dodge(1))

#Replication Complex

rep <- subset(zscore_long, Type %in% c("Replication complex"))

stat_test <- rep %>%
  group_by(Condition) %>%
  wilcox_test(value ~ Age, p.adjust.method = "bonferroni")
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

# Nafa & 0.6 Virus Mean Heatmap #### 
tpm <- read.table(file = "viral_tpm_mockinc.txt", sep = "\t", header = T, row.names = 1)
tpm <- tpm %>% select(contains("V6") | contains("Nafa") | contains("V2.")) #%>%
  select(contains(".BM.") | contains(".HK.") | contains(".LT.") | contains("HL") |
           contains("021CsH") | contains(".MB.") | contains(".SLW.") |
           contains("RB") | contains("PT") | contains("RG") | contains("AT"))

virusgenes <- paste("GU280_gp", str_pad(1:11,2,pad=0),sep = "")

heatmap_average2(genes = virusgenes, mean_only = T)
wilcoxon_test2(genes = virusgenes)

#NSPs - Replication

nsps <- interactome_genes %>%
  subset(CoV2 %in% c("nsp8","nsp12")) %>%
  select(Protein)

heatmap_average2(genes = nsps$Protein,mean_only = T)
wilcoxon_test2(genes = nsps$Protein)


# Nata & 0.6 ISG cluster plots - Figure S3F ####
tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm %>% select(contains("V6") | contains("Nafa") | contains("V2.")) %>%
  select(contains(".BM.") | contains(".HK.") | contains(".LT.") | contains("HL") |
           contains("021CsH") | contains(".MB.") | contains(".SLW.") |
           contains("RB") | contains("PT") | contains("RG") | contains("AT"))

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
  
  C12_mean <- tpm %>%
    filter(row.names(tpm) %in% C12) %>%
    mutate(Virus0.2_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
    mutate(Virus0.6_OA = rowMeans(select(.,contains("V6.OA.")))) %>%
    
    mutate(Virus0.2_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
    mutate(Virus0.6_YA = rowMeans(select(.,contains("V6.YA.")))) %>%
    
    mutate(Virus0.2_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
    mutate(Virus0.6_CH = rowMeans(select(.,contains("V6.CH.")))) %>%
    select(Virus0.2_OA:Virus0.6_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C1/2 Mean")
  
  C3_mean <- tpm %>%
    filter(row.names(tpm) %in% C3) %>%
    mutate(Virus0.2_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
    mutate(Virus0.6_OA = rowMeans(select(.,contains("V6.OA.")))) %>%
    
    mutate(Virus0.2_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
    mutate(Virus0.6_YA = rowMeans(select(.,contains("V6.YA.")))) %>%
    
    mutate(Virus0.2_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
    mutate(Virus0.6_CH = rowMeans(select(.,contains("V6.CH.")))) %>%
    select(Virus0.2_OA:Virus0.6_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C3 Mean")
  
  C4_mean <- tpm %>%
    filter(row.names(tpm) %in% C4) %>%
    mutate(Mock_OA = rowMeans(select(.,contains("M.OA.")))) %>%
    mutate(Virus0.2_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
    mutate(Virus0.6_OA = rowMeans(select(.,contains("V6.OA.")))) %>%
    
    mutate(Virus0.2_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
    mutate(Virus0.6_YA = rowMeans(select(.,contains("V6.YA.")))) %>%
    
    mutate(Virus0.2_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
    mutate(Virus0.6_CH = rowMeans(select(.,contains("V6.CH.")))) %>%
    select(Virus0.2_OA:Virus0.6_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C4 Mean")
  
  C5_mean <- tpm %>%
    filter(row.names(tpm) %in% C5) %>%
    mutate(Virus0.2_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
    mutate(Virus0.6_OA = rowMeans(select(.,contains("V6.OA.")))) %>%
    
    mutate(Virus0.2_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
    mutate(Virus0.6_YA = rowMeans(select(.,contains("V6.YA.")))) %>%
    
    mutate(Virus0.2_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
    mutate(Virus0.6_CH = rowMeans(select(.,contains("V6.CH.")))) %>%
    select(Virus0.2_OA:Virus0.6_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C5 Mean")
  
  C35_mean <- tpm %>%
    filter(row.names(tpm) %in% C35) %>%
    mutate(Virus0.2_OA = rowMeans(select(.,contains("V2.OA.")))) %>%
    mutate(Nafa_OA = rowMeans(select(.,contains("Nafa.OA.")))) %>%
    mutate(Virus0.6_OA = rowMeans(select(.,contains("V6.OA.")))) %>%
    
    mutate(Virus0.2_YA = rowMeans(select(.,contains("V2.YA.")))) %>%
    mutate(Nafa_YA = rowMeans(select(.,contains("Nafa.YA.")))) %>%
    mutate(Virus0.6_YA = rowMeans(select(.,contains("V6.YA.")))) %>%
    
    mutate(Virus0.2_CH = rowMeans(select(.,contains("V2.CH.")))) %>%
    mutate(Nafa_CH = rowMeans(select(.,contains("Nafa.CH.")))) %>%
    mutate(Virus0.6_CH = rowMeans(select(.,contains("V6.CH.")))) %>%
    select(Virus0.2_OA:Virus0.6_CH) %>%
    summarise_all(.,mean) %>%
    set_rownames("C3-5 Mean")
  
  cluster_means <- bind_rows(C12_mean,C3_mean,C4_mean,C5_mean,C35_mean)
  
  zscore <- data.frame(t(scale(t(cluster_means))))
  
  col_fun = colorRamp2(c(-1.5,-0.75,0,0.75,1.5), 
                       c("#440154FF","#3B528BFF","#21908CFF","#5DC863FF","#FDE725FF"))
  
  Heatmap(as.matrix(zscore),
          cluster_rows = F, 
          cluster_columns = F, 
          col = col_fun, 
          top_annotation = HeatmapAnnotation("Age Group" = c("OA","OA","OA",
                                                             "YA","YA","YA",
                                                             "CH","CH","CH"), 
                                             annotation_name_side = "left", 
                                             col = list("Age Group" = c("OA"="#BEBADA","OA"="#BEBADA","OA"="#BEBADA",
                                                                        "YA"="#FFFFB3", "YA"="#FFFFB3", "YA"="#FFFFB3",
                                                                        "CH"="#8DD3C7", "CH"="#8DD3C7", "CH"="#8DD3C7")),
                                             gp = gpar(col = "black"),
                                             show_legend = F),
          heatmap_legend_param = list(
            title = "Gene Expression Z-score",
            direction = "horizontal",
            title_position = "topleft",
            at = c(-2,-1,0,1,2),
            legend_width = unit(4, "cm")))
  
}

#Stats test 
#C1
wilcoxon_test2(genes = isgs$Gene[isgs$Cluster == "1" | isgs$Cluster == "2"])
wilcoxon_test2(genes = isgs$Gene[isgs$Cluster == "3"])
wilcoxon_test2(genes = isgs$Gene[isgs$Cluster == "4"])
wilcoxon_test2(genes = isgs$Gene[isgs$Cluster == "5"])
wilcoxon_test2(genes = isgs$Gene[isgs$Cluster == "3" | isgs$Cluster == "4" | isgs$Cluster == "5"])

# Nafa & 0.6 MOI PCoA Analysis - Figure S3K ####

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
tpm <- tpm %>% select(contains("V6") | contains("Nafa") | contains("V2.")) %>%
  select(contains(".BM.") | contains(".HK.") | contains(".LT.") | contains("HL") |
           contains("021CsH") | contains(".MB.") | contains(".SLW.") |
           contains("RB") | contains("PT") | contains("RG") | contains("AT"))

sinfo <- read.table(file = "sinfo.txt", sep = "\t", header = T)
viral_percent <- read.table(file = "viral_percent_incMock.txt", sep = "\t", header = T)

group <- data.frame(names(tpm))

group <- sinfo %>%
  select(sname,individualid,agegroup,treatment,technicalrepid,tissue) %>%
  mutate_at("sname",str_replace_all,"-",".") %>%
  merge(group,.,by.x="names.tpm.",by.y="sname") %>%
  mutate(group,Condition = paste0(agegroup,"_",treatment))

group <- merge(group,viral_percent,by.x="names.tpm.",by.y="sname")

group <- group[ order(group$names.tpm.), ]
mds <- plotMDS(tpm,plot = T,gene.selection = "common", top = 20000)

mds_plot <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  set_colnames(c("Dim1", "Dim2")) %>%
  rownames_to_column("SampleID")

gg <- merge(mds_plot,group,by.x="SampleID",by.y="names.tpm.")
gg$agegroup <- factor(gg$agegroup, levels = c("OA","YA","CH"))
gg <- gg %>%
  filter(!treatment == "Mock") %>%
  filter(tissue == "N")

ggplot(gg, aes(x = Dim1, y = Dim2, colour = agegroup, fill = agegroup)) +
  geom_point(size = 5,aes(shape=treatment)) + 
  scale_colour_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  scale_fill_manual(values=c("#BEBADA","#FFFFB3","#8DD3C7")) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_text(aes(label=paste(individualid,treatment,tissue,sep=" ")),hjust=0, vjust=0)

#Vegan ADONIS

counts <- tpm %>%
  select(gg$SampleID)

adon.results<-adonis(t(tpm) ~ gg$agegroup + gg$treatment + gg$percentage,perm=999)
print(adon.results)

