#Code is used for correlation analyses of RNA-seq data for COVID-19 2020 project
#Code developed by Alexander Capraro, Septmeber 30 2020

library(tidyverse)
library(dplyr)
library(magrittr)
library(vegan)
library(ggplot2)
library(grid)
library(edgeR)
library(RColorBrewer)
library(usedist)

#### Mock + Virus  RNA-seq ####

tpm <- read.table(file = "log_cpm_1.txt", sep = "\t", header = T) #change to tpm.txt for MDS plot
tpm_plots <- tpm[,c(1:16,18:40,42:44)] #Remove BC_B1

#EdgeR MDS plot 
groups <- c("Bronchial_Mock","Bronchial_Mock","Bronchial_Mock","Bronchial_Mock","Bronchial_Mock",
            "Child_Mock","Child_Mock","Child_Mock","Child_Mock","Child_Mock",
            "Older Adult_Mock","Older Adult_Mock","Older Adult_Mock","Older Adult_Mock","Older Adult_Mock",
            "Young Adult_Mock","Young Adult_Mock","Young Adult_Mock","Young Adult_Mock",
            "Bronchial_Virus","Bronchial_Virus","Bronchial_Virus","Bronchial_Virus","Bronchial_Virus",
            "Child_Virus","Child_Virus","Child_Virus","Child_Virus","Child_Virus","Child_Virus","Child_Virus","Child_Virus",
            "Older Adult_Virus","Older Adult_Virus","Older Adult_Virus","Older Adult_Virus","Older Adult_Virus","Older Adult_Virus",
            "Young Adult_Virus","Young Adult_Virus","Young Adult_Virus","Young Adult_Virus")

mds <- plotMDS(tpm_plots,plot = T,col=c(rep("black",5),rep("red",5),rep("blue",5),rep("green",4),rep("grey",5),rep("orange",8),rep("purple",6),rep("yellow",4)),gene.selection = "common", top = 20000)

mds_plot <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  set_colnames(c("Dim1", "Dim2")) %>%
  rownames_to_column("SampleID") %>%
  add_column(groups)

gg <- merge(mds_plot,aggregate(cbind(mean.x=Dim1,mean.y=Dim2)~groups,mds_plot,mean),by="groups")

ggplot(gg, aes(x = Dim1, y = Dim2, colour = groups, fill = groups)) +
  geom_point(size = 3) + 
  scale_colour_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  ylim(-0.3,0.4) +
  xlim(-0.4,0.8)
ggplot(gg, aes(x = Dim1, y = Dim2, colour = groups, fill = groups)) +
  geom_point(size = 3) + 
  scale_colour_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  geom_point(aes(x=mean.x,y=mean.y),size=5) +
  geom_segment(aes(x=mean.x, y=mean.y, xend=Dim1, yend=Dim2)) +
  stat_ellipse(aes(fill = groups),level = 0.4,geom="polygon",alpha=0.2)#export 6x4.5

#Vegan ADONIS

age <- c("Child","Child","Child","Child","Child","Child","Child","Child","Child","Child","Older Adult","Older Adult","Older Adult","Older Adult","Older Adult","Young Adult","Young Adult","Young Adult","Young Adult",
         "Child","Child","Child","Child","Child","Child","Child","Child","Child","Child","Child","Child","Child","Older Adult","Older Adult","Older Adult","Older Adult","Older Adult","Older Adult","Young Adult","Young Adult","Young Adult","Young Adult")
treatment <- c("Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock",
               "Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus")
tissue <- c("B","B","B","B","B","N","N","N","N","N","N","N","N","N","N","N","N","N","N",
            "B","B","B","B","B","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N")
load <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          14,46,17,3.4,0.3,
          47,35,38,12,43,33,30,19,
          3,8,5,5,31,16,38,6,39,12)

adon.results<-adonis(t(tpm_plots) ~ age * treatment * load,perm=999)
#adon.results<-adonis(t(tpm_plots) ~ tissue, perm=999)
print(adon.results)

#betadispr

dst <- vegdist(t(tpm_plots))
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

dist_between_centroids(dst, c("Mock_CH_B_376AC_B2_t1","Mock_CH_B_AM_B2_t1","Mock_CH_B_IB_B2_t1","Mock_CH_B_NT_B2_t1","Mock_CH_B_TJ_B2_t1"),
                       c("Virus_CH_B_376AC_B2_t1","Virus_CH_B_AM_B2_t1","Virus_CH_B_IB_B2_t1","Virus_CH_B_NT_B2_t1","Virus_CH_B_TJ_B2_t1"))

dist_between_centroids(dst, c("Mock_CH_N_376AC_B2_t1","Mock_CH_N_AC_B1_t1","Mock_CH_N_IB_B2_t1","Mock_CH_N_NT_B1_t1","Mock_CH_N_TJ_B1_t1"), 
                       c("Virus_CH_N_376AC_B2_t1","Virus_CH_N_AC_B1_t1","Virus_CH_N_AC_B1_t2","Virus_CH_N_AM_B2_t2","Virus_CH_N_IB_B2_t1","Virus_CH_N_NT_B1_t1","Virus_CH_N_NT_B1_t2","Virus_CH_N_TJ_B1_t1"))


dist_between_centroids(dst,c("Mock_YA_N_7AC_B2_t1","Mock_YA_N_BC_B2_t1","Mock_YA_N_NA_B2_t1","Mock_YA_N_SB_B2_t1"), 
                       c("Virus_YA_N_7AC_B2_t1","Virus_YA_N_BC_B2_t1","Virus_YA_N_NA_B2_t1","Virus_YA_N_SB_B2_t1"))


dist_between_centroids(dst,c("Mock_OA_N_JW_B1_t1","Mock_OA_N_JW_B2_t1","Mock_OA_N_NiT_B1_t1","Mock_OA_N_NiT_B2_t1","Mock_OA_N_WB_B2_t1"),
                       c("Virus_OA_N_JW_B1_t1","Virus_OA_N_JW_B2_t1","Virus_OA_N_NiT_B1_t1","Virus_OA_N_NiT_B2_t1","Virus_OA_N_WB_B2_t1","Virus_OA_N_WB_B2_t2")) # Older adults

df <- data.frame(dst)
write.table(as.matrix(dst), file = "/Users/z3416833/Downloads/dst.txt", sep = "\t", quote = F)

#Vegan ADONIS

age <- c("Child","Child","Child","Child","Child","Older Adult","Older Adult","Older Adult","Older Adult","Older Adult","Young Adult","Young Adult","Young Adult","Young Adult",
         "Child","Child","Child","Child","Child","Child","Child","Child","Older Adult","Older Adult","Older Adult","Older Adult","Older Adult","Older Adult","Young Adult","Young Adult","Young Adult","Young Adult")
treatment <- c("Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock",
               "Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus")


adon.results<-adonis(t(tpm_plots) ~ treatment * age,perm=999)
print(adon.results)

#betadispr

dst <- vegdist(t(tpm_plots))
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

dist_between_centroids(dst, c("Mock_CH_N_376AC_B2_t1","Mock_CH_N_AC_B1_t1","Mock_CH_N_IB_B2_t1","Mock_CH_N_NT_B1_t1","Mock_CH_N_TJ_B1_t1"), 
                       c("Virus_CH_N_376AC_B2_t1","Virus_CH_N_AC_B1_t1","Virus_CH_N_AC_B1_t2","Virus_CH_N_AM_B2_t2","Virus_CH_N_IB_B2_t1","Virus_CH_N_NT_B1_t1","Virus_CH_N_NT_B1_t2","Virus_CH_N_TJ_B1_t1"))


dist_between_centroids(dst,c("Mock_YA_N_7AC_B2_t1","Mock_YA_N_BC_B2_t1","Mock_YA_N_NA_B2_t1","Mock_YA_N_SB_B2_t1"), 
                       c("Virus_YA_N_7AC_B2_t1","Virus_YA_N_BC_B2_t1","Virus_YA_N_NA_B2_t1","Virus_YA_N_SB_B2_t1"))


dist_between_centroids(dst,c("Mock_OA_N_JW_B1_t1","Mock_OA_N_JW_B2_t1","Mock_OA_N_NiT_B1_t1","Mock_OA_N_NiT_B2_t1","Mock_OA_N_WB_B2_t1"),
                       c("Virus_OA_N_JW_B1_t1","Virus_OA_N_JW_B2_t1","Virus_OA_N_NiT_B1_t1","Virus_OA_N_NiT_B2_t1","Virus_OA_N_WB_B2_t1","Virus_OA_N_WB_B2_t2")) # Older adults



df <- data.frame(dst)
write.table(as.matrix(dst), file = "dst.txt", sep = "\t", quote = F)


#### Mock + Virus  Proteome ####
prot <- read.table(file = "/Users/z3416833/Documents/miCF/SARS-CoV2_Study/Proteomics/6.10.20_Lysate/txt/log2_lfq_matrix.txt", sep = "\t", header = T, row.names = 1)

groups <- c("Mock_Bronchial","Mock_Bronchial",
            "Mock_Child","Mock_Child",
            "Mock_Older Adult","Mock_Older Adult","Mock_Older Adult",
            "Mock_Young Adult","Mock_Young Adult","Mock_Young Adult","Mock_Young Adult",
            "Virus_Bronchial","Virus_Bronchial","Virus_Bronchial",
            "Virus_Child","Virus_Child","Virus_Child",
            "Virus_Older Adult","Virus_Older Adult","Virus_Older Adult",
            "Virus_Young Adult","Virus_Young Adult","Virus_Young Adult","Virus_Young Adult")

mds <- plotMDS(prot,plot = T,col=c(rep("black",2),rep("red",2),rep("blue",3),rep("green",4),rep("grey",3),rep("orange",3),rep("purple",3),rep("yellow",4)), gene.selection = "common", top = 20000)

mds_plot <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  set_colnames(c("Dim1", "Dim2")) %>%
  rownames_to_column("SampleID") %>%
  add_column(groups)

ggplot(mds_plot, aes(x = Dim1, y = Dim2, colour = groups)) +
  geom_point(size = 3) + 
  scale_colour_brewer(palette = "Set3") +
  theme_bw()

#Vegan ADONIS

age <- c("Child","Child","Child","Child","Older Adult","Older Adult","Older Adult","Young Adult","Young Adult","Young Adult","Young Adult",
         "Child","Child","Child","Child","Child","Child","Older Adult","Older Adult","Older Adult","Young Adult","Young Adult","Young Adult","Young Adult")
treatment <- c("Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock","Mock",
               "Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus","Virus")

adon.results<-adonis(t(prot) ~ age * treatment ,perm=999)
print(adon.results)
