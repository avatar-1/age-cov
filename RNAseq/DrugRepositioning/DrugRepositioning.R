library(data.table)
library(igraph)
library(stringr)
library(tidyverse)
library(magrittr)
library(circlize)
library(ComplexHeatmap)

# Load in files - Transcriptome ####

drugs <- read.csv("./input/drugs.csv", header = T)
sars_human_ppi <- read.csv("./input/sars-cov2_human_ppi.csv", header = T)
disease <- unique(sars_human_ppi$Human)
human_ppi <- read.csv("./input/human_ppi.csv", header = T)
prot2name_map <- unique(rbind(data.frame(UniProt.ID = human_ppi$uniprot1, Gene = human_ppi$symbol1),
                              data.frame(UniProt.ID = human_ppi$uniprot2, Gene = human_ppi$symbol2)))
DrugBank <- read.csv("./input/uniprot links.csv", header = T)
db <- merge(as.data.table(DrugBank),as.data.table(prot2name_map), by.x = "UniProt.ID", by.y = "UniProt.ID", all.x = TRUE, all.y = F)

tpm <- read.table(file = "tpm.txt", sep = "\t", header = T)
expr <- tpm %>%
  select(contains("V2")) %>%
  select(-contains(".B."))
expr <- 2^expr #un-log tpm

# Proximity measure & Z-score functions - Transcriptome ####

getShortestDistance <- function(drug, net, targets, disease, iter = 100){
  print(drug)
  
  d <- disease
  keep <- which(d %in% V(net)$name)
  d <- unique(d[keep])
  
  t <-  targets[which(targets$DrugBank.ID %in% drug),]$Gene
  keep <- which(t %in% V(net)$name)
  t <- unique(t[keep])

  if(length(t) > 0){
    tmp <- shortest.paths(net, v = as.character(unique(t)), to=unique(d))
    tmp[is.infinite(tmp)] <- NA
    if(!is.null(dim(tmp))){
      score <- mean(c(colMeans(tmp, na.rm = T), rowMeans(tmp, na.rm = T)), na.rm = T)
      
    }else{
      score <- mean(tmp, na.rm = T)
    }
    z <- getzScore(drug, net, length(t), length(d), t, d, iter = iter, s = score)
    
    return(c(score, z, pnorm(z), length(t), mean(degree(net, v = as.character((t))))))
  }else{
    return(c(NA, NA, NA, NA, NA))
  }
  
}

getzScore <- function (drug, net, tl, dl, tr, de, iter, s){
  null <- c()
  for(i in 1:iter){
    t <- names(sample(V(net), tl))
    d <- names(sample(V(net), dl))
    
    tmp <- shortest.paths(net, v = t, to=d)
    tmp[is.infinite(tmp)] <- NA
    if(!is.null(dim(tmp))){
      score <- mean(c(colMeans(tmp, na.rm = T), rowMeans(tmp, na.rm = T)), na.rm = T)
      
    }else{
      score <- mean(tmp, na.rm = T)
    }

    null[i] <- score
    
  }

  z <- (s - mean(null, na.rm = T))/sd(null, na.rm = T)
  return(z)
}

# All ages - Transcriptome####

expr.all <- expr %>%
  rowMeans(.)
all <- names(which(expr.all > mean(expr.all)))
ppi.all <- subset(human_ppi, human_ppi$symbol1 %in% all & human_ppi$symbol2 %in% all)

dt.a <- subset(db, db$Gene %in% all)
drugs.a <- unique(dt.a$DrugBank.ID)

net.all <- graph_from_data_frame(ppi.all[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

all.dist <- sapply(drugs.a, getShortestDistance, net.all, dt.a, disease, iter = 100)
all.dist <- all.dist %>%
  t(.) %>%
  set_rownames(drugs.a)
colnames(all.dist) <- c("mean_proximity", "z-score", "p-value", "No.of.Targets", "Mean_degrees")
write.csv(all.dist, file = "./output/All.csv", quote = F)



# Child - Transcriptome ####

expr.ch <- expr %>%
  select(contains(".CH.")) %>%
  rowMeans(.)
ch <- names(which(expr.ch > mean(expr.ch)))
ppi.ch <- subset(human_ppi, human_ppi$symbol1 %in% ch & human_ppi$symbol2 %in% ch)

dt.c <- subset(db, db$Gene %in% ch)
drugs.c <- unique(dt.c$DrugBank.ID)

net.ch <- graph_from_data_frame(ppi.ch[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

child.dist <- sapply(drugs.c, getShortestDistance, net.ch, dt.c, disease, iter = 100)
child.dist <- child.dist %>%
  t(.) %>%
  set_rownames(drugs.c)
colnames(child.dist) <- c("mean_proximity", "z-score", "p-value", "No.of.Targets", "Mean_degrees")
write.csv(child.dist, file = "./output/Child.csv", quote = F)

# Young Adult - Transcriptome ####

expr.ya <- expr %>%
  select(contains(".YA."))%>%
  rowMeans(.)
ya <- names(which(expr.ya > mean(expr.ya)))
ppi.ya <- subset(human_ppi, human_ppi$symbol1 %in% ya & human_ppi$symbol2 %in% ya)
dt.y <- subset(db, db$Gene %in% ya)
drugs.y <- unique(dt.y$DrugBank.ID)

net.ya <- graph_from_data_frame(ppi.ya[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

young.dist <- sapply(drugs.y, getShortestDistance, net.ya, dt.y, disease, iter = 100)
young.dist <- young.dist %>%
  t(.) %>%
  set_rownames(drugs.y)
colnames(young.dist) <- c("mean_proximity", "z-score", "p-value", "No.of.Targets", "Mean_degrees")
write.csv(young.dist, file = "./output/YoungerAdult.csv", quote = F)

# Older Adult - Transcriptome ####

expr.oa <- expr %>%
  select(contains(".OA."))%>%
  rowMeans(.)

oa <- names(which(expr.oa > mean(expr.oa)))
ppi.oa <- subset(human_ppi, human_ppi$symbol1 %in% oa & human_ppi$symbol2 %in% oa)
dt.o <- subset(db, db$Gene %in% oa)
drugs.o <- unique(dt.o$DrugBank.ID)

net.oa <- graph_from_data_frame(ppi.oa[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

old.dist   <- sapply(drugs.o, getShortestDistance, net.oa, dt.o, disease, iter = 100)
old.dist <- old.dist %>%
  t(.) %>%
  set_rownames(drugs.o)
colnames(old.dist) <- c("mean_proximity", "z-score", "p-value", "No.of.Targets", "Mean_degrees")
write.csv(old.dist, file = "./output/OlderAdult2.csv", quote = F)
# Load in files - Proteomics ####

drugs <- read.csv("./input/drugs.csv", header = T)
sars_human_ppi <- read.csv("./input/sars-cov2_human_ppi.csv", header = T)
disease <- unique(sars_human_ppi$Human)
human_ppi <- read.csv("./input/human_ppi.csv", header = T)
prot2name_map <- unique(rbind(data.frame(UniProt.ID = human_ppi$uniprot1, Gene = human_ppi$symbol1),
                              data.frame(UniProt.ID = human_ppi$uniprot2, Gene = human_ppi$symbol2)))
DrugBank <- read.csv("./input/uniprot links.csv", header = T)
db <- merge(as.data.table(DrugBank),as.data.table(prot2name_map), by.x = "UniProt.ID", by.y = "UniProt.ID", all.x = TRUE, all.y = F)

log2lfq <- read.table(file="log2_lfq.txt", sep = "\t", header=T)
expr <- log2lfq %>%
  select(contains("V2")) %>%
  select(-contains(".B."))
expr <- 2^expr #un-log tpm

# Proximity measure & Z-score functions - Proteomics####

getShortestDistance <- function(drug, net, targets, disease, iter = 100){
  print(drug)
  
  d <- disease
  keep <- which(d %in% V(net)$name)
  d <- unique(d[keep])
  
  t <-  targets[which(targets$DrugBank.ID %in% drug),]$Gene
  keep <- which(t %in% V(net)$name)
  t <- unique(t[keep])
  
  if(length(t) > 0){
    tmp <- shortest.paths(net, v = as.character(unique(t)), to=unique(d))
    tmp[is.infinite(tmp)] <- NA
    if(!is.null(dim(tmp))){
      score <- mean(c(colMeans(tmp, na.rm = T), rowMeans(tmp, na.rm = T)), na.rm = T)
      
    }else{
      score <- mean(tmp, na.rm = T)
    }
    z <- getzScore(drug, net, length(t), length(d), t, d, iter = iter, s = score)
    
    return(c(score, z, pnorm(z), length(t), mean(degree(net, v = as.character((t))))))
  }else{
    return(c(NA, NA, NA, NA, NA))
  }
  
}

getzScore <- function (drug, net, tl, dl, tr, de, iter, s){
  null <- c()
  for(i in 1:iter){
    t <- names(sample(V(net), tl))
    d <- names(sample(V(net), dl))
    
    tmp <- shortest.paths(net, v = t, to=d)
    tmp[is.infinite(tmp)] <- NA
    if(!is.null(dim(tmp))){
      score <- mean(c(colMeans(tmp, na.rm = T), rowMeans(tmp, na.rm = T)), na.rm = T)
      
    }else{
      score <- mean(tmp, na.rm = T)
    }
    
    null[i] <- score
    
  }
  
  z <- (s - mean(null, na.rm = T))/sd(null, na.rm = T)
  return(z)
}

# All ages - Proteomics####

expr.all <- expr %>%
  rowMeans(.)
all <- names(which(expr.all > mean(expr.all)))
ppi.all <- subset(human_ppi, human_ppi$symbol1 %in% all & human_ppi$symbol2 %in% all)

dt.a <- subset(db, db$Gene %in% all)
drugs.a <- unique(dt.a$DrugBank.ID)

net.all <- graph_from_data_frame(ppi.all[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

all.dist <- sapply(drugs.a, getShortestDistance, net.all, dt.a, disease, iter = 100)
all.dist <- all.dist %>%
  t(.) %>%
  set_rownames(drugs.a)
colnames(all.dist) <- c("mean_proximity", "z-score", "p-value", "No.of.Targets", "Mean_degrees")
write.csv(all.dist, file = "./output/Proteomics/All.csv", quote = F)

# Child - Proteomics####

expr.ch <- expr %>%
  select(contains(".CH.")) %>%
  rowMeans(.)
ch <- names(which(expr.ch > mean(expr.ch)))
ppi.ch <- subset(human_ppi, human_ppi$symbol1 %in% ch & human_ppi$symbol2 %in% ch)

dt.c <- subset(db, db$Gene %in% ch)
drugs.c <- unique(dt.c$DrugBank.ID)

net.ch <- graph_from_data_frame(ppi.ch[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

child.dist <- sapply(drugs.c, getShortestDistance, net.ch, dt.c, disease, iter = 100)
child.dist <- child.dist %>%
  t(.) %>%
  set_rownames(drugs.c)
colnames(child.dist) <- c("mean_proximity", "z-score", "p-value", "No.of.Targets", "Mean_degrees")
write.csv(child.dist, file = "./output/Proteomics/Child.csv", quote = F)

# Young Adult - Proteomics####

expr.ya <- expr %>%
  select(contains(".YA."))%>%
  rowMeans(.)
ya <- names(which(expr.ya > mean(expr.ya)))
ppi.ya <- subset(human_ppi, human_ppi$symbol1 %in% ya & human_ppi$symbol2 %in% ya)
dt.y <- subset(db, db$Gene %in% ya)
drugs.y <- unique(dt.y$DrugBank.ID)

net.ya <- graph_from_data_frame(ppi.ya[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

young.dist <- sapply(drugs.y, getShortestDistance, net.ya, dt.y, disease, iter = 100)
young.dist <- young.dist %>%
  t(.) %>%
  set_rownames(drugs.y)
colnames(young.dist) <- c("mean_proximity", "z-score", "p-value", "No.of.Targets", "Mean_degrees")
write.csv(young.dist, file = "./output/Proteomics/YoungerAdult.csv", quote = F)

# Older Adult - Proteomics####

expr.oa <- expr %>%
  select(contains(".OA."))%>%
  rowMeans(.)

oa <- names(which(expr.oa > mean(expr.oa)))
ppi.oa <- subset(human_ppi, human_ppi$symbol1 %in% oa & human_ppi$symbol2 %in% oa)
dt.o <- subset(db, db$Gene %in% oa)
drugs.o <- unique(dt.o$DrugBank.ID)

net.oa <- graph_from_data_frame(ppi.oa[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

old.dist   <- sapply(drugs.o, getShortestDistance, net.oa, dt.o, disease, iter = 100)
old.dist <- old.dist %>%
  t(.) %>%
  set_rownames(drugs.o)
colnames(old.dist) <- c("mean_proximity", "z-score", "p-value", "No.of.Targets", "Mean_degrees")
write.csv(old.dist, file = "./output/Proteomics/OlderAdult2.csv", quote = F)




# Drug class bubble plot ####

plot <- read.table(file = "./Plots/drugs_plots_combined.txt", sep = "\t", header = T)
plot <- read.table(file = "./Plots/drugs_plots_combined.txt", sep = "\t", header = T)
plot$AgeGroup <- factor(x = plot$AgeGroup, levels = c("All","OA", "YA", "CH"))

filt_plot <- plot %>%
  filter(Therap_class != "Adrenergic and dopaminergic agents" &
           Therap_class != "Amino acids" & 
           Therap_class != "Antibiotics" &
           Therap_class != "Antioxidants" & 
           Therap_class != "Bisphosphonates" & 
           Therap_class != "Fungicide" & 
           Therap_class != "Insulin" & 
           Therap_class != "Lipids" & 
           Therap_class != "Nucleotides" &
           Therap_class != "Plasminogen activator" & 
           Therap_class != "Pyrans" & 
           Therap_class != "Sulfur compounds" & 
           Therap_class != "Unknown") %>%
  mutate(Type = factor(.$Type, levels = c("RNA","Protein")))

ggplot(filt_plot, aes(x = AgeGroup, y = Therap_class, color = AgeGroup)) +
  geom_count()  +
  scale_size(range = c(1,10), breaks = seq(2,20,2),name = "Number of drugs") +
  scale_color_manual(values = c("#E896C1","#BEBADA","#F8F5B7","#8DD3C7")) +
  theme_bw() +
  scale_y_discrete(limits=rev) +
  facet_wrap(~Type)#4,5x5,5

# Mean z-score plots ####

mean_zscore_rna <- filt_plot %>%
  filter(Type == "RNA") %>%
  group_by(AgeGroup) %>%
  group_by(Therap_class, add = TRUE) %>%
  summarise(count = mean(z.score)) %>%
  spread(AgeGroup, count) %>%
  set_rownames(.$Therap_class) %>%
  mutate(Therap_class = NULL)

mean_zscore_prot <- filt_plot %>%
  filter(Type == "Protein") %>%
  group_by(AgeGroup) %>%
  group_by(Therap_class, add = TRUE) %>%
  summarise(count = mean(z.score)) %>%
  spread(AgeGroup, count) %>%
  set_rownames(.$Therap_class) %>%
  mutate(Therap_class = NULL)

col_fun = colorRamp2(c(-14,-11,-8,-5,-2.5), 
                     c("#440154FF","#3B528BFF","#21908CFF","#5DC863FF","#FDE725FF"))

Heatmap(as.matrix(mean_zscore_rna), 
        cluster_rows = F, 
        cluster_columns = F, 
        col=col_fun, 
        top_annotation = HeatmapAnnotation("Age Group" = c("All","OA","YA","CH"), 
                                           annotation_name_side = "left", 
                                           col = list("Age Group" = c("All"="#E896C1","OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7")),
                                           gp = gpar(col = "black"),
                                           show_legend = F),
        heatmap_legend_param = list(
          title = "Z-score",
          direction = "horizontal",
          title_position = "topleft",
          at = c(-15,-10,-5,0),
          legend_width = unit(4, "cm")))

Heatmap(as.matrix(mean_zscore_prot), 
        cluster_rows = F, 
        cluster_columns = F, 
        col=col_fun, 
        top_annotation = HeatmapAnnotation("Age Group" = c("All","OA","YA","CH"), 
                                           annotation_name_side = "left", 
                                           col = list("Age Group" = c("All"="#E896C1","OA"="#BEBADA","YA"="#FFFFB3","CH"="#8DD3C7")),
                                           gp = gpar(col = "black"),
                                           show_legend = F),
        heatmap_legend_param = list(
          title = "Z-score",
          direction = "horizontal",
          title_position = "topleft",
          at = c(-15,-10,-5,0),
          legend_width = unit(4, "cm")))


# Cytoscape file manipulation ####

drug_targets <- read.table(file="D:/OneDrive - UNSW/miCF/SARS-CoV2_Study/DrugRepositioning/Cytoscape/drug_targets.txt",sep="\t",header=T)
sars_human_ppi <- read.table(file="D:/OneDrive - UNSW/miCF/SARS-CoV2_Study/DrugRepositioning/Cytoscape/sars-cov2_human.txt",sep="\t",header=T)
sig <- read.table(file="D:/OneDrive - UNSW/miCF/SARS-CoV2_Study/DrugRepositioning/Cytoscape/DrugsSignificantList.txt",sep="\t",header=T)

sig_targets <- drug_targets %>%
  filter(DrugBank.ID %in% sig$dbID)

sig_sars <- sars_human_ppi %>%
  filter(Human %in% sig_targets$Protein.name)

write.table(sig_targets, file="D:/OneDrive - UNSW/miCF/SARS-CoV2_Study/DrugRepositioning/Cytoscape/sig_targets.txt",sep="\t",quote=F,row.names = F)
write.table(sig_sars, file="D:/OneDrive - UNSW/miCF/SARS-CoV2_Study/DrugRepositioning/Cytoscape/sig_sars.txt",sep="\t",quote=F,row.names=F)
