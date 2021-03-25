library(data.table)
library(igraph)
library(stringr)


# Proximity measure -------------------------------------------------------

getShortestDistance <- function(drug, net, targets, disease, iter = 100){
  print(drug)
  
  d <- disease
  keep <- which(d %in% V(net)$name)
  d <- unique(d[keep])
  
  t <-  targets[which(targets$DrugBank.ID %in% drug),]$Gene
  keep <- which(t %in% V(net)$name)
  t <- unique(t[keep])

  if(length(t) > 0){
    tmp <- shortest.paths(net, v = unique(t), to=unique(d))
    tmp[is.infinite(tmp)] <- NA
    if(!is.null(dim(tmp))){
      score <- mean(c(colMeans(tmp, na.rm = T), rowMeans(tmp, na.rm = T)), na.rm = T)
      
    }else{
      score <- mean(tmp, na.rm = T)
    }
    z <- getzScore(drug, net, length(t), length(d), t, d, iter = iter, s = score)
    
    return(c(score, z, pnorm(z), length(t), mean(degree(net, v = t))))
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


sars_human_ppi <- read.csv("./input/sars-cov2_human_ppi.csv", header = T)
disease <- unique(sars_human_ppi$Human)

human_ppi <- read.csv("./input/human_ppi.csv", header = T)
prot2name_map <- as.data.frame(unique(rbind(cbind(human_ppi$uniprot1,human_ppi$symbol1), cbind(human_ppi$uniprot2, human_ppi$symbol2))))
colnames(prot2name_map) <- c("UniProt.ID", "Gene")

# Add Gene Names to DrugBank
DrugBank <- read.csv("./input/uniprot links.csv", header = T)
db <- merge(as.data.table(DrugBank),as.data.table(prot2name_map), by.x = "UniProt.ID", by.y = "UniProt.ID", all.x = TRUE, all.y = F)

expr <- read.csv("./input/expr.csv", header = T, row.names = 1)
expr <- expr[,which(is.na(str_extract(colnames(expr), "Mock")))] # Exclude Mock
expr <- expr[,which(is.na(str_extract(colnames(expr), "_B_")))] # Exclude _B_
expr <- 2^expr

expr.ch <- expr[,which(!is.na(str_extract(colnames(expr), "_CH_")))]
expr.ya <- expr[,which(!is.na(str_extract(colnames(expr), "_YA_")))]
expr.oa <- expr[,which(!is.na(str_extract(colnames(expr), "_OA_")))]

expr.ch <- rowMeans(expr.ch)
ch <- names(which(expr.ch > mean(expr.ch)))

expr.ya <- rowMeans(expr.ya)
ya <- names(which(expr.ya > mean(expr.ya)))

expr.oa <- rowMeans(expr.oa)
oa <- names(which(expr.oa > mean(expr.oa)))

# Build the PPI network

ppi.ch <- subset(human_ppi, human_ppi$symbol1 %in% ch & human_ppi$symbol2 %in% ch)
ppi.ya <- subset(human_ppi, human_ppi$symbol1 %in% ya & human_ppi$symbol2 %in% ya)
ppi.oa <- subset(human_ppi, human_ppi$symbol1 %in% oa & human_ppi$symbol2 %in% oa)

# Children:
dt.c <- subset(db, db$Gene %in% ch)
drugs.c <- unique(dt.c$DrugBank.ID)

net.ch <- graph_from_data_frame(ppi.ch[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

child.dist <- sapply(drugs.c[1:10], getShortestDistance, net.ch, dt.c, disease, iter = 10)
child.dist <- t(child.dist)
colnames(child.dist) <- c("mean_proximity", "z-score", "p-value", "No.of.Targets", "Mean_degrees")
write.csv(child.dist, file = "./output/Child.csv", quote = F)


# Younger Adault:
dt.y <- subset(db, db$Gene %in% ya)
drugs.y <- unique(dt.y$DrugBank.ID)

net.ya <- graph_from_data_frame(ppi.ya[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

young.dist <- sapply(drugs.y, getShortestDistance, net.ya, dt.y, disease, iter = 100)
young.dist <- t(young.dist)
colnames(young.dist) <- c("mean_proximity", "z-score", "p-value", "No.of.Targets", "Mean_degrees")
write.csv(young.dist, file = "./output/YoungerAdult.csv", quote = F)


# Older Adault:
dt.o <- subset(db, db$Gene %in% oa)
drugs.o <- unique(dt.o$DrugBank.ID)

net.oa <- graph_from_data_frame(ppi.oa[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

old.dist   <- sapply(drugs.o, getShortestDistance, net.oa, dt.o, disease, iter = 100)
old.dist <- t(old.dist)
colnames(old.dist) <- c("mean_proximity", "z-score", "p-value", "No.of.Targets", "Mean_degrees")
write.csv(old.dist, file = "./output/OlderAdult2.csv", quote = F)

