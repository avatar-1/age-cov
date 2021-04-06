getFunctions <- function(expr, lncRNAs, pop, fdr = 0.05, fname){
  c <- cor(t(expr), method = "pearson")
  diag(c) <- 0
  #z <-corrplot(c, type = "full",  order = "hclust", method = "color", tl.col = "black", tl.cex = 0.45)# tl.srt = 45)
  out <- c()
  out2 <- c()
  lincs <- intersect(rownames(expr),lncRNAs)
  for(i in lincs){
    print(i)
    tmp <- try(OverRep(names(which(c[i,]>0.8)), pops = pop))
    if(inherits(tmp, "try-error"))
    {
      next
    }
    tmp <- subset(tmp, FDR < fdr)
    tmp <- subset(tmp, no.overlap > 1)
    
    out <- c(out, tmp[,"Gene.set"])
    if(length(tmp)>0){
      if(nrow(tmp) >0){
        for(j in 1:nrow(tmp)){
          out2 <- rbind(out2, c(i, tmp[j,1], tmp[j,"FDR"], trimws(paste(tmp[j, 11:ncol(tmp)], collapse = " "))))}
      }}}
  res <- list()
  res[[1]] <- unique(out)
  colnames(out2) <- c("lncRNA", "Gene Set", "FDR", "Overlaping Genes")
  res[[2]] <- out2
  #write.csv(res[[1]], file = paste("./results/",fname,"_genesets.csv", sep = ""), row.names = F, col.names = F)
  write.csv(res[[2]], file = paste("./results/",fname,"_lncRNAs.csv", sep = ""), row.names = F, col.names = F)
  
  return(res)
}
