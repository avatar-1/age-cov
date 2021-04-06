OverRep<- function(goi.matrix, pops, URL = FALSE, pop.genes = 45956, k){
#make pop file vector to iterate through
# all pops -> "c1.all.gmt", "c2.all.gmt", "c2.cgp.gmt", "c2.cp.all.gmt","c2.cp.biocarta.gmt", "c2.cp.kegg.gmt", "c2.cp.reactome.gmt", "c3.all.gmt", "c3.mir.gmt", "c3.tft.gmt", "c4.all.gmt", "c4.cgn.gmt", "c4.cm.gmt", "c5.all.gmt", "c5.bp.gmt", "c5.cc.gmt", "c5.mf.gmt", "c6.all.gmt", "c7.all.gmt", "h.all.gmt"
#pops = c( "c1.all.gmt", "c2.cgp.gmt", "c2.cp.biocarta.gmt", "c2.cp.kegg.gmt", "c2.cp.reactome.gmt", "c4.cgn.gmt", "c4.cm.gmt", "c5.bp.gmt", "c5.cc.gmt", "c5.mf.gmt", "c6.all.gmt", "c7.all.gmt", "h.all.gmt", "coag.gmt")
#pops = "c2.cp.kegg.gmt"#c("c5.bp.gmt", "c5.cc.gmt", "c5.mf.gmt")
l.pops = length(pops)

#select the genes of interest (differential in our case)
#goi.file = "LFC2p001bon.gmt"
#goi.filepath = paste("./OverRep/GOI/", goi.file, sep = "")

#cal no genes of interest = n
#read in genes of interest gmt (could be any sort of delineated file with gene names only)
#goi.matrix = read.table(goi.filepath, sep = '\t', header = FALSE, stringsAsFactors = FALSE, fill = TRUE, colClasses = "character")

#sum of no empty cells = num of genes drawn BUT note this will depend on how the gene of interest file is set out!!!
goi.genes =  sum(goi.matrix !='')

for (h in 1:l.pops) {

#select the population
pop.select = pops[h]
pop.filepath = paste("./OverRep/Pop/", pop.select, sep = "")

#calc no. genes in pop of genesets = N
#read in gmt file of genesets
num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
pop.matrix = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
# if(!URL){
#   cat(dim(pop.matrix),"\t")
#   pop.matrix = cbind(pop.matrix[,1],replicate(nrow(pop.matrix), "URL"),pop.matrix[,2:num.col])
#   cat(dim(pop.matrix),"\n")
# }
n.gs = nrow(pop.matrix)
max.gs.size = max(rowSums(pop.matrix!=""))

#calc non empty cells, minus first two columns (description) = no. of genes
popvector = unique(c(as.matrix(pop.matrix)))
pop.genes = length(popvector) - nrow(pop.matrix)*2 - 1 #45956

#Calc num of goi in the gene population
drawn.genes = length(goi.matrix[match(popvector, goi.matrix,nomatch = 0)])


#create the dataframe for storing results
gene.overrep = data.frame(matrix(c(12), nrow = n.gs, ncol = 10, dimnames = list(pop.matrix[,1], c("Gene set","URL","universe size","goi", "goi in universe", "pathway size", "no.overlap", "k over K", "pvalue", "FDR"))), stringsAsFactors = FALSE)
gene.overrep[,1:2] = pop.matrix[,1:2]

#create data frame for storing gene overlaps
gene.lap = data.frame(matrix(c(""), nrow = n.gs, ncol = drawn.genes), stringsAsFactors = FALSE)

#loop through every gs in the pop file
for (x in 1:n.gs) {
  #calc no. genes in Gene set of interest = K
  #calc number of cells - 2 = no. of genes in gs x
  gs.genes = sum(pop.matrix[x,] !='', na.rm = TRUE) - 2
  #calc number of overlapping genes
  match.genes = goi.matrix[match(pop.matrix[x,],goi.matrix,nomatch = 0)]
  overlap.genes = length(match.genes)
  #fill gene overlap matrix
  if(overlap.genes>0) {
    gene.lap[x,] = c(match.genes, rep("",times = drawn.genes - overlap.genes))
  }
  #fill the over rep matrix
  #universe size
  gene.overrep[x,3] = pop.genes
  #goi genes from de
  gene.overrep[x,4] = goi.genes
  #goi genes in universe
  gene.overrep[x,5] = drawn.genes
 # pathway size
  gene.overrep[x,6] = gs.genes
  #no.overlap
  gene.overrep[x,7] = overlap.genes
  #k/K
  gene.overrep[x,8] = overlap.genes/gs.genes
  #pvalue
  gene.overrep[x,9] = phyper(overlap.genes, gs.genes, pop.genes-gs.genes, drawn.genes, lower.tail = FALSE) #He wrote : overlap.genes-1
  #insert FDR val
  #n.gs = length(which(gene.overrep[x,9]> 0))
  gene.overrep[x,10] = p.adjust(gene.overrep[x,9], method = "bonferroni", n = n.gs)
}

gene.overrep = cbind(gene.overrep,gene.lap)
colnames(gene.overrep)[11] = "Overlapping Genes"

gene.overrep = gene.overrep[order(gene.overrep$FDR),]

#write results file
output.or = paste("./OverRep/Output/",pop.select, sep = "")
output.or = paste(output.or,".csv", sep = "")
write.table(gene.overrep, output.or, sep = ",",row.names = FALSE, col.names = colnames(gene.overrep))

}
return(gene.overrep)
}