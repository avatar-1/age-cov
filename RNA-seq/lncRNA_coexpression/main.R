

# Load data ---------------------------------------------------------------

data  <- read.table("./data/tpm.txt", header = T, sep = "\t", row.names = 1)
pheno <- read.table("./data/pheno.txt", header = T, sep = "\t")
data  <- data[,which(pheno$Tissue == "nasal")]
pheno <- pheno[which(pheno$Tissue == "nasal"),]


de   <- read.csv("./data/All_MockvsVirusDEGs.csv", row.names = 1)
c.de <- read.csv("./data/Child_MockvsVirusDEGs.csv", row.names = 1)
y.de <- read.csv("./data/Young_MockvsVirusDEGs.csv", row.names = 1)
o.de <- read.csv("./data/Old_MockvsVirusDEGs.csv", row.names = 1)

lncRNAs <- read.table("./data/lncRNAs.txt", header = T)[,1]


# Data - age groups -------------------------------------------------------

data.c  <- data[rownames(c.de),which(pheno$Age == "Child")]
pheno.c <- pheno[which(pheno$Age == "Child"),c("SampleName", "Class")]

data.y  <- data[rownames(y.de),which(pheno$Age == "Young")]
pheno.y <- pheno[which(pheno$Age == "Young"),c("SampleName", "Class")]



# Correlation analysis ----------------------------------------------------


source("./OverRep.R")

# Enrichment by up regulated lncRNAs --------------------------------------

up.c.only <- data.c[setdiff(rownames(subset(c.de, logFC>0)), rownames(subset(y.de, logFC>0))),]
bp.up.c.only <- getFunctions(up.c.only, lncRNAs, "GO_Biological_Process_2018.txt", fname = "Gene_Ontology/Up_Child_only_GO_BP")
path.up.c.only <- getFunctions(up.c.only, lncRNAs, "KEGG_2019_Human.txt", fname = "KEGG_Pathway/Up_Child_only_Pathway")


up.y.only <- data.y[setdiff(rownames(subset(y.de, logFC>0)), rownames(subset(c.de, logFC>0))),]
bp.up.y.only <- getFunctions(up.y.only, lncRNAs, "GO_Biological_Process_2018.txt", fname = "Gene_Ontology/Up_Young_only_GO_BP")
path.up.y.only <- getFunctions(up.y.only, lncRNAs, "KEGG_2019_Human.txt", fname = "KEGG_Pathway/Up_Young_only_Pathway")


up.y.c.both <- data.c[intersect(rownames(subset(c.de, logFC>0)), rownames(subset(y.de, logFC>0))),]
bp.up.y.c.both <- getFunctions(up.y.c.both, lncRNAs, "GO_Biological_Process_2018.txt", fname = "Gene_Ontology/Up_C&Y_GO_BP")
path.up.y.c.both <- getFunctions(up.y.c.both, lncRNAs, "KEGG_2019_Human.txt", fname = "KEGG_Pathway/Up_C&Y_pathway")


# Enrichment by down regulated lncRNAs --------------------------------------

down.c.only <- data.c[setdiff(rownames(subset(c.de, logFC<0)), rownames(subset(y.de, logFC<0))),]
bp.down.c.only <- getFunctions(down.c.only, lncRNAs, "GO_Biological_Process_2018.txt", fname = "Gene_Ontology/Down_Child_only_GO_BP")
path.down.c.only <- getFunctions(down.c.only, lncRNAs, "KEGG_2019_Human.txt", fname = "KEGG_Pathway/Down_Child_only_Pathway")


down.y.only <- data.y[setdiff(rownames(subset(y.de, logFC<0)), rownames(subset(c.de, logFC<0))),]
bp.down.y.only <- getFunctions(down.y.only, lncRNAs, "GO_Biological_Process_2018.txt", fname = "Gene_Ontology/Down_Young_only_GO_BP")
path.down.y.only <- getFunctions(down.y.only, lncRNAs, "KEGG_2019_Human.txt", fname = "KEGG_Pathway/Down_Young_only_Pathway")


down.y.c.both <- data.c[intersect(rownames(subset(c.de, logFC<0)), rownames(subset(y.de, logFC<0))),]
bp.down.y.c.both <- getFunctions(down.y.c.both, lncRNAs, "GO_Biological_Process_2018.txt", fname = "Gene_Ontology/Down_C&Y_only_GO_BP")
path.down.y.c.both <- getFunctions(down.y.c.both, lncRNAs, "KEGG_2019_Human.txt", fname = "KEGG_Pathway/Down_C&Y_only_Pathway")


# targets = c("ALB", "EIF2AK2", "EIF4E", "IMPDH1", "IMPDH2")
# res <- OverRep(targets, pops = "GO_Biological_Process_2018.txt")


