source("./OverRep.R")
source("./getFunctions.R")

data  <- read.table("./data/tpm.txt", header = T, sep = "\t", row.names = 1)
pheno <- read.table("./data/pheno.txt", header = T, sep = "\t")
data  <- data[,which(pheno$Tissue == "bronchial")]
pheno <- pheno[which(pheno$Tissue == "bronchial"),]


b.de   <- read.csv("./data/Bronchial_MockvsVirusDEGs.csv", row.names = 1)
lncRNAs <- read.table("./data/lncRNAs.txt", header = T)[,1]

up <- data[rownames(subset(b.de, logFC>0)),]
covid.up <- getFunctions(up, lncRNAs, "COVID-19_Related_Gene_Sets.txt", fname = "COVID_19/Up_Bronchial_COVID_19")
bp.up <- getFunctions(up, lncRNAs, "GO_Biological_Process_2018.txt", fname = "Gene_Ontology/Up_Bronchial_GO_BP")
path.up <- getFunctions(up, lncRNAs, "KEGG_2019_Human.txt", fname = "KEGG_Pathway/Up_Bronchial_Pathway")


down <- data[rownames(subset(b.de, logFC<0)),]
covid.down <- getFunctions(down, lncRNAs, "COVID-19_Related_Gene_Sets.txt", fname = "COVID_19/Down_Bronchial_COVID_19")
bp.down <- getFunctions(down, lncRNAs, "GO_Biological_Process_2018.txt", fname = "Gene_Ontology/Down_Bronchial_GO_BP")
path.down <- getFunctions(down, lncRNAs, "KEGG_2019_Human.txt", fname = "KEGG_Pathway/Down_Bronchial_Pathway")
