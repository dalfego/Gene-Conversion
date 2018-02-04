#---------------------------------------------------------
# Script Prep
#---------------------------------------------------------

setwd("~/Documents/Resume:CV/Insight/Insight")
ge <- read.csv(file = "qui_example.csv", header = TRUE, sep = ',')
ge <- as.data.frame(ge)

#---------------------------------------------------------
# Selection of Highest and Lowest Expressed Genes
#---------------------------------------------------------

# Sort by Fold value (low to high) and concatenate to gene name
ge_ascend <- ge[order(ge$fold, decreasing = FALSE), c(13,5,14)];
ge_descend <- ge[order(ge$fold, decreasing = TRUE),c(13,5,14)];

# Select 75 lowest and 75 highest expressed genes
low <- ge_ascend[1:75,1:3];
high <- ge_descend[1:75,1:3];
gene_selection <- rbind(high,low)


#---------------------------------------------------------
# Gene Name Conversion to RefSeq for use in pscan
#---------------------------------------------------------

# load library
library("biomaRt")

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
refseq <- getBM(filters = "ensembl_gene_id",
                attributes=c("ensembl_gene_id", "refseq_mrna"), 
                values = gene_selection[1], 
                mart = ensembl)

## Conversion often yields multiple RefSeq per Ensembl Gene or sometimes no conversion (so empty cell)
# Remove rows with empty cell
refseq[refseq==""] <- NA
refseq_final <- refseq[complete.cases(refseq),]

# Only first RefSeq is necessary - remove duplicates based on Ensembl produced by biomaRt
refseq_final <- refseq_final[!duplicated(refseq$ensembl_gene_id),]
refseq_final <- refseq_final[complete.cases(refseq_final),]

