library(tidyverse)

# Simulated data:
DATA_FILE <- 'demo_files/demo_datafile.txt'             # methylation levels file path
PHENO_FILE <- 'demo_files/demo_phenotype.txt'           # phenotype file path
CELL_COMP_FILE <- 'demo_files/demo_cellproportions.txt' # cell composition file path


# methylation lvl
O <- read_delim('demo_files/demo_datafile.txt', delim = "\t")
G <- t(as.matrix(O[,-1]))
colnames(G) <- O$ID

saveRDS(G, "betanormalized_metylationlvl.rds")


# phenotype
phenotype <- read_delim('demo_files/demo_phenotype.txt', delim = "\t", col_names = FALSE)
X <-  as.matrix(phenotype[, -1])
rownames(X) <- phenotype$X1
saveRDS(X, "phenotype.rds")
