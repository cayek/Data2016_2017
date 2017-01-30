library(tidyverse)

# read getoype
G.file <- "./data/Hgdp_Li_SNP.lfmm"
G <- LEA::read.lfmm(G.file)

## NAs
G[G == 9] <- NA

# read inviv and locus name
snp.name <- read_delim("snp_infos/Hgdp_Lie_name.txt", delim = " ",
                       col_names = FALSE)
indiv.name <- read_delim("ind/Lie_ind_pop.txt", delim = " ", col_names = FALSE)
rownames(G) <- indiv.name$X1
colnames(G) <- snp.name$X1
n <- nrow(G)
L <- ncol(G)


# indiv meta
indiv.meta <- read_delim("ind/Lie_ind_pop.txt", delim = " ", col_names = FALSE)
colnames(indiv.meta) <- c("indiv.name","pop")

# read covariable
indiv.meta <- indiv.meta %>%
  mutate(X_tmp = read_delim("cov/covariable_hgdp_Lie_PC1-tmp.txt", delim = " ", col_names = FALSE)$X1)
indiv.meta <- indiv.meta %>%
  mutate(X_prec = read_delim("cov/covariable_hgdp_Lie_PC1-prec.txt", delim = " ", col_names = FALSE)$X1)
X_tmp <- matrix(indiv.meta$X_tmp, n,1)
rownames(X_tmp) <- rownames(G)
X_prec <- matrix(indiv.meta$X_prec, n,1)
rownames(X_prec) <- rownames(G)

# save a sample
sample.locus <- sample.int(L, size = 50000)
saveRDS(G[,sample.locus], "Hgdp_Li.sample.rds")

# save
saveRDS(X_tmp, "X_tmp.rds")
saveRDS(X_prec, "X_prec.rds")
saveRDS(G, "Hgdp_Li.rds")
