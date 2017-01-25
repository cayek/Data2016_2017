# We want to extract a snips matrix with only european

# libs
library(tidyverse)

# read indiv informations
indiv <- read_delim("./integrated_call_samples_v3.20130502.ALL.panel",
                    delim = "\t",
                    skip = 1,
                    col_names = FALSE)

names(indiv) <- c("sample", "pop", "super_pop","gender")

unique(indiv %>% select(super_pop))
EUR.index <- which(indiv$super_pop == "EUR")
length(EUR.index)

# convert vcf into geno, il parrait que ca marche avec LEA::vcf2lfmm pour les 1000 genomes                                        ;)
LEA::vcf2geno(
  "./ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf")
# Rmk: LEA::vcf2geno filter not SNPs, means not
geno <- LEA::read.geno("./ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.geno")
# there is no col and name of course....
# il parrait que c'est dans le .vcfsnp ...
snps.info <- read_delim("./ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcfsnp",
                        delim = " ",
                        col_names = FALSE,
                        progress = FALSE)

# keep only europe
# geno.EUR <- matrix(as.raw(geno[EUR.index,]),
#                    nrow = length(EUR.index),
#                    ncol = ncol(geno))
geno.EUR <- geno[EUR.index,]

# Filter
## maf > 0.05 %
maf <- apply(geno.EUR, 2, function(locus) {p <- mean(locus); min(p, 1 - p)})
maf.threshold <- 0.03
mean(maf < maf.threshold)
filtered.index <- (maf >= maf.threshold)
geno.EUR.filtered <- geno.EUR[,filtered.index]

# names con and row
rownames(geno.EUR.filtered) <- indiv$sample[EUR.index]
colnames(geno.EUR.filtered) <- snps.info$X3[filtered.index]

# save
saveRDS(geno.EUR.filtered, file = "European_Chrm22.rds")



