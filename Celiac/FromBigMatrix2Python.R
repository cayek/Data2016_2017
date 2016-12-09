################################################################################
library(pacman)
p_load("bigmemory")
p_load_gh("privefl/bigsnpr")
p_load("purrr")
p_load(readr)
# p_load("RcppCNPy")
################################################################################
celiac <- AttachBigSNP("celiac", backingpath = "./")

celiac.data <- list(genotype = celiac$genotypes[],
  indiv.covariable = celiac$fam,
  locus.covariable = celiac$map)

## info on the size
dim(celiac.data$genotype)
typeof(celiac.data$genotype)
object.size(celiac.data$genotype) %>% format(units = "GB")

## Save the RData sctruct
save(celiac.data, file = "celiac.RData")

## Save genotype into numpy array data
# np_save... can not be installed

write_delim(as.data.frame(celiac.data$genotype), path = "celiac_geno.txt")
