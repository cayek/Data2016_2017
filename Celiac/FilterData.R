################################################################################
library(pacman)
p_load(tidyverse)
p_load_gh("privefl/bigsnpr")
p_load(purrr)
p_load(readr)
p_load("RcppCNPy")
################################################################################
# Read data
load("celiac.RData")

################################################################################
# filter

#na.prop <- colMeans(is.na(celiac.data$genotype))
# max missing data for a snips is 0.1590795
#na.prop <- is.na(celiac.data$genotype) %>% rowMeans()
# max missing data prop for indiv is 0.01580244
# So we keep everyone !!


################################################################################
# Only chrm 6

celiac.data.chrm6 <- list()
celiac.data.chrm6$indiv.covariable <- celiac.data$indiv.covariable
celiac.data.chrm6$chrm6.index <- which(celiac.data$locus.covariable$chromosome == 6)
celiac.data.chrm6$genotype <- celiac.data$genotype[,celiac.data.chrm6$chrm6.index]
celiac.data.chrm6$locus.covariable <- celiac.data$locus.covariable[celiac.data.chrm6$chrm6.index,]

str(celiac.data.chrm6)
object.size(celiac.data.chrm6)  %>% format(units = "GB")
save(celiac.data.chrm6, file = "celiac_chrm6.RData")
# npySave("celiac_chrm6_spnsmatrix.npy", celiac.data.chrm6$genotype) # use a lot of memory ... WHY !!
# load("celiac_chrm6.RData")
write_delim(as.data.frame(celiac.data.chrm6$genotype), "celiac_chrm6_spnsmatrix.txt")
