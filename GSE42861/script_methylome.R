## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("Biobase")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")


require(Biobase)
require(GEOquery)

## get le jeu de données dans le format biobase
obj861 <- getGEO("GSE42861",GSEMatrix = T)

## extrait les phenotypes (factors)
disease.state <- pData(phenoData(obj861[[1]]))[,11]

## extrait les covariables (subject, age, gender, smocking.status)
## age est converti en numeric

subject <- pData(phenoData(obj861[[1]]))[,12]

age.f <- pData(phenoData(obj861[[1]]))[,13]
write.table(file = "age.txt", as.character(age.f))
age <- as.numeric(read.table(file = "age.txt")[,1])


gender <- pData(phenoData(obj861[[1]]))[,14]

smocking.status <- pData(phenoData(obj861[[1]]))[,15]

## download la matrice d'expression. Attention elle est transposée (individus en colonnes)
expmat861 <- exprs(obj861[[1]])



