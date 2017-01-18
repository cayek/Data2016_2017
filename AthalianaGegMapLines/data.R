# ddl data
# http://bergelson.uchicago.edu/wp-content/uploads/2015/04/call_method_75.tar.gz

################################################################################
## read data
library(data.table)
data.dir <- "./call_method_75/"
ancestral.alleles <- fread(paste0(data.dir,"../ancestral_alleles.txt"),
  data.table = TRUE)
accessions_info <- fread(paste0(data.dir,"/call_method_75_info.tsv"),
  data.table = TRUE)

TAIR8.array_id <- fread(paste0(data.dir,"call_method_75_TAIR8.csv"), nrows = 1, sep = ",", header=FALSE, skip = 0, data.table = TRUE)
TAIR8.snps <- fread(paste0(data.dir,"call_method_75_TAIR8.csv"), sep = ",", header=TRUE, skip = 1, data.table = TRUE)

################################################################################
## keep europe
bb <- c(min(accessions_info$longitude, na.rm = TRUE), min(accessions_info$latitude, na.rm = TRUE), max(accessions_info$longitude, na.rm = TRUE), max(accessions_info$latitude, na.rm = TRUE))

################################################################################
## format data
## we take col0 as reference genotype
call_method_75_TAIR8.europe = list()

# find col0
col0 <- accessions_info[grepl("Col", accessions_info$nativename)]
genome.ref <- TAIR8.snps[,as.character(col0$ecotype_id), with=FALSE]

# order in accessions_info
ecotype_id <- as.character(accessions_info$ecotype_id)

call_method_75_TAIR8.europe$X <- apply(TAIR8.snps[,ecotype_id, with=FALSE], 2, function(c) as.integer(genome.ref == c))
call_method_75_TAIR8.europe$X <- t(call_method_75_TAIR8.europe$X)

call_method_75_TAIR8.europe$coord <- as.matrix(accessions_info[,c("longitude","latitude"), with=FALSE])

# filter european
european.id <- which(call_method_75_TAIR8.europe$coord[,1] > -10 & call_method_75_TAIR8.europe$coord[,1] < 40 & call_method_75_TAIR8.europe$coord[,2] > 35 & call_method_75_TAIR8.europe$coord[,2] < 65)


call_method_75_TAIR8.europe$coord <- call_method_75_TAIR8.europe$coord[european.id,]
colnames(call_method_75_TAIR8.europe$coord) <- c("long", "lat")
rownames(call_method_75_TAIR8.europe$coord) <- accessions_info$nativename[european.id]
call_method_75_TAIR8.europe$locus.coord <- TAIR8.snps[,c("Chromosome", "Positions"), with = FALSE]
call_method_75_TAIR8.europe$X <- call_method_75_TAIR8.europe$X[european.id,]
aux <- apply(call_method_75_TAIR8.europe$locus.coord, 1, function(x) paste('chr:', x[1], 'pos:', x[2]))
colnames(call_method_75_TAIR8.europe$X) <- aux
rownames(call_method_75_TAIR8.europe$X) <- accessions_info$nativename[european.id]

save(call_method_75_TAIR8.europe, file = paste0(data.dir,"/call_method_75_TAIR8.RData"))

# in .geno
library(LEA)
write.geno(call_method_75_TAIR8.europe$X, paste0(data.dir,"/call_method_75_TAIR8.geno"))

################################################################################
## TAIR9 data
rm(list = ls()[grepl(".*TAIR8.*", ls())])

TAIR9.array_id <- fread(paste0(data.dir,"call_method_75_TAIR9.csv"), nrows = 1, sep = ",", header=FALSE, skip = 0, data.table = TRUE)
TAIR9.snps <- fread(paste0(data.dir,"call_method_75_TAIR9.csv"), sep = ",", header=TRUE, skip = 1, data.table = TRUE)

call_method_75_TAIR9.europe = list()

# find col0
genome.ref <- TAIR9.snps[,as.character(col0$ecotype_id), with=FALSE]

# order in accessions_info

call_method_75_TAIR9.europe$X <- apply(TAIR9.snps[,ecotype_id, with=FALSE], 2, function(c) as.integer(genome.ref == c))
call_method_75_TAIR9.europe$X <- t(call_method_75_TAIR9.europe$X)

call_method_75_TAIR9.europe$coord <- as.matrix(accessions_info[,c("longitude","latitude"), with=FALSE])

# filter european
european.id <- which(call_method_75_TAIR9.europe$coord[,1] > -10 & call_method_75_TAIR9.europe$coord[,1] < 40 & call_method_75_TAIR9.europe$coord[,2] > 35 & call_method_75_TAIR9.europe$coord[,2] < 65)


call_method_75_TAIR9.europe$coord <- call_method_75_TAIR9.europe$coord[european.id,]
colnames(call_method_75_TAIR9.europe$coord) <- c("long", "lat")
rownames(call_method_75_TAIR9.europe$coord) <- accessions_info$nativename[european.id]
call_method_75_TAIR9.europe$locus.coord <- TAIR9.snps[,c("Chromosome", "Positions"), with = FALSE]
call_method_75_TAIR9.europe$X <- call_method_75_TAIR9.europe$X[european.id,]
aux <- apply(call_method_75_TAIR9.europe$locus.coord, 1, function(x) paste('chr:', x[1], 'pos:', x[2]))
colnames(call_method_75_TAIR9.europe$X) <- aux
rownames(call_method_75_TAIR9.europe$X) <- accessions_info$nativename[european.id]

save(call_method_75_TAIR9.europe, file = paste0(data.dir,"/call_method_75_TAIR9.RData"))
