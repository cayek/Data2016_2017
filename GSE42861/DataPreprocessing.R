X <- readRDS("X.rds")
G <- readRDS("betanormalized_metylationlvl.rds")


# filter maf !
maf <- apply(G, 2, function(l){p <- mean(l);min(p, 1 - p)})
out.index <- which(maf <= 0.2)

G.filtered <- G[,-out.index]
dim(G.filtered)

saveRDS(G.filtered, "betanormalized_metylationlvl.filtered.rds")
