rm(list = ls())

setwd("c:/Users/L1000_NMF_test/")


#=== https://maayanlab.cloud/L1000FWD/download_page
#===  CD_signatures_binary_42809.gmt	
#===  CD signatures (up/down gene sets) in the full space in gmt format.

file1 <- "c:/Users/m092469/Downloads/CD_signatures_binary_42809.gmt"

#require(qusage)
#r1 <- read.gmt(file1)
#name.vec <- names(r1)

require(GSA)
r1 <- GSA.read.gmt(file1)

pert.name.vec <- r1$geneset.names
name.vec <- paste(r1$geneset.names, r1$geneset.descriptions, sep = "@")

sel.idx <- grep("MDAMB231",name.vec)
head(name.vec[sel.idx])

r1$genesets[head(sel.idx,3)]

union.sel.gene.vec <- unique(unlist(r1$genesets[sel.idx]))
union.UpDW.gene.vec <- c(paste0(union.sel.gene.vec,"@up"),
                         paste0(union.sel.gene.vec,"@down"))

N.pertb <- length(unique(pert.name.vec[sel.idx]))
N.gene  <- length(union.UpDW.gene.vec)

X.mtx <- mat.or.vec(N.gene, N.pertb)
rownames(X.mtx) <- union.UpDW.gene.vec
colnames(X.mtx) <- unique(pert.name.vec[sel.idx])

for(k in 1 : length(sel.idx)){
  # k <- 1
  tmp.pertb.idx <- sel.idx[k]
  tmp.gene.vec <- unlist(r1$genesets[tmp.pertb.idx])
  if(r1$geneset.descriptions[tmp.pertb.idx]=="up"){
    tmp.gene.vec <- paste0(tmp.gene.vec,"@up")
  } else {
    tmp.gene.vec <- paste0(tmp.gene.vec,"@down")
  }
  cat(k, name.vec[k], "\n",pert.name.vec[tmp.pertb.idx],"\n")
  X.mtx[tmp.gene.vec,pert.name.vec[tmp.pertb.idx]] <- 1
}

N.times.affected <- apply(X.mtx, 1, sum) 

save(list = c("X.mtx"), file = "MDAMB231_all.Rdata")