library(dada2)
setwd("/work/users/c/a/cabscott/long_coral/SymPortal_framework/Nunez2003/fastqs/")

load("Nunez_OTUtable.ITS2.Rdata")

get.taxa = assignTaxonomy(colnames(seqtab.nochim), refFasta = "/work/users/c/a/cabscott/long_coral/SymportalDADA2_RefDB.fasta", multithread = 12)
save(get.taxa, file = "get.taxa.nunez.Rdata")
