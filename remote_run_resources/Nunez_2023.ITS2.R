library(dada2)
setwd("/work/users/c/a/cabscott/long_coral/SymPortal_framework/Nunez2003/fastqs/")

path <- "."
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


library(ShortRead)

counts <- data.frame(
  Forward = fnFs,
  Reverse = fnRs,
  nF = sapply(fnFs, function(f) length(readFastq(f))),
  nR = sapply(fnRs, function(f) length(readFastq(f)))
)

# Find mismatched pairs
mismatched <- counts$nF != counts$nR
cat("Mismatched pairs:\n")
print(counts[mismatched, ])

# Keep only matched pairs
fnFs <- fnFs[!mismatched]
fnRs <- fnRs[!mismatched]
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# filtFs <- na.omit(filtFs[!mismatched])
# filtRs <- na.omit(filtRs[!mismatched])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 20, trimRight = 50,
                     maxN=0, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=1) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=1)
errR <- learnErrors(filtRs, multithread=1)

dadaFs <- dada(filtFs, err=errF, multithread=1)
dadaRs <- dada(filtRs, err=errR, multithread=1)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=1, verbose=TRUE)
save(seqtab.nochim, file = "Nunez_OTUtable.ITS2.Rdata")

load("Nunez_OTUtable.ITS2.Rdata")
library(dada2)

get.taxa = assignTaxonomy(colnames(seqtab.nochim), refFasta = "/work/users/c/a/cabscott/long_coral/SymportalDADA2_RefDB.fasta", multithread = 1)
save(get.taxa, file = "get.taxa.nunez.Rdata")