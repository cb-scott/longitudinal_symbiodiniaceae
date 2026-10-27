library(tximport)
library(DESeq2)
library(fgsea)
library(dplyr)
library(tidyr)
library(rtracklayer)

#this script is missing the download of the brite files
load("data/combinedMeta.Rdata")

#### ----- BRITE hierarchy (download once) ----- ####
brite.rds <- "data/kegg_brite_ko00001.rds"
if (file.exists(brite.rds)) {
  brite <- readRDS(brite.rds)
} else {
  js <- fromJSON("https://rest.kegg.jp/get/br:ko00001/json", simplifyVector = FALSE)
  
  # Walk the tree; each leaf is a KO with its A/B/C path
  walk <- function(node, path = character(0)) {
    if (is.null(node$children)) return(list(c(path, node$name)))
    unlist(lapply(node$children, walk, path = c(path, node$name)), recursive = FALSE)
  }
  leaves <- walk(js)
  leaves <- leaves[lengths(leaves) == 5]  # root, A, B, C, KO
  
  brite <- as.data.frame(do.call(rbind, leaves))[, 2:5]
  colnames(brite) <- c("A", "B", "C", "D")
  brite <- brite %>%
    mutate(KO = sub("^(K[0-9]{5}).*", "\\1", D),
           is.pathway = !grepl("^091[89]", A),  # 0918x/0919x = BRITE families, not pathways
           across(c(A, B, C), ~ sub("^[0-9]+ ", "", .x)),
           C = sub(" \\[.*\\]$", "", C)) %>%
    filter(grepl("^K[0-9]{5}$", KO)) %>%
    distinct(KO, A, B, C, is.pathway)
  saveRDS(brite, brite.rds)
}



comp     <- "rhizoc"
out.dir  <- "results/fgsea_rhizoc"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
coef     <- "treatment_MOCK_vs_CHAMBER"   # positive = higher in MOCK
#loosing too much power w/full model; use reduced without microbiome
padj.cut <- 0.1

#### ----- Salmon files ----- ####
lol.samples <- list.files("data/FINAL_salmon/", pattern = ".filt.gz.rhizoc.quant")
nice.nam    <- gsub(".filt.gz.rhizoc.quant", "", lol.samples)
lol.files   <- file.path("data", "FINAL_salmon", lol.samples, "quant.sf")
names(lol.files) <- nice.nam

#### ----- tx2gene from GFF ----- ####
gff.df <- read.delim("data/annotations/S1W6_Primary_Evidences_Final_Annotation_5_18_26.gff3", header = FALSE,
                     comment.char = "#", quote = "")
table(gff.df$V3)  # feature types

tx.rows <- gff.df[gff.df$V3 %in% c("mRNA", "transcript"), ]
tx2gene <- data.frame(
  TXNAME = sub(".*ID=([^;]+).*",     "\\1", tx.rows$V9),
  GENEID = sub(".*Parent=([^;]+).*", "\\1", tx.rows$V9)
)
head(tx2gene)


all.conts <- read.table(lol.files[1], header = TRUE)
mean(all.conts$Name %in% tx2gene$TXNAME)  # should be ~1
tx2gene <- tx2gene[tx2gene$TXNAME %in% all.conts$Name, ]

#### ----- Gene-level import ----- ####
txi <- tximport(lol.files, type = "salmon", tx2gene = tx2gene)
dim(txi$counts)  # genes x samples

# colData must match txi column order
stopifnot(all(colnames(txi$counts) %in% rownames(meta.gen)))
dds <- DESeqDataSetFromTximport(txi, colData = meta.gen[colnames(txi$counts), ],
                                design = ~ lol_genMDS1 + lol_genMDS2 + rhiz_genMDS1 + rhiz_genMDS2 + treatment)

dds$isrep <- factor(gsub("_TR", "", colnames(dds)))
dds <- collapseReplicates(dds, dds$isrep)

#### ----- Blank and prevalence filters ----- ####
blank.q <- quantile(counts(dds)[, "EXTRACT_BLANK"], probs = seq(0, 1, .01))
dds <- dds[counts(dds)[, "EXTRACT_BLANK"] <= blank.q[95], ]
dds <- dds[, colnames(dds) != "EXTRACT_BLANK"]

smallestGroupSize <- 4
dds  <- estimateSizeFactors(dds, type = "poscounts")
dds  <- dds[rowSums(counts(dds, normalized = TRUE) > 0) >= smallestGroupSize, ]
dds  <- estimateSizeFactors(dds, type = "poscounts")
dim(dds)

dds <- DESeq(dds)
res <- results(dds)
data.frame(res) %>% arrange(padj)

#res_no_filter <- results(dds, independentFiltering = FALSE)

resLFC <- lfcShrink(dds, coef="treatment_MOCK_vs_CHAMBER", type = "apeglm") #almost no effect of treatment after we account for host + patho genotype
data.frame(resLFC) %>% arrange(padj)
plotMA(resLFC) #not sure whay this looks bad, many genes look fine

data.frame(res) %>%
  ggplot(aes(x=log2FoldChange, y = -1*log10(padj))) + geom_point() + 
 geom_hline( yintercept = -1*log10(.05)) + geom_vline(xintercept = c(-2, 2)) +
  theme_classic()

res.df = data.frame(res) %>% filter(is.na(padj) == F) %>%
  mutate(signif = ifelse(padj <= 0.05, 1, 0),
         updown = case_when(log2FoldChange > 2~"up",
                            log2FoldChange < -2~"down",
                            .default = NA),
         tocol = case_when(signif == 1 & updown == "up" ~ plotcols[1],
                           signif == 1 & updown == "down" ~ plotcols[2],
                           signif == 0 | is.na(updown) == T ~ "gray80"))
res.df %>%  #mutate(log2FoldChange = ifelse(log2FoldChange < -10, -10, ifelse(log2FoldChange > 10, 10, log2FoldChange))) %>%
  ggplot(aes(x=log2FoldChange, y = -1*log10(padj))) +
  geom_point(aes(col = tocol), alpha = .5, size = 1) + 
  geom_vline(xintercept = c(-2, 2), col = "gray40", lty = "dashed") + geom_hline( yintercept = -1*log10(.05), col = "gray40") +
  scale_color_manual(values = c(plotcols[c(1)], "gray80")) +
  theme_classic() + theme(legend.position =  "none", text = element_text(size = 8)) + 
  ylab(expression("log"[10]*"(p"["adj"]*")"))+ labs(x = "log2FoldChange") +
  annotate(
    "text", x = -Inf, y = Inf,
    label = "← Mock",
    hjust = -.1, vjust = 2, size = 8/2.84,
    col = plotcols[2]
  ) +
  annotate(
    "text", x = Inf, y = Inf,
    label = "Chamber →",
    hjust = 1, vjust = 2, size = 8/2.84,
    col = plotcols[1]
  ) +
  coord_cartesian(clip = "off")# +
ggsave("results/DEG_Rhizoc.svg", width = 2.1, height = 2, units = "in")


#### ----- try a ko level DESeq analysis -------  #######

#FOR LOLIUM, DON'T NECESSARILY WANT THE KO-LEVEL ANALYSIS; only to look for enrichment
#Keep for consistency here, ultimately will want the KO-LEVEL analysis for the metatranscriptome

#### ----- gene2ko ----- ####
read_emapper <- function(file) {
  hdr <- grep("^#query", readLines(file, n = 200))
  ann <- read.delim(file, skip = hdr - 1, header = TRUE, quote = "",
                    comment.char = "", na.strings = "-", fill = TRUE,
                    check.names = FALSE)
  colnames(ann)[1] <- "query"
  ann[!grepl("^##", ann$query), ]
}
rhiz.ann <- read_emapper("data/annotations/eggnog_rhizoc.emapper.annotations")  # your path
dim(rhiz.ann)
head(rhiz.ann[, c("query", "Preferred_name", "KEGG_ko", "Description")])

rhiz.ann$transcript <- rhiz.ann$query

gene2ko <- rhiz.ann %>%
  filter(!is.na(KEGG_ko)) %>%
  select(transcript, KEGG_ko) %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KO = sub("^ko:", "", KEGG_ko)) %>%
  inner_join(tx2gene, by = c("transcript" = "TXNAME")) %>%
  distinct(gene = GENEID, KO)

table(table(gene2ko$gene))                 # KOs per gene
mean(rownames(dds) %in% gene2ko$gene)      # KO coverage of tested genes


#### ----- Gene-to-KO mapping for tested genes ----- ####
g2k <- gene2ko %>% filter(gene %in% rownames(dds))

# Genes with several KOs: summing them into each KO counts their reads twice
multi <- g2k %>% dplyr::count(gene) %>% filter(n > 1)
c(genes.with.KO = n_distinct(g2k$gene), multi.KO.genes = nrow(multi))

g2k1 <- g2k %>% filter(!gene %in% multi$gene)  # drop multi-KO genes

#### ----- Sum gene counts to KOs ----- ####
ko.cts <- rowsum(counts(dds)[g2k1$gene, ], g2k1$KO)
dim(ko.cts)

# Start fresh: don't carry over gene-level size factors
cd <- colData(dds)
cd$sizeFactor <- NULL

dds.ko <- DESeqDataSetFromMatrix(round(ko.cts), colData = cd, design = design(dds))
dds.ko <- dds.ko[rowSums(counts(dds.ko) > 0) >= smallestGroupSize, ]
dds.ko <- DESeq(dds.ko)

res.ko <- results(dds.ko, name = coef)
summary(res.ko)
resLFC.ko <- lfcShrink(dds.ko, coef = coef, type = "apeglm")

#### ----- Tidy KO results ----- ####
#where is ko.names?
ko.de <- as.data.frame(res.ko) %>%
  mutate(KO = rownames(res.ko),
         lfc.shrunk = resLFC.ko[KO, "log2FoldChange"],
         name  = ko.names[KO],
         label = paste(KO, sub(";.*", "", name))) %>%
  left_join(ko.cat, by = "KO") %>%
  mutate(A = coalesce(A, "Unclassified"), B = coalesce(B, "Unclassified")) %>%
  arrange(padj)
write.csv(ko.de, file.path(out.dir, paste0(comp, "_KO_DESeq2.csv")), row.names = FALSE)

#### ----- Volcano colored by category ----- ####
ggplot(filter(ko.de, !is.na(padj)), aes(lfc.shrunk, -log10(padj), col = A)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(padj.cut), lty = 2) +
  labs(x = "log2 fold change (positive = higher in MOCK)", y = "-log10(padj)",
       col = "KEGG level A") +
  theme_classic()

#### ----- Category enrichment on KO stats ----- ####
brite.c <- brite %>%
  filter(!grepl("Human Diseases|Drug Development", A)) %>%
  distinct(C, KO)
brite.sets <- split(brite.c$KO, brite.c$C)

#this is grabbing the heirarchichal category for each C-level BRITE term
#also is just grabbing one representative function per group
c2B <- brite %>% dplyr::distinct(C, B, A) %>% group_by(C) %>% dplyr::slice(1) %>% ungroup()

ok <- !is.na(res.ko$stat) & !is.na(res.ko$pvalue)
ko.ranks <- setNames(res.ko$stat[ok], rownames(res.ko)[ok])

set.seed(1)
fg.cat     <- fgsea(brite.sets, ko.ranks, minSize = 10, maxSize = 1000)


cat.res <- as.data.frame(fg.cat) %>%
  left_join(c2B, by = c("pathway" = "C")) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
  arrange(padj)
#write.csv(cat.res, file.path(out.dir, paste0(comp, "_category_fgsea.csv")), row.names = FALSE)
cat.res %>% filter(padj < padj.cut) %>% select(pathway, B, size, NES, padj)

#--------CORRECT the ANALYSIS AT THE GENE LEVEL ------
# Let's do an fgsea of genes to Brite Level C hierarchy terms. For each gene, we need it's BRITE C hierarchy, and then 
# we can use C2B to color the volcano plot

#only use genes with one ko term -> the problem is that some KOs have multiple brite annotations. Maybe that's okay?
red.brite.c = brite.c %>% filter(C %in% c2B$C) %>%
  left_join(c2B) %>%
  filter(A %in% c("Cellular Processes", "Environmental Information Processing", 
                  "Genetic Information Processing", "Metabolism", "Organismal Systems")) 

#maybe just let it be a many to many join? let's see what happens with the fgsea
gene2B <- g2k1 %>% left_join(red.brite.c)
uqcs = gene2B %>% dplyr::select(gene, C) %>% distinct()
gene.sets.c <- split(uqcs$gene, uqcs$C)

ok <- !is.na(res$stat) & !is.na(res$pvalue)
gene.ranks <- setNames(res$stat[ok], rownames(res)[ok])

#need lists of categories for fgsea
fg.gene.c    <- fgsea(gene.sets.c, gene.ranks, minSize = 10, maxSize = 1000)
data.frame(fg.gene.c) %>% arrange(padj)

###B level
uqbs = gene2B  %>% dplyr::select(gene, B) %>% distinct()
gene.sets.B <- split(uqbs$gene, uqbs$B)

ok <- !is.na(res$stat) & !is.na(res$pvalue)
gene.ranks <- setNames(res$stat[ok], rownames(res)[ok])

#need lists of categories for fgsea
fg.gene.B    <- fgsea(gene.sets.B, gene.ranks, minSize = 10, maxSize = 1000)
data.frame(fg.gene.B) %>% arrange(padj)

sig.cs = data.frame(fg.gene.c) %>% arrange(padj) %>% left_join(c2B, by = c("pathway" = "C")) %>%
  mutate(ppi = ifelse(pathway %in% c("Protein processing in endoplasmic reticulum", "Photosynthesis", "Betalain biosynthesis", "IgSF CAM signaling"), pathway, ""))

sig.cs %>%
  ggplot(aes(x=NES, y=-1*log10(padj))) + geom_point(aes(color = A)) +
  geom_hline(yintercept = -1*log10(.05)) + geom_text(aes(label = ppi), size = 3) + 
  theme_classic() 



highcs = sig.cs %>% filter(NES > 0, padj <= 0.05) %>% arrange(padj) %>% head(10)
lowcs = sig.cs %>% filter(NES < 0, padj <= 0.05) %>% arrange(padj) %>% head(10)
#NOTHING TO PLOT HERE
