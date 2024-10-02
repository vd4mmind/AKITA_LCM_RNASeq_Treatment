## Differential expression analysis
## Comparison: Normal_Dapagliflozin_KO_vs_Normal_Vehicle_WT

Library setup
Get personal library location:

R_LIBS_USER <- Sys.getenv("R_LIBS_USER")

dir.create(R_LIBS_USER, recursive = TRUE)
Install required CRAN packages:

if (!require("DT")) install.packages("DT", lib = R_LIBS_USER)

if (!require("ggplot2")) install.packages("ggplot2", lib = R_LIBS_USER)

if (!require("gtools")) install.packages("gtools", lib = R_LIBS_USER)

if (!require("pheatmap")) install.packages("pheatmap", lib = R_LIBS_USER)
Load required CRAN packages:

library(DT, lib.loc = R_LIBS_USER)

library(ggplot2, lib.loc = R_LIBS_USER)

library(gtools, lib.loc = R_LIBS_USER)

library(pheatmap, lib.loc = R_LIBS_USER)
Install required BIOC packages:

if (!require("DESeq2")) BiocManager::install("DESeq2", lib = R_LIBS_USER)

if (!require("EnsDb.Mmusculus.v79")) BiocManager::install("EnsDb.Mmusculus.v79", lib = R_LIBS_USER)

if (!require("GO.db")) BiocManager::install("GO.db", lib = R_LIBS_USER)

if (!require("limma")) BiocManager::install("limma", lib = R_LIBS_USER)

if (!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db", lib = R_LIBS_USER)

if (!require("scater")) BiocManager::install("scater", lib = R_LIBS_USER)
Load required BIOC packages:

library(DESeq2, lib.loc = R_LIBS_USER)

library(EnsDb.Mmusculus.v79, lib.loc = R_LIBS_USER)

library(GO.db, lib.loc = R_LIBS_USER)

library(limma, lib.loc = R_LIBS_USER)

library(org.Mm.eg.db, lib.loc = R_LIBS_USER)

library(scater, lib.loc = R_LIBS_USER)
Set theme and colour palettes:

theme_set(theme_bw())

options(ggplot2.continuous.colour = "viridis")
Data setup
Read counts matrix:

counts <- read.delim("/sbgenomics/project-files/Group_Comparisons/Data_aggregation_outputs/Group_Comparisons.gene.numreads", row.names = "Gene_ID")
Read FPKM matrix:

fpkm <- read.delim("/sbgenomics/project-files/Group_Comparisons/Data_aggregation_outputs/Group_Comparisons.gene.fpkm", row.names = "Gene_ID")
Read TPM matrix:

tpm <- read.delim("/sbgenomics/project-files/Group_Comparisons/Data_aggregation_outputs/Group_Comparisons.gene.tpm", row.names = "Gene_ID")
Ensure consistent gene order:

ids <- Reduce(intersect, list(rownames(counts), rownames(fpkm), rownames(tpm)))

counts <- counts[ids, ]

fpkm <- fpkm[ids, ]

tpm <- tpm[ids, ]
Read manifest table:

manifest <- read.csv("/sbgenomics/project-files/Manifest.csv")
Create sample data from manifest table:

samples <- manifest

samples$file_name <- NULL

samples$paired_end <- NULL

samples <- unique(samples)

rownames(samples) <- samples$sample_id

samples$sample_id <- NULL
Ensure consistent sample order:

ids <- Reduce(intersect, list(rownames(samples), colnames(counts), colnames(fpkm), colnames(tpm)))

ids <- gtools::mixedsort(ids)

samples <- samples[ids, ]

counts <- counts[, ids]

fpkm <- fpkm[, ids]

tpm <- tpm[, ids]
Remove outlier samples:

keep <- setdiff(rownames(samples), c("LCM10", "LCM44", "LCM87", "LCM58", "LCM85", "LCM43"))

samples <- samples[keep, ]

counts <- counts[, keep]

fpkm <- fpkm[, keep]

tpm <- tpm[, keep]
Make group and disease levels syntactically valid names:

samples$group <- sub("Non-Diabetic", "Normal", samples$group)

samples$disease <- sub("Non-Diabetic", "Normal", samples$disease)
Create DDS object:

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = samples,
  design = ~ group
)
Add FPKM and TPM assays:

assay(dds, "fpkm") <- fpkm

assay(dds, "tpm") <- tpm
Filter DDS object to samples in test and reference groups:

ind <- dds$group %in% c(params$test, params$reference)

dds <- dds[, ind]
Drop unused levels from all factors:

colData(dds) <- droplevels(colData(dds))
Retrieve gene symbols:

rowData(dds)$GENEID <- rownames(dds)

rowData(dds)$SYMBOL <- mapIds(EnsDb.Mmusculus.v79, keys = rownames(dds), column = "SYMBOL", keytype = "GENEID")

rowData(dds)$ENTREZID <- mapIds(EnsDb.Mmusculus.v79, rowData(dds)$GENEID, "ENTREZID", "GENEID")

rownames(dds) <- scater::uniquifyFeatureNames(rowData(dds)$GENEID, rowData(dds)$SYMBOL)
Run VSD transformation:

vsd <- vst(dds, blind = TRUE)
Data exploration
Define factors and covariates of interest:

FACTORS <- c(
  "group"     = "Group",
  "disease"   = "Disease",
  "treatment" = "Treatment",
  "genotype"  = "Genotype"
)

COVARIATES <- c(
  "age"                  = "Age",
  "weight"               = "Weight",
  "blood_glucose_before" = "Blood Glucose (Before)",
  "blood_glucose_after"  = "Blood Glucose (After)",
  "uacr_before"          = "uACR (Before)",
  "uacr_after"           = "uACR (After)",
  "ugcr_before"          = "uGCR (Before)",
  "ugcr_after"           = "uGCR (After)"
)
Sample distance
Calculate sample distance matrix:

dist.dist <- dist(t(assay(vsd)))

dist.matrix <- as.matrix(dist.dist)
Factors
Plot sample distances with factors displayed:

ann.row <- as.data.frame(colData(vsd)[ , names(FACTORS)])

colnames(ann.row) <- FACTORS

pheatmap(
  mat = dist.matrix,
  clustering_distance_rows = dist.dist,
  clustering_distance_cols = dist.dist,
  annotation_row = ann.row,
  show_colnames = FALSE
)


Covariates
Plot sample distances with covariates displayed:

ann.row <- as.data.frame(colData(vsd)[ , names(COVARIATES)])

colnames(ann.row) <- COVARIATES

pheatmap(
  mat = dist.matrix,
  border_color = NA,
  clustering_distance_rows = dist.dist,
  clustering_distance_cols = dist.dist,
  annotation_row = ann.row,
  annotation_legend = FALSE,
  show_colnames = FALSE
)


Variance explained
Genes
Compute, for each gene, the percentage of variance that is explained by one or more variables of interest:

var.fct <- setdiff(names(FACTORS), "group")

var.cov <- names(COVARIATES)

var.dat <- colData(vsd)[, c(var.fct, var.cov)]

var.exp <- getVarianceExplained(assay(vsd), variables = var.dat)
Plot explanatory variables by percentage of variance explained:

var.ord <- names(sort(colSums(var.exp, na.rm = TRUE), decreasing = TRUE)) # order by variance explained

plotExplanatoryVariables(var.exp) + 
  scale_colour_discrete(
    breaks = var.ord,
    labels = c(FACTORS, COVARIATES)[var.ord]
  ) + 
  labs(colour = "Variable") + 
  theme_bw()


PCs
Compute, for each principal component, the percentage of variance that is explained by one or more variables of interest:

row.var <- rowVars(assay(vsd))

row.ind <- order(row.var, decreasing=TRUE)[seq_len(min(500, length(row.var)))]

pca <- prcomp(t(assay(vsd)[row.ind, ]))

var.fct <- setdiff(names(FACTORS), "group")

var.cov <- names(COVARIATES)

var.dat <- colData(vsd)[, c(var.fct, var.cov)]

var.exp <- getVarianceExplained(t(pca$x), variables = var.dat)
Plot the explanatory PCs for each variable:

var.top <- head(var.exp, n = 10)

var.dat <- reshape2::melt(var.top)

colnames(var.dat) <- c("PC", "Variable", "R2")

var.ord <- names(sort(colSums(var.top, na.rm = TRUE), decreasing = TRUE)) # order by variance explained

ggplot(var.dat, aes(PC, Variable, fill = R2)) + 
  geom_tile() + 
  scale_fill_viridis_c(name = bquote(R^2)) +
  scale_y_discrete(
    breaks = var.ord,
    labels = c(FACTORS, COVARIATES)[var.ord]
  )


PCA plot
Factors
Visualize factors on PCA plot:

for (n in names(FACTORS)) {

  cat("#####", FACTORS[n], "\n")
  
  pca <- plotPCA(vsd, intgroup = n) + labs(colour = FACTORS[n])
  
  print(pca)
  
  cat("\n\n")
  
}
Group
Disease
Treatment
Genotype


Covariates
Visualize covariates on PCA plot:

for (n in names(COVARIATES)) {

  cat("#####", COVARIATES[n], "\n")
  
  pca <- plotPCA(vsd, intgroup = n) + labs(colour = COVARIATES[n])
  
  print(pca)
  
  cat("\n\n")
  
}
Age
Weight
Blood Glucose (Before)
Blood Glucose (After)
uACR (Before)
uACR (After)
uGCR (Before)
uGCR (After)


MDS plot
Calculate MDS matrix:

dist.matrix <- as.matrix(dist(t(assay(vsd))))

dist.data <- data.frame(cmdscale(dist.matrix))

colnames(dist.data) <- c("MDS1", "MDS2")

dist.data <- cbind(dist.data, as.data.frame(colData(vsd)))
Factors
Visualize factors on MDS plot:

for (n in names(FACTORS)) {

  cat("#####", FACTORS[n], "\n")
  
  mds <- ggplot(dist.data, aes_string("MDS1", "MDS2", colour = n)) + geom_point(size = 3) + labs(colour = FACTORS[n])
  
  print(mds)
  
  cat("\n\n")
  
}
Group
Disease
Treatment
Genotype


Covariates
Visualize covariates on MDS plot:

for (n in names(COVARIATES)) {

  cat("#####", COVARIATES[n], "\n")
  
  mds <- ggplot(dist.data, aes_string("MDS1", "MDS2", colour = n)) + geom_point(size = 3) + labs(colour = COVARIATES[n])
  
  print(mds)
  
  cat("\n\n")
  
}