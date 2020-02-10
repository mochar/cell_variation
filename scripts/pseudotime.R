library(loomR)
library(dplyr)
library(SingleCellExperiment)
library(tidyr)

options(mc.cores = parallel::detectCores())

## Need to install Ouija using devtools 
if (!require('ouija')) {
    library(devtools)
    devtools::install_github('kieranrcampbell/ouija', local = FALSE, 
                             args = "--preclean", build_vignettes = TRUE)
}
library(ouija)

## Load in counts and metadata
lfile <- connect(filename = snakemake@input[['normalized']], mode = 'r+')
rowData <- lfile$get.attribute.df(MARGIN = 1, row.names = 'Accession')
colData <- lfile$get.attribute.df(MARGIN = 2, col.names = 'CellID')
counts <- t(lfile$matrix[, ])
rownames(counts) <- rowData$GENE
colnames(counts) <- colData$CellID
lfile$close_all()

## Load in scHPF factors to retrieve list of genes
cell_scores <- read.csv(snakemake@input[['cell_scores']], 
                        sep = '\t', header = FALSE)
rownames(cell_scores) <- colData$CellID

gene_names_ids <- read.csv(snakemake@input[['genes']],
                           sep = '\t', header = FALSE)
gene_names_ids <- gene_names_ids[, c(1, 2, 5)]
colnames(gene_names_ids) <- c('Accession', 'Gene', 'GENE')

gene_scores <- read.csv(snakemake@output[['gene_scores']],
                        sep = '\t', header = FALSE)
rownames(gene_scores) <- gene_names_ids$Accession

## Retrieve top genes
factors <- snakemake@params[['factors']]
ngenes <- snakemake@params[['ngenes']]
genes <- sapply(factors, 
                function(k) gene_names_ids[order(gene_scores[, k+1]) %>% tail(ngenes), 'GENE'] %>% as.character)
genes <- as.character(genes)

## Run Ouija
sce <- SingleCellExperiment(assays = list(logcounts = counts[genes, ]))
oui <- ouija(sce, inference_type = snakemake@params[['inference']])

## Add pseudotime column to loom file
tmap <- map_pseudotime(oui)
lfile <- connect(filename = snakemake@input[['normalized']], mode = 'r+')
lfile$add.col.attribute(list(pseudotime = tmap), overwrite = TRUE)
lfile$close_all()

## Export ouija object
save(oui, file = snakemake@output[[1]])