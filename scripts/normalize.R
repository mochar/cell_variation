library(loomR)
library(Seurat)
library(dplyr)

## Load in data
lfile <- connect(filename = snakemake@input[[1]], mode = "r+")
rowData <- lfile$get.attribute.df(MARGIN = 1, row.names = 'Accession')
colData <- lfile$get.attribute.df(MARGIN = 2, col.names = 'CellID')

# Spliced
counts <- t(lfile$layers$spliced[, ])
rownames(counts) <- rowData$Accession
colnames(counts) <- colData$CellID

# Unspliced
counts_u <- t(lfile$layers$unspliced[, ])
rownames(counts_u) <- rowData$Accession
colnames(counts_u) <- colData$CellID

# Ambiguous
counts_a <- t(lfile$layers$ambiguous[, ])
rownames(counts_a) <- rowData$Accession
colnames(counts_a) <- colData$CellID

lfile$close_all()

## Unspliced normalization
srt_unspliced <- CreateSeuratObject(counts = counts_u)
srt_unspliced <- SCTransform(srt_unspliced, return.only.var.genes = TRUE)

## Spliced normalization
srt <- CreateSeuratObject(counts = counts)
srt <- SCTransform(srt, return.only.var.genes = TRUE)

## Add gene attributes of normalization to row data
rowData$residual_variance_unspliced <- Misc(srt_unspliced[['SCT']], slot = 'vst.out')$gene_attr[rownames(rowData), 'residual_variance']
rowData$residual_variance <- Misc(srt[['SCT']], slot = 'vst.out')$gene_attr[rownames(rowData), 'residual_variance']
rowData$varying <- rowData$Accession %in% VariableFeatures(srt)

## Save sctransform objects
sct_spliced <- Misc(srt[['SCT']], slot = 'vst.out')
sct_unspliced <- Misc(srt_unspliced[['SCT']], slot = 'vst.out')
save(sct_spliced, sct_unspliced, file = snakemake@output[['sctransform']])

## Unspliced matrix not same dim as spliced because rowsums == 0 which
## leads to error in loomR::create. So: change values in count_u to that in
## srt_unspliced@count, and pass that to loomR::create.
counts_u_export <- counts_u
counts_u_export[rownames(srt_unspliced[['SCT']]@counts), ] <- srt_unspliced[['SCT']]@counts %>% as.matrix

## Export normalized loom
lfile <- loomR::create(filename = snakemake@output[['normalized']],
                    #    max.size = '3gb',
                       feature.attrs = rowData,
                       cell.attrs = colData,
                       transpose = TRUE,
                    #    data = (srt[['SCT']]@counts + counts_u_export + counts_a) %>% unname,
                       data = srt[['SCT']]@counts %>% unname,
                       layers = list(spliced = srt[['SCT']]@counts %>% unname,
                                     unspliced = counts_u_export %>% unname,
                                     ambiguous = counts_a %>% unname))