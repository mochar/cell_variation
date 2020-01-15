library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

## Need to install loomR using devtools as conda package is outdated
if (!require('loomR')) {
    library(devtools)
    devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
}

library(loomR)

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

## For technical reasons downstream useful to have uppercase gene names
rowData$GENE <- toupper(rowData$Gene)

## Measures of cell quality
cell_sizes <- data.frame(spliced = apply(counts, 2, sum), 
                         unspliced = apply(counts_u, 2, sum))
cell_ngenes <- data.frame(spliced = apply(counts > 0, 2, sum), 
                          unspliced = apply(counts_u > 0, 2, sum))

colData$zero_fraction <- apply(counts == 0, 2, sum) / nrow(counts)
colData$log_total_umi <- log(apply(counts, 2, sum))

## Cutoff for low quality cells
cutoff <- quantile(cell_sizes$unspliced, probs = snakemake@params[['cell_min_percentile']])

## Plot cell quanitty distributions
cell_sizes_stacked <- stack(cell_sizes)
p_count <- ggplot(cell_sizes_stacked, aes(x = values)) + 
    geom_histogram(binwidth = 100) +
    facet_wrap(~ind) +
    geom_vline(data = filter(cell_sizes_stacked, ind=='unspliced'), 
               aes(xintercept = cutoff), colour = 'red') +
    labs(x = 'Total counts per droplet')
p_feature <- ggplot(cell_ngenes %>% stack, aes(x = values)) + 
    geom_histogram(binwidth = 100) +
    facet_wrap(~ind) +
    labs(x = 'Total features per droplet')
p <- plot_grid(p_count, p_feature, nrow = 2)
save_plot(snakemake@output[['figure']], p, base_height = 5.0)

## Filter low quality cells based on cutoff
print(paste(ncol(counts), 'cells before filtering'))
keep_cells <- cell_sizes$unspliced > cutoff
counts <- counts[, keep_cells]
counts_u <- counts_u[, keep_cells]
counts_a <- counts_a[, keep_cells]
colData <- colData[keep_cells, ]
print(paste(ncol(counts), 'cells after filtering'))

## Filter low quality genes based on spliced matrix
gene_counts <- data.frame(total_counts = apply(counts, 1, sum),
                          num_cells = apply(counts > 0, 1, sum))

print(paste(nrow(counts), 'genes before filtering (spliced)'))
keep_genes <- (gene_counts$total_counts >= snakemake@params[['min_expr_counts_s']]) & 
              (gene_counts$num_cells > snakemake@params[['min_cells_express_s']])
counts <- counts[keep_genes, ]
counts_u <- counts_u[keep_genes, ]
counts_a <- counts_a[keep_genes, ]
rowData <- rowData[keep_genes, ]
print(paste(nrow(counts), 'genes after filtering (spliced)'))

## Filter low quality genes based on unspliced matrix
gene_counts_u <- data.frame(total_counts = apply(counts_u, 1, sum),
                            num_cells = apply(counts_u > 0, 1, sum))

print(paste(nrow(counts), 'genes before filtering (unspliced)'))
keep_genes <- (gene_counts_u$total_counts >= snakemake@params[['min_expr_counts_u']]) & 
              (gene_counts_u$num_cells > snakemake@params[['min_cells_express_u']])
counts <- counts[keep_genes, ]
counts_u <- counts_u[keep_genes, ]
counts_a <- counts_a[keep_genes, ]
rowData <- rowData[keep_genes, ]
print(paste(nrow(counts), 'genes after filtering (spliced)'))

## Export
lfile <- loomR::create(filename = snakemake@output[['cleaned']],
                    #    max.size = '3gb',
                       feature.attrs = rowData,
                       cell.attrs = colData,
                       transpose = TRUE,
                    #    data = (counts + counts_u + counts_a) %>% unname,
                       data = counts %>% unname,
                       layers = list(spliced = counts %>% unname,
                                     unspliced = counts_u %>% unname,
                                     ambiguous = counts_a %>% unname))
lfile$close_all()