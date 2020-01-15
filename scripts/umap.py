import umap
import loompy
import seaborn

with loompy.connect(snakemake.input[0]) as ds:
    counts = ds[:, :]
    # gene_info = pd.DataFrame(dict(ds.row_attrs.items()))
    # cell_info = pd.DataFrame(dict(ds.col_attrs.items()))

umap_model = umap.UMAP(n_neighbors=snakemake.params['n_neighbors'],
                       metric=snakemake.params['metric'],
                       min_dist=snakemake.params['min_dist'])
embedding = umap_model.fit_transform(counts.T)

with loompy.connect(snakemake.input[0]) as ds:
    ds.ca.umap1 = embedding[:, 0]
    ds.ca.umap2 = embedding[:, 1]

p = sns.scatterplot(embedding[:, 0], embedding[:, 1], 
                hue=snakemake.params['color'],
                legend=snakemake.params['legend'],
                s=snakemake.params['psize'],
                linewidth=0)
plt.axis('off')
p.savefig(snakemake.output['figure'])
