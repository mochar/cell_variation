import umap
import loompy
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

with loompy.connect(snakemake.input[0]) as ds:
    counts = ds[:, :]
    varying = ds.ra['varying']
    color = ds.ca[snakemake.params['color']]

umap_model = umap.UMAP(n_neighbors=snakemake.params['n_neighbors'],
                       metric=snakemake.params['metric'],
                       min_dist=snakemake.params['min_dist'])
embedding = umap_model.fit_transform(counts.T[:, varying == 1])

with loompy.connect(snakemake.input[0]) as ds:
    ds.ca.umap1 = embedding[:, 0]
    ds.ca.umap2 = embedding[:, 1]

p = sns.scatterplot(embedding[:, 0], embedding[:, 1], 
                hue=color, s=snakemake.params['psize'],
                legend='full' if snakemake.params['legend'] else False,
                linewidth=0)
plt.axis('off')
p.get_figure().savefig(snakemake.output[0], format='svg')
