import velocyto as vcy

vlm = vcy.VelocytoLoom(snakemake.input[0])

# Set the clusters
vlm.set_clusters(vlm.ca[snakemake.params['cluster']])
vlm.ts = np.column_stack([vlm.ca['umap1'], vlm.ca['umap2']])

# Normalize and transformation
vlm.S_sz = vlm.S
vlm.U_sz = vlm.U
vlm.S_norm = np.log1p(vlm.S_sz)
vlm.U_norm = np.log1p(vlm.U_sz)

# PCA
vlm.perform_PCA(n_components=100)
n_comps = np.where(np.diff(np.diff(np.cumsum(vlm.pca.explained_variance_ratio_))>0.002))[0][0]
# plt.plot(np.cumsum(vlm.pca.explained_variance_ratio_)[:100])
# plt.axvline(n_comps, c="k")

# Imputation / pooling
k = snakemake.params['k']
vlm.knn_imputation(n_pca_dims=n_comps, k=k, n_jobs=16)

# Estimate velocity
vlm.fit_gammas(fit_offset=snakemake.params['fit_offset'])
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.) 

# Transition probablities
vlm.estimate_transition_prob(psc=1, transform=snakemake.params['transform'], 
                             n_neighbors=snakemake.params['nn'])
vlm.calculate_embedding_shift(expression_scaling=snakemake.params['scaling'])
vlm.calculate_grid_arrows(smooth=snakemake.params['smooth'], 
                          steps=(snakemake.params['steps'], snakemake.params['steps']), 
                          n_neighbors=snakemake.params['arrows_nn'])

# Store vlm
vlm.to_hdf5(snakemake.output['vlm'])