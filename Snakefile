localrules: all, whitelist, schpf_install

org = 'human' if config['human'] else 'mouse'
whitelist = 'gencode.vM19.annotation.gene_l1l2.pc_TRC_IGC.stripped.txt' if org == 'mouse' else 'gencode.v29.annotation.gene_l1l2.pc_TRC_IGC.stripped.txt'
data_dir = '{}/{}'.format(config['data_dir'], config['name'])

rule clean:
    input: config['loom_file']
    output: 
        figure=f'{data_dir}/figures/cell_distributions.png',
        cleaned=f'{data_dir}/data/cleaned.loom'
    params:
        cell_min_percentile=config['cell_min_percentile'],
        min_expr_counts_s=config['min_expr_counts_s'],
        min_expr_counts_u=config['min_expr_counts_u'],
        min_cells_express_s=config['min_cells_express_s'],
        min_cells_express_u=config['min_cells_express_u']
    conda: 'envs/r_cellvariation.yaml'
    script: 'scripts/clean.R'

rule normalize:
    input: f'{data_dir}/data/cleaned.loom'
    output:
        sctransform=f'{data_dir}/data/sctransform.RData',
        normalized=f'{data_dir}/data/normalized.loom'
    conda: 'envs/r_cellvariation.yaml'
    script: 'scripts/normalize.R'

rule umap:
    input: f'{data_dir}/data/normalized.loom',
    output: f'{data_dir}/figures/umap.svg'
    params:
        n_neighbors=config['n_neighbors'],
        metric=config['metric'],
        min_dist=config['min_dist'],
        color=config['umap_color'],
        legend=config['umap_legend'],
        psize=config['umap_psize']
    conda: 'envs/cellvariation.yaml'
    script: 'scripts/umap.py'

rule schpf_install:
    output: directory(f'{config["data_dir"]}/scHPF')
    conda: 'envs/schpf.yaml'
    shell:
        """
        git clone https://github.com/simslab/scHPF.git {config[data_dir]}/scHPF
        pip install --user {config[data_dir]}/scHPF
        """

rule schpf:
    input: 
        schpf=f'{config["data_dir"]}/scHPF',
        cleaned=f'{data_dir}/data/cleaned.loom'
    output:
        prep_dir=directory(f'{data_dir}/schpf/prep'),
        out_dir=directory(f'{data_dir}/schpf/k{config["n_factors"]}')
    params:
        schpf=f'{config["data_dir"]}/scHPF/bin/scHPF',
        whitelist=f'{config["data_dir"]}/scHPF/resources/{whitelist}',
        k=config['n_factors']
    conda: 'envs/schpf.yaml'
    shell:
        """
        {params.schpf} prep -i "{input.cleaned}" -o "{output.prep_dir}" -w "{params.whitelist}"
        {params.schpf} train -i "{output.prep_dir}/filtered.mtx" -o "{output.out_dir}" -k "{params.k}" -t 5
        {params.schpf} score -m `ls {output.out_dir}/*.joblib` -o "{output.out_dir}" -g "{output.prep_dir}/genes.txt"
        """

rule velocity:
    input: f'{data_dir}/data/normalized.loom'
    output:
        vlm=f'{data_dir}/data/vlm.hdf5'
    params:
        k=config['k_imputation'],
        cluster=config['cluster_name'],
        fit_offset=config['fit_offset'],
        transform=config['trans_prob_transform'],
        nn=config['trans_prob_nn'],
        scaling=config['shift_scaling'],
        smooth=config['arrows_smooth'],
        steps=config['arrows_steps']
    conda: 'envs/velocity.yaml'
    script: 'scripts/velocity.py'

rule pseudotime:
    input: 
        normalized=f'{data_dir}/data/normalized.loom',
        genes=f'{data_dir}/schpf/prep/genes.txt',
        gene_scores=f'{data_dir}/schpf/k{config["n_factors"]}/gene_score.txt',
        cell_scores=f'{data_dir}/schpf/k{config["n_factors"]}/cell_score.txt'
    output: f'{data_dir}/data/ouija.RData'
    params:
        factors=config['trajectory_factors'],
        ngenes=config['ngenes_per_factor'],
        inference=config['ouija_inference']
    conda: 'envs/r_cellvariation.yaml'
    script: 'scripts/ouija.R'