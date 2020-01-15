configfile: "config.yaml"

localrules: all, whitelist, schpf_install

org = 'human' if config['human'] else 'mouse'
whitelist = 'gencode.vM19.annotation.gene_l1l2.pc_TRC_IGC.stripped.txt'  if org == 'mouse' else 'gencode.v29.annotation.gene_l1l2.pc_TRC_IGC.stripped.txt'
data_dir = '{}/{}'.format(config['data_dir'], config['name'])

rule all:
    input: 
        normalized=expand("{data_dir}/data/normalized.loom", data_dir=data_dir)

rule clean:
    input: config['loom_file']
    output: 
        figure=expand("{data_dir}/figures/cell_distributions.png", data_dir=data_dir),
        cleaned=expand("{data_dir}/data/cleaned.loom", data_dir=data_dir)
    params:
        cell_min_percentile=config['cell_min_percentile'],
        min_expr_counts_s=config['min_expr_counts_s'],
        min_expr_counts_u=config['min_expr_counts_u'],
        min_cells_express_s=config['min_cells_express_s'],
        min_cells_express_u=config['min_cells_express_u']
    conda: "envs/r_cellvariation.yaml"
    script: "scripts/clean.R"

rule normalize:
    input: expand("{data_dir}/data/cleaned.loom", data_dir=data_dir)
    output:
        sctransform=expand("{data_dir}/data/sctransform.RData", data_dir=data_dir),
        normalized=expand("{data_dir}/data/normalized.loom", data_dir=data_dir)
    conda: "envs/r_cellvariation.yaml"
    script: "scripts/normalize.R"

rule umap:
    input: expand("{data_dir}/data/normalized.loom", data_dir=data_dir)
    output: 
        figure=expand("{data_dir}/figures/umap.png", data_dir=data_dir)
    params:
        n_neighbors=config['n_neighbors'],
        metric=config['metric'],
        min_dist=config['min_dist'],
        color=config['umap_color'],
        legend=config['umap_legend'],
        psize=config['umap_psize']
    conda: "envs/cellvariation.yaml"
    script: "scripts/umap.py"

rule schpf_install:
    output: directory(expand('{dir}/scHPF', dir=config['data_dir']))
    conda: 'envs/schpf.yaml'
    shell:
        """
        git clone https://github.com/simslab/scHPF.git {config[data_dir]}/scHPF
        pip install --user {config[data_dir]}/scHPF
        """

rule schpf:
    input: 
        schpf=expand('{dir}/scHPF', dir=config['data_dir']),
        cleaned=expand("{data_dir}/data/cleaned.loom", data_dir=data_dir)
    output:
        prep_dir=directory(expand("{data_dir}/schpf/prep", data_dir=data_dir)),
        out_dir=directory(expand("{data_dir}/schpf/k{k}", data_dir=data_dir, 
            k=config['n_factors']))
    params:
        schpf=expand('{dir}/scHPF/bin/scHPF', dir=config['data_dir']),
        whitelist=expand("{dir}/scHPF/resources/{whitelist}", dir=config['data_dir'], whitelist=whitelist),
        k=config['n_factors']
    conda: 'envs/schpf.yaml'
    shell:
        """
        {params.schpf} prep -i "{input.cleaned}" -o "{output.prep_dir}" -w "{params.whitelist}"
        {params.schpf} train -i "{output.prep_dir}/filtered.mtx" -o "{output.out_dir}" -k "{params.k}" -t 5
        {params.schpf} score -m "{output.out_dir}/*.joblib" -o "{output.out_dir}" -g "{output.prep_dir}/genes.txt"
        """