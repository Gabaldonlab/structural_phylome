import glob
import pandas as pd
import snakemake.utils

# User-specified dataset-specific config file
dataset_config_file = config.get("dataset_config_file", None)

# Load fixed config parameters
configfile: "config/params.yaml"

# Check if the dataset-specific config file is provided
if dataset_config_file is not None:
    # Load dataset-specific parameters
    configfile: dataset_config_file

include: 'rules/prepare.smk'
include: 'rules/structure.smk'
include: 'rules/sequence.smk'


outdir=config['outdir']+'homology/'+config['homology_dataset']

plots = ["homology", "saturation", "singletons"]

rule all:
    input:
        # expand(outdir+"/{seed}_blast_brh.tsv", seed=config['seed']),
        # expand(outdir+"/{seed}_fs_brh.tsv", seed=config['seed']),
        expand(outdir+"/plots/{seed}_{plot}.pdf", seed=config['seed'], plot=plots),
        outdir+"/db/all_seqs_fsdb_ss.fa",
        outdir+'/db/db_num.tsv',
        outdir+"/db/gene_species.map",
        outdir+"/db/taxidmap_sps",
        expand(outdir+"/{seed}_{method}_filtered.tsv", seed=config['seed'], method = ["blast", "fs"])


rule plot_evalues:
    input:
        table=config['data_dir']+'meta/'+config["homology_dataset"]+'_uniprot_genomes.tsv',
        taxidmap=outdir+"/db/taxidmap",
        blast=outdir+"/{seed}_blast.tsv",
        # blast_brh=rules.blast_brh.output,
        self_blast=outdir+"/allvall/{seed}_{seed}_blast.tsv",
        fs=outdir+"/{seed}_fs.tsv",
        # fs_brh=rules.foldseek_brh.output,
        self_fs=outdir+"/allvall/{seed}_{seed}_fs.tsv"
    params: 
        eval_both=config['eval_both'],
        max_seqs=config['max_seqs']
    output: 
        eda=outdir+"/plots/{seed}_homology.pdf",
        saturation=outdir+"/plots/{seed}_saturation.pdf"
    conda: "./envs/sp_R.yaml"
    script: "./scripts/compare_sampling.R"


rule plot_evalues_trees:
    input:
        table=config['data_dir']+'meta/'+config["homology_dataset"]+'_uniprot_genomes.tsv',
        # groups=config['taxons_file'],
        taxidmap=outdir+"/db/taxidmap",
        blast=outdir+"/{seed}_blast.tsv",
        fs=outdir+"/{seed}_fs.tsv"
        # ids=outdir+"/ids/{seed}_common.ids"
    params: 
        eval_both=config['eval_both'],
        coverage=config['coverage'],
        max_seqs=config['max_seqs']
    output: outdir+"/plots/{seed}_singletons.pdf"
    conda: "./envs/sp_R.yaml"
    script: "./scripts/compare_sampling_trees.R"
