input_table = pd.read_csv(config['taxids'], sep='\t')
input_table.columns = ['uniprot', 'taxid', 'mnemo']
input_dict = input_table.set_index('uniprot').T.to_dict()

codes = list(input_table['uniprot'])



rule plot_blens:
    input:
        blast=homodir+"/{seed}_blast.tsv",
        fs=homodir+"/{seed}_fs.tsv",
        trees=rules.get_unrooted_trees.output
    params: 
        eval_both=config['eval_both'],
        coverage=config['coverage'],
        max_seqs=config['max_seqs']
    output: outdir+"/plots/{seed}_distance.pdf"
    conda: "../envs/sp_R.yaml"
    script: "../scripts/analyze_bl.R"


checkpoint get_examples_ids:
    input: outdir+"/ids/{seed}_aln.ids"
    output: outdir+"/ids/{seed}_examples.ids"
    params: config["n_examples"]
    shell: '''
shuf -n {params} {input} > {output}    
'''

def seeds_examples(wildcards):
    checkpoint_output = checkpoints.get_examples_ids.get(**wildcards).output[0]
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.html", 
                      seed=wildcards.seed, i=seed_genes, mode=config["foldmason_set"])
    return outfiles

rule get_examples_report:
    input: seeds_examples
    output: directory(outdir+"/plots/{seed}_foldmason")
    shell: 'mkdir {output}; cp {input} {output}'

rule plot_examples:
    input: 
        trees=rules.get_unrooted_trees.output,
        ids=rules.get_examples_ids.output,
        taxidmap=homodir+"/db/taxidmap",
        reco=rules.merge_Notung.output,
        sptree=config['species_tree'],
        table=config['data_dir']+'meta/'+config["homology_dataset"]+'_uniprot_genomes.tsv',
        gff=expand(config['data_dir']+'gffs/{code}.gff', code=codes)
    params: 
        seed=config["seed"],
        struct_dir=config["data_dir"]+"structures/"
    output: outdir+"/plots/{seed}_examples.pdf"
    conda: "../envs/sp_R.yaml"
    script: "../scripts/plot_example.R"

### Comparisons

rule get_lineage:
    input: config['data_dir']+'meta/'+config["homology_dataset"]+'_uniprot_genomes.tsv'
    output: temp(outdir+"/reco/{seed}_lineage.tsv")
    conda: "../envs/sp_utils.yaml"
    shell:'''
cut -f2,3 {input} | taxonkit reformat -I 1 | sed 's/;/,/g' | awk 'NR>1' > {output}
'''

# rule get_verticality:
#     input: 
#         trees=rules.get_unrooted_trees.output,
#         lineage=rules.get_lineage.output,
#         taxidmap=rules.make_blastdb.output.mapid
#         # sptree=config['species_tree']
#     output:
#         # ranger=temp(outdir+"/reco/{seed}_ranger.tsv"),
#         outdir+"/reco/{seed}_scores.tsv"
#     conda: "../envs/sp_python.yaml"
#     script:"../scripts/foldtree/compute_scores.py"

rule support_astral_pro:
    input:
        sptree=config['species_tree'],
        genemap=homodir+"/db/gene_species.map",
        trees=rules.get_unrooted_trees.output,
    output: 
        gt=outdir+"/reco/apro/{seed}_{mode}_{alphabet}_{model}_apro_input.nwk",
        st=outdir+"/reco/apro/{seed}_{mode}_{alphabet}_{model}_apro_support.nwk"
    log: outdir+"/log/apro/{seed}_{mode}_{alphabet}_{model}_apro.log"
    params: config['root']
    conda: "../envs/reco.yaml"
    shell:'''
awk '$2=="{wildcards.mode}" && $3=="{wildcards.alphabet}" && $4=="{wildcards.model}"' {input} | \
cut -f5 > {output.gt}
astral-pro -c {input.sptree} -a {input.genemap} -u 2 -i {output.gt} -o {output.st} -C --root {params} 2> {log}
'''

rule disco:
    input:
        # sptree=config['species_tree'],
        genemap=homodir+"/db/gene_species.map",
        trees=rules.get_unrooted_trees.output,
    output: 
        gt=outdir+"/reco/disco/{seed}_{mode}_{alphabet}_{model}_disco_input.nwk",
        gt_out=outdir+"/reco/disco/{seed}_{mode}_{alphabet}_{model}_disco_output.nwk"
    log: outdir+"/log/disco/{seed}_{mode}_{alphabet}_{model}_disco.log"
    params: config['root']
    conda: "../envs/reco.yaml"
    shell:'''
awk '$2=="{wildcards.mode}" && $3=="{wildcards.alphabet}" && $4=="{wildcards.model}"' {input} | \
cut -f5 | nw_rename - {input.genemap} > {output.gt}
disco.py -i {output.gt} -o {output.gt_out}
'''


rule plot_trees:
    input: 
        sptree=config['species_tree'],
        # scores=rules.get_verticality.output,
        reco=rules.merge_Notung.output,
        trees=outdir+"/trees/{seed}_unrooted_trees.txt",
        mltrees=outdir+"/trees/{seed}_mltrees.txt",
        table=config['data_dir']+'meta/'+config["homology_dataset"]+'_uniprot_genomes.tsv',
        taxidmap=homodir+"/db/taxidmap",
        disco=expand(outdir+"/reco/disco/{seed}_{mode}_{comb}_disco_output.nwk", 
                     seed=config['seed'], mode=config["modes"], comb=config["combinations"]),
        apro_trees=expand(outdir+"/reco/apro/{seed}_{mode}_{comb}_apro_support.nwk", 
                     seed=config['seed'], mode=config["modes"], comb=config["combinations"])
    output:
        # apro=outdir+"/plots/{seed}_astral_pro.pdf",
        # model=outdir+"/plots/{seed}_trees.pdf",
        reco=outdir+"/plots/{seed}_discordance.pdf"
    conda: "../envs/sp_R.yaml"
    script: "../scripts/compare_sptree.R"


# rule plot_trees:
#     input: 
#         sptree=config['species_tree'],
#         # scores=rules.get_verticality.output,
#         reco=rules.merge_Notung.output,
#         trees=outdir+"/trees/{seed}_unrooted_trees.txt",
#         mltrees=outdir+"/trees/{seed}_mltrees.txt",
#         table=config['data_dir']+'meta/'+config["homology_dataset"]+'_uniprot_genomes.tsv',
#         taxidmap=homodir+"/db/taxidmap",
#         disco=expand(outdir+"/reco/disco/{seed}_{mode}_{comb}_disco_output.nwk", 
#                      seed=config['seed'], mode=config["modes"], comb=config["combinations"]),
#         apro_trees=expand(outdir+"/reco/{seed}_{mode}_{comb}_apro_support.nwk", 
#                      seed=config['seed'], mode=config["modes"], comb=config["combinations"])
#     output:
#         apro=outdir+"/plots/{seed}_astral_pro.pdf",
#         model=outdir+"/plots/{seed}_trees.pdf",
#         reco=outdir+"/plots/{seed}_discordance.pdf"
#     conda: "../envs/sp_R.yaml"
#     script: "../scripts/compare_trees.R"


# rule get_runstats:
#     input: seeds_unrooted_trees
#     output:
#         aln=outdir+"/stats/{seed}_aln.stats",
#         time=outdir+"/stats/{seed}_runtime.stats"
#     threads: 12
#     conda: "../envs/sp_utils.yaml"
#     shell:'''
# basedir=$(dirname {input} | rev | cut -f2- -d'/' | rev | sort -u)
# seqkit stats $basedir/*/*alg* -j {threads} -a -b -T | cut -f1,4,6,12 > {output.aln}

# datasetdir=$(echo $basedir | rev | cut -f3- -d'/' | rev)
# echo -e "dirname\\tbasename\\ts\\th:m:s\\tmax_rss\\tmax_vms\\tmax_uss\\tmax_pss\\tio_in\\tio_out\\tmean_load\\tcpu_time" > {output.time}
# find $datasetdir/benchmarks/ -type f -name "*txt" -printf "%h\\t%f\\t" -exec tail -1 {{}} \; >> {output.time}
# '''

# rule plot_runstats:
#     input: 
#         aln=rules.get_runstats.output.aln,
#         time=rules.get_runstats.output.time
#     output: outdir+"/plots/{seed}_runtime.pdf"
#     conda: "../envs/sp_R.yaml"
#     script: "../scripts/analyze_runtimes.R"


# rule plot_astral_pro:
#     input:
#         sptree=config['species_tree'],
#         table=config['data_dir']+'/meta/uniprot_genomes.tsv',
#         trees=expand(outdir+"/reco/{seed}_{mode}_{comb}_apro_support.nwk", 
#                      seed=config['seed'], mode=config["modes"], comb=config["combinations"])
#     output: outdir+"/plots/{seed}_astral_pro.pdf"
#     conda: "../envs/sp_R.yaml"
#     script: "../scripts/analyze_apro.R"
# reco_dir=$(dirname {input} | sort -u)

