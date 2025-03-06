homodir=config['outdir']+'homology/'+config['homology_dataset']
outdir=config['outdir']+'phylogeny/'+config['phylo_dataset']

seeds = pd.read_csv(config["test_seeds"], sep='\t', header=None)
seeds = list(set(seeds[0].tolist()))

rule get_ids:
    input:
        blast=homodir+"/{seed}_blast_filtered.tsv",
        fs=homodir+"/{seed}_fs_filtered.tsv",
        seed_file=config["test_seeds"]
    params:
        modes=config["modes"],
        max_seqs=config['max_seqs']
    output:
        expand(outdir+"/seeds/{{seed}}/{i}/{i}_{mode}.ids", i=seeds, mode=config["modes"])
    conda: "../envs/sp_R.yaml"
    script: "../scripts/get_sets.R"

checkpoint check_orphans:
    input: expand(outdir+"/seeds/{{seed}}/{i}/{i}_common.ids", i=seeds, mode=config["modes"])
    output:
        exclude=outdir+"/ids/{seed}_orphans.exclude",
        continue_aln=outdir+"/ids/{seed}_aln.ids"
    params: min_common=config["min_common"]
    shell:'''
> {output.exclude}
> {output.continue_aln}

for file in {input}; do
    n_hits=$(wc -l < $file)
    if [ $n_hits -lt {params.min_common} ]; then
        echo -e "$(basename $file | cut -f1 -d'_')\\tless than {params.min_common} common hits" >> {output.exclude}
    else
        echo "$(basename $file | cut -f1 -d'_')" >> {output.continue_aln}
    fi
done
'''

##### FINAL CHECK IF TREE WORKED #######

def seeds_treefiles(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_{mode}_{comb}.iqtree", 
                      seed=wildcards.seed, i=seed_genes, mode=config["modes"], comb=config["combinations_ML"])
    return outfiles


checkpoint check_ML_trees:
    input: seeds_treefiles
    output: outdir+"/trees/{seed}_mltrees.txt"
    shell:'''
echo -e "gene\\ttargets\\tModel\\tLogL\\tAIC\\tw-AIC\\tAICc\\tw-AICc\\tBIC\\tw-BIC" > {output}

for file in {input}; do
    gene=$(basename $file ".iqtree" | cut -f1 -d'_')
    targets=$(basename $file ".iqtree" | cut -f2 -d'_')
    if [[ $file =~ "comb_part" ]]; then
        continue
    fi
    model_line=$(grep -E '^(LG|GTR20|3DI|resources)' $file | awk 'NR==1' | awk -v OFS="\\t" '$1=$1' | sed 's/\\t[+-]\\t/\\t/g')
    echo -e "$gene\\t$targets\\t$(echo $model_line | tr ' ' '\\t')"
done | \
sed 's/resources\\/subst_matrixes\\/Q_mat_//g ; s/_Garg.txt//g' >> {output}
'''


def seeds_notung(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/reco/notung/{seed}/{i}_{mode}_{comb}_rnm.nwk.rooting.0.parsable.txt", 
                      seed=wildcards.seed, i=seed_genes, mode=config["modes"], 
                      comb=config["combinations"])
    return outfiles

rule merge_Notung:
    input: seeds_notung
    output: outdir+"/reco/{seed}_notung.tsv"
    shell: '''
for file in {input}; do
echo -e "$(basename $file | cut -f1-4 -d'_' |tr '_' '\\t')\\t$(awk 'NR==1' $file | cut -f2,5)"
done > {output}
'''


def common_trees(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_common_{comb}.nwk", 
                      seed=wildcards.seed, i=seed_genes, comb=config["combinations"])
    return outfiles

def fs_trees(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_fs_{comb}.nwk", 
                      seed=wildcards.seed, i=seed_genes, comb=config["combinations"])
    return outfiles

def blast_trees(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_blast_{comb}.nwk", 
                      seed=wildcards.seed, i=seed_genes, comb=config["combinations"])
    return outfiles


def union_trees(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_union_{comb}.nwk", 
                      seed=wildcards.seed, i=seed_genes, comb=config["combinations"])
    return outfiles


def seeds_unrooted_trees(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/trees/{seed}_{mode}_trees.txt", 
                      seed=wildcards.seed, mode=config["modes"])
    return outfiles


checkpoint check_union_trees:
    input: union_trees
    output: temp(outdir+"/trees/{seed}_union_trees.tmp")
    shell:'''
> {output}
for file in {input}; do
    echo -e "$(basename $file ".nwk" | tr '_' '\\t')\\t$(cat $file)" >> {output}
done
'''

rule union_trees:
    input: outdir+"/trees/{seed}_union_trees.tmp"
    output: outdir+"/trees/{seed}_union_trees.txt"
    shell: 'mv {input} {output}'


checkpoint check_common_trees:
    input: common_trees
    output: temp(outdir+"/trees/{seed}_common_trees.tmp")
    shell:'''
> {output}
for file in {input}; do
    echo -e "$(basename $file ".nwk" | tr '_' '\\t')\\t$(cat $file)" >> {output}
done
'''

rule common_trees:
    input: outdir+"/trees/{seed}_common_trees.tmp"
    output: outdir+"/trees/{seed}_common_trees.txt"
    shell: 'mv {input} {output}'


checkpoint check_fs_trees:
    input: fs_trees
    output: temp(outdir+"/trees/{seed}_fs_trees.tmp")
    shell:'''
> {output}
for file in {input}; do
    echo -e "$(basename $file ".nwk" | tr '_' '\\t')\\t$(cat $file)" >> {output}
done
'''

rule fs_trees:
    input: outdir+"/trees/{seed}_fs_trees.tmp"
    output: outdir+"/trees/{seed}_fs_trees.txt"
    shell: 'mv {input} {output}'


checkpoint check_blast_trees:
    input: blast_trees
    output: temp(outdir+"/trees/{seed}_blast_trees.tmp")
    shell:'''
> {output}
for file in {input}; do
    echo -e "$(basename $file ".nwk" | tr '_' '\\t')\\t$(cat $file)" >> {output}
done
'''

rule blast_trees:
    input: outdir+"/trees/{seed}_blast_trees.tmp"
    output: outdir+"/trees/{seed}_blast_trees.txt"
    shell: 'mv {input} {output}'


rule get_unrooted_trees:
    input: seeds_unrooted_trees
    output: outdir+"/trees/{seed}_unrooted_trees.txt"
    shell:'''
cat {input} > {output}
'''
# for file in {input}; do
#     echo -e "$(basename $file ".nwk" | tr '_' '\\t')\\t$(cat $file)" >> {output}
# done
