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

def seeds_unrooted_trees(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_{mode}_{comb}.nwk", 
                      seed=wildcards.seed, i=seed_genes, mode=config["modes"], comb=config["combinations"])
    return outfiles


checkpoint get_unrooted_trees:
    input: 
        seeds_unrooted_trees
    output:
        trees=outdir+"/trees/{seed}_unrooted_trees.txt"
    shell:'''
> {output.trees}
for file in {input}; do
    echo -e "$(basename $file ".nwk" | tr '_' '\\t')\\t$(cat $file)" >> {output.trees}
done
'''


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
echo -e "$(basename $file | cut -f1,2,3,4 -d'_')\\t$(awk 'NR==1' $file | cut -f2,5)"
done > {output}
'''


# def seeds_rangers(wildcards):
#     checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
#     with open(checkpoint_output) as all_genes:
#         seed_genes = [gn.strip() for gn in all_genes]
#     outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_{mode}_{comb}_ranger.txt", 
#                       seed=wildcards.seed, i=seed_genes, mode=config["modes"], 
#                       comb=config["combinations"])
#     return outfiles

# rule merge_Ranger:
#     input: seeds_rangers
#     output: 
#         DTLs=outdir+"/reco/{seed}_DTL.tsv",
#         mrca=outdir+"/reco/{seed}_mrca.tsv"
#     shell: '''
# echo -e "id\\ttargets\\talphabet\\tmodel\\treco_score\\tdups\\ttransfers\\tlosses" > {output.DTLs}
# echo -e "id\\tn\\tmrca" > {output.mrca}

# for file in {input}; do
#     DTL=$(grep reco $file | grep -E "[0-9]+" -o | paste -s -d'\\t')
#     echo -e "$(basename $file "_ranger.txt" | tr '_' '\\t')\\t$DTL" >> {output.DTLs}
#     grep "Duplication, Mapping" $file | awk 'NF>1{{print $NF}}' | \
#     sort | uniq -c | sed -E 's/^ *//; s/ /\\t/' | \
#     awk -v file=$(basename $file "_ranger.txt") '{{print file"\\t"$0}}' >> {output.mrca} || true 
# done
# '''

# # first checkpoint and get the uniq ids in the blast or foldseek result
# checkpoint geneids_todo:
#     input:
#         blast=homodir+"/{seed}_blast_filtered.tsv",
#         fs=homodir+"/{seed}_fs_filtered.tsv"
#     output:
#         blast=outdir+"/ids/{seed}_blast.ids",
#         fs=outdir+"/ids/{seed}_fs.ids",
#         common=outdir+"/ids/{seed}_common.ids"
#     params: config['test_seeds']
#     shell:'''
# cut -f1 {input.blast} | sort -u | sort > {output.blast}
# cut -f1 {input.fs} | sort -u | sort > {output.fs}
# comm -12 {output.blast} {output.fs} | shuf -n {params} > {output.common}
# '''


# # get the targets ids and the subset homology table for fs and blast 
# rule get_ids:
#     input: homodir+"/{seed}_{method}_filtered.tsv"
#     output: 
#         txt=outdir+"/seeds/{seed}/{i}/{i}.{method}",
#         ids=outdir+"/seeds/{seed}/{i}/{i}_{method}.ids"
#         # top_ids=outdir+"/seeds/{seed}/{i}/{i}_{method}.top"
#     params: 
#         eval_both=config['eval_both'],
#         coverage=config['coverage'],
#         max_seqs=config['max_seqs']
#     shell:'''
# mkdir -p $(dirname {output.txt})

# awk '$1=="AF-{wildcards.i}-F1"' {input} > {output.txt}

# seed_gene=$(cut -f1 {output.txt} | sort -u)
# echo $seed_gene > {output.ids}

# n_hits=$(wc -l < {output.txt})
# if [ $n_hits -gt 1 ]; then
#     cut -f2 {output.txt} | grep -v -w $seed_gene | awk 'NR < {params.max_seqs}' >> {output.ids}
# fi
# sort {output.ids} -o {output.ids}
# '''
# foldseek results are sorted by bitScore * sqrt(alnlddt * alntmscore) https://github.com/steineggerlab/foldseek/issues/85

# rule get_common_ids:
#     input: 
#         blast=outdir+"/seeds/{seed}/{i}/{i}_blast.ids",
#         fs=outdir+"/seeds/{seed}/{i}/{i}_fs.ids"
#     output: 
#         outdir+"/seeds/{seed}/{i}/{i}_common.ids"
#     shell:'''
# comm -12 {input.blast} {input.fs} > {output}
# '''


# rule get_union_ids:
#     input: 
#         blast=outdir+"/seeds/{seed}/{i}/{i}_blast.ids",
#         fs=outdir+"/seeds/{seed}/{i}/{i}_fs.ids"
#     output: 
#         outdir+"/seeds/{seed}/{i}/{i}_union.ids"
#     shell:'''
# cat {input.blast} {input.fs} | sort | uniq > {output}
# '''

# def seeds_homology(wildcards):
#     checkpoint_output = checkpoints.geneids_todo.get(**wildcards).output.common
#     with open(checkpoint_output) as all_genes:
#         seed_genes = [gn.strip() for gn in all_genes]
#         parsed_seed_genes = [gn.split('-')[1] for gn in seed_genes]
#     return expand(outdir+"/seeds/{seed}/{i}/{i}_common.ids", seed=wildcards.seed, i=parsed_seed_genes)

# def common_trees(wildcards):
#     checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
#     with open(checkpoint_output) as all_genes:
#         seed_genes = [gn.strip() for gn in all_genes]
#     outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_common_{comb}.nwk", 
#                       seed=wildcards.seed, i=seed_genes, comb=config["combinations"])
#     return outfiles

# def fs_trees(wildcards):
#     checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
#     with open(checkpoint_output) as all_genes:
#         seed_genes = [gn.strip() for gn in all_genes]
#     outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_fs_{comb}.nwk", 
#                       seed=wildcards.seed, i=seed_genes, comb=config["combinations"])
#     return outfiles

# def blast_trees(wildcards):
#     checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
#     with open(checkpoint_output) as all_genes:
#         seed_genes = [gn.strip() for gn in all_genes]
#     outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_blast_{comb}.nwk", 
#                       seed=wildcards.seed, i=seed_genes, comb=config["combinations"])
#     return outfiles


# def union_trees(wildcards):
#     checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
#     with open(checkpoint_output) as all_genes:
#         seed_genes = [gn.strip() for gn in all_genes]
#     outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_union_{comb}.nwk", 
#                       seed=wildcards.seed, i=seed_genes, comb=config["combinations"])
#     return outfiles

# checkpoint check_union_trees:
#     input: union_trees
#     output:
#         trees=temp(outdir+"/trees/{seed}_union_trees.tmp")
#     shell:'''
# > {output.trees}
# for file in {input}; do
#     echo -e "$(basename $file ".nwk" | tr '_' '\\t')\\t$(cat $file)" >> {output.trees}
# done
# '''

# rule union_trees:
#     input: outdir+"/trees/{seed}_union_trees.tmp"
#     output: outdir+"/trees/{seed}_union_trees.txt"
#     shell: 'mv {input} {output}'

# checkpoint check_common_trees:
#     input: common_trees
#     output:
#         trees=temp(outdir+"/trees/{seed}_common_trees.tmp")
#     shell:'''
# > {output.trees}
# for file in {input}; do
#     echo -e "$(basename $file ".nwk" | tr '_' '\\t')\\t$(cat $file)" >> {output.trees}
# done
# '''

# rule common_trees:
#     input: outdir+"/trees/{seed}_common_trees.tmp"
#     output: outdir+"/trees/{seed}_common_trees.txt"
#     shell: 'mv {input} {output}'


# checkpoint check_fs_trees:
#     input: fs_trees
#     output:
#         trees=temp(outdir+"/trees/{seed}_fs_trees.tmp")
#     shell:'''
# > {output.trees}
# for file in {input}; do
#     echo -e "$(basename $file ".nwk" | tr '_' '\\t')\\t$(cat $file)" >> {output.trees}
# done
# '''

# rule fs_trees:
#     input: outdir+"/trees/{seed}_fs_trees.tmp"
#     output: outdir+"/trees/{seed}_fs_trees.txt"
#     shell: 'mv {input} {output}'


# checkpoint check_blast_trees:
#     input: blast_trees
#     output:
#         trees=temp(outdir+"/trees/{seed}_blast_trees.tmp")
#     shell:'''
# > {output.trees}
# for file in {input}; do
#     echo -e "$(basename $file ".nwk" | tr '_' '\\t')\\t$(cat $file)" >> {output.trees}
# done
# '''

# rule blast_trees:
#     input: outdir+"/trees/{seed}_blast_trees.tmp"
#     output: outdir+"/trees/{seed}_blast_trees.txt"
#     shell: 'mv {input} {output}'
