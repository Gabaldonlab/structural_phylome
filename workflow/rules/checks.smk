outdir=config['outdir']+config['dataset']
alphabets = ['aa', '3Di']
models=["3Di", "QT", "LG", "GTR", "FT"]
combinations=["3Di_3Di", "aa_QT", "aa_LG", "3Di_GTR", "3Di_FT"]
combinations_ML=["3Di_3Di", "aa_LG", "3Di_GTR"]
# first checkpoint and get the uniq ids in the blast or foldseek result
checkpoint geneids_todo:
    input:
        blast=rules.blast.output,
        fs=rules.foldseek.output
    output:
        blast=outdir+"/ids/{seed}_blast.ids",
        fs=outdir+"/ids/{seed}_fs.ids",
        common=outdir+"/ids/{seed}_common.ids"
    params: config['test_seeds']
    shell:'''
cut -f1 {input.blast} | sort -u | sort > {output.blast}
cut -f1 {input.fs} | sort -u | sort > {output.fs}
comm -12 {output.blast} {output.fs} | shuf -n {params} > {output.common}
'''


# get the targets ids and the subset homology table for fs and blast 
rule get_ids:
    input: outdir+"/homology/{seed}_{method}.tsv"
    output: 
        txt=outdir+"/seeds/{seed}/{i}/{i}.{method}",
        ids=outdir+"/seeds/{seed}/{i}/{i}_{method}.ids"
    params: 
        eval_both=config['eval_both'],
        coverage=config['coverage'],
        max_seqs=config['max_seqs']
    shell:'''
mkdir -p $(dirname {output.txt})

if [ {wildcards.method} == "blast" ]; then
    awk '$1=="AF-{wildcards.i}-F1"' {input} | awk '$11<{params.eval_both} && $4/$15*100>{params.coverage} && $4/$16*100>{params.coverage}' > {output.txt}
else
    awk '$1=="AF-{wildcards.i}-F1"' {input} | awk '$11<{params.eval_both} && $17*100>{params.coverage} && $18*100>{params.coverage}' > {output.txt}
fi

seed_gene=$(cut -f1 {output.txt} | sort -u)
echo $seed_gene > {output.ids}
n_hits=$(wc -l < {output.txt})
if [ $n_hits -gt 1 ]; then
    sort -gr -k12 {output.txt} | cut -f2 | grep -v -w $seed_gene | sort -u | awk 'NR < {params.max_seqs}' >> {output.ids}
fi
sort {output.ids} -o {output.ids}
'''

rule get_common_ids:
    input: 
        blast=outdir+"/seeds/{seed}/{i}/{i}_blast.ids",
        fs=outdir+"/seeds/{seed}/{i}/{i}_fs.ids"
    output: 
        outdir+"/seeds/{seed}/{i}/{i}_common.ids"
    shell:'''
comm -12 {input.blast} {input.fs} > {output}
'''

def seeds_homology(wildcards):
    checkpoint_output = checkpoints.geneids_todo.get(**wildcards).output.common
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
        parsed_seed_genes = [gn.split('-')[1] for gn in seed_genes]
    return expand(outdir+"/seeds/{seed}/{i}/{i}_common.ids", seed=wildcards.seed, i=parsed_seed_genes, mode=modes)

checkpoint check_orphans:
    input: seeds_homology
    output:
        exclude=outdir+"/ids/{seed}_orphans.exclude",
        continue_aln=outdir+"/ids/{seed}_aln.ids",
    shell:'''
> {output.exclude}
> {output.continue_aln}

for file in {input}; do
    n_hits=$(wc -l < $file)
    if [ $n_hits -lt 4 ]; then
        echo -e "$(basename $file | cut -f1 -d'_')\tless than 4 common hits" >> {output.exclude}
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
                      seed=wildcards.seed, i=seed_genes, mode=modes, comb=combinations)
    return outfiles

def common_trees(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_common_{comb}.nwk", 
                      seed=wildcards.seed, i=seed_genes, mode=modes, comb=combinations)
    return outfiles

def fs_trees(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_fs_{comb}.nwk", 
                      seed=wildcards.seed, i=seed_genes, mode=modes, comb=combinations)
    return outfiles

def blast_trees(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_blast_{comb}.nwk", 
                      seed=wildcards.seed, i=seed_genes, mode=modes, comb=combinations)
    return outfiles


checkpoint check_common_trees:
    input: common_trees
    output:
        trees=temp(outdir+"/trees/{seed}_common_trees.tmp")
    shell:'''
> {output.trees}
for file in {input}; do
    id=$(basename $file ".nwk" | cut -f1 -d'_')
    targets=$(basename $file ".nwk" | cut -f2 -d'_')
    alphabet=$(basename $file ".nwk" | cut -f3 -d'_')
    model=$(basename $file ".nwk" | cut -f4 -d'_')
    echo -e "$id\\t$targets\\t$alphabet\\t$model\\t$(cat $file)" >> {output.trees}
done
'''

rule common_trees:
    input: outdir+"/trees/{seed}_common_trees.tmp"
    output: outdir+"/trees/{seed}_common_trees.txt"
    shell: 'mv {input} {output}'


checkpoint check_fs_trees:
    input: fs_trees
    output:
        trees=temp(outdir+"/trees/{seed}_fs_trees.tmp")
    shell:'''
> {output.trees}
for file in {input}; do
    id=$(basename $file ".nwk" | cut -f1 -d'_')
    targets=$(basename $file ".nwk" | cut -f2 -d'_')
    alphabet=$(basename $file ".nwk" | cut -f3 -d'_')
    model=$(basename $file ".nwk" | cut -f4 -d'_')
    echo -e "$id\\t$targets\\t$alphabet\\t$model\\t$(cat $file)" >> {output.trees}
done
'''

rule fs_trees:
    input: outdir+"/trees/{seed}_fs_trees.tmp"
    output: outdir+"/trees/{seed}_fs_trees.txt"
    shell: 'mv {input} {output}'


checkpoint check_blast_trees:
    input: blast_trees
    output:
        trees=temp(outdir+"/trees/{seed}_blast_trees.tmp")
    shell:'''
> {output.trees}
for file in {input}; do
    id=$(basename $file ".nwk" | cut -f1 -d'_')
    targets=$(basename $file ".nwk" | cut -f2 -d'_')
    alphabet=$(basename $file ".nwk" | cut -f3 -d'_')
    model=$(basename $file ".nwk" | cut -f4 -d'_')
    echo -e "$id\\t$targets\\t$alphabet\\t$model\\t$(cat $file)" >> {output.trees}
done
'''

rule blast_trees:
    input: outdir+"/trees/{seed}_blast_trees.tmp"
    output: outdir+"/trees/{seed}_blast_trees.txt"
    shell: 'mv {input} {output}'


rule get_unrooted_trees:
    input: 
        rules.common_trees.output,
        rules.fs_trees.output,
        rules.blast_trees.output
    output:
        trees=outdir+"/trees/{seed}_unrooted_trees.txt"
    shell:'''
cat {input} > {output}
'''


def seeds_rooted_trees(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_{mode}_{comb}.nwk.rooted", 
                      seed=wildcards.seed, i=seed_genes, mode=modes, comb=combinations)
    return outfiles


checkpoint check_rooted_trees:
    input: seeds_rooted_trees
    output:
        ids=outdir+"/ids/{seed}_finished.ids",
        trees=outdir+"/trees/{seed}_rooted_trees.txt"
    shell:'''
> {output.ids}
> {output.trees}
for file in {input}; do
    id=$(basename $file ".nwk.rooted" | cut -f1 -d'_')
    targets=$(basename $file ".nwk.rooted" | cut -f2 -d'_')
    alphabet=$(basename $file ".nwk.rooted" | cut -f3 -d'_')
    echo $id >> {output.ids}
    echo -e "$id\\t$targets\\t$alphabet\\t$(cat $file)" >> {output.trees}
done
'''


def seeds_treefiles(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_{mode}_{comb}.treeline", 
                      seed=wildcards.seed, i=seed_genes, mode=modes, comb=combinations_ML)
    return outfiles


checkpoint check_ML_trees:
    input: seeds_treefiles
    output: outdir+"/trees/{seed}_mltrees.txt"
    shell:'''
for file in {input}; do
    targets=$(basename $file ".treeline" | cut -f2 -d'_')
    echo -e "$(cat $file)\\t$targets"
done > {output}
'''

def seeds_rangers(wildcards):
    checkpoint_output = checkpoints.check_orphans.get(**wildcards).output.continue_aln
    with open(checkpoint_output) as all_genes:
        seed_genes = [gn.strip() for gn in all_genes]
    outfiles = expand(outdir+"/seeds/{seed}/{i}/{i}_{mode}_{comb}_DTL.txt", 
                      seed=wildcards.seed, i=seed_genes, mode=modes, 
                      comb=combinations)
    return outfiles

rule merge_Ranger:
    input: seeds_rangers
    output: outdir+"/reco/{seed}_DTL.tsv"
    shell: '''
echo -e "gene\\treco_score\\tdups\\ttransfers\\tlosses" > {output}

for file in {input}; do
    gene=$(basename $file | cut -f1 -d'_')
    echo -e "$gene\\g$(cat $file)"
done >> {output}
'''
