outdir=config['outdir']+config['dataset']

rule trim_aln:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    shell: '''
trimal -in {input} -out /dev/stdout -gappyout | seqtk seq -A | \
awk '!/^[X-]+$/' | seqtk seq -L 1 -l 60 > {output}
'''
# trimal -in {input} -out {output} -cons {trimal_cons} -gt {trimal_gt}

rule iqtree_3Di:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    output: 
        tree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_3Di.nwk",
        treeline=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_3Di.treeline",
        ufboot=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_3Di.ufboot"
    wildcard_constraints:
        alphabet="3Di"
    params: 
        submat=config['subst_matrix_tree'],
        ufboot=config['UF_boot']
    log: outdir+"/log/iqtree/{seed}_{i}_{mode}_{alphabet}_3Di.log"
    benchmark: outdir+"/benchmarks/iqtree/{seed}_{i}_{mode}_{alphabet}_3Di.txt"
    threads: 4
    shell: '''
tree_prefix=$(echo {output.tree} | sed 's/.nwk//')

iqtree2 -s {input} --prefix $tree_prefix -B {params.ufboot} -T {threads} --boot-trees --quiet \
--mem 4G --cmin 4 --cmax 10 --mset 3DI -mdef {params.submat}

best_model=$(grep "Model of substitution:" ${{tree_prefix}}.iqtree | cut -f2 -d':' | sed 's/ //')
loglik=$(grep "Log-likelihood of the tree:" ${{tree_prefix}}.iqtree | cut -f2 -d':' | cut -f2 -d' ')

mv ${{tree_prefix}}.treefile {output.tree}

echo -e "{wildcards.i}\\t$best_model\\t$loglik\\t$(cat {output.tree})" > {output.treeline}

mv ${{tree_prefix}}.log {log}
rm -f ${{tree_prefix}}.iqtree ${{tree_prefix}}.model.gz ${{tree_prefix}}.splits.nex ${{tree_prefix}}.contree ${{tree_prefix}}.ckp.gz
'''

rule iqtree_GTR:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    output: 
        tree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_GTR.nwk",
        treeline=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_GTR.treeline",
        ufboot=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_GTR.ufboot"
    wildcard_constraints:
        alphabet="3Di"
    params:
        ufboot=config['UF_boot']
    log: outdir+"/log/iqtree/{seed}_{i}_{mode}_{alphabet}_GTR.log"
    benchmark: outdir+"/benchmarks/iqtree/{seed}_{i}_{mode}_{alphabet}_GTR.txt"
    threads: 4
    shell: '''
tree_prefix=$(echo {output.tree} | sed 's/.nwk//')

iqtree2 -s {input} --prefix $tree_prefix -B {params.ufboot} -T {threads} --boot-trees --quiet \
--mem 4G --cmin 4 --cmax 10 --mset GTR20

best_model=$(grep "Model of substitution:" ${{tree_prefix}}.iqtree | cut -f2 -d':' | sed 's/ //')
loglik=$(grep "Log-likelihood of the tree:" ${{tree_prefix}}.iqtree | cut -f2 -d':' | cut -f2 -d' ')

mv ${{tree_prefix}}.treefile {output.tree}

echo -e "{wildcards.i}\\t$best_model\\t$loglik\\t$(cat {output.tree})" > {output.treeline}

mv ${{tree_prefix}}.log {log}
rm -f ${{tree_prefix}}.iqtree ${{tree_prefix}}.model.gz ${{tree_prefix}}.splits.nex ${{tree_prefix}}.contree ${{tree_prefix}}.ckp.gz
'''


rule iqtree_LG:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    output: 
        tree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_LG.nwk",
        treeline=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_LG.treeline",
        ufboot=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_LG.ufboot"
    wildcard_constraints:
        alphabet="aa"
    params: 
        ufboot=config['UF_boot']
    log: outdir+"/log/iqtree/{seed}_{i}_{mode}_{alphabet}_LG.log"
    benchmark: outdir+"/benchmarks/iqtree/{seed}_{i}_{mode}_{alphabet}_LG.txt"
    threads: 4
    shell: '''
tree_prefix=$(echo {output.tree} | sed 's/.nwk//')

iqtree2 -s {input} --prefix $tree_prefix -B {params.ufboot} -T {threads} --boot-trees --quiet \
--mem 4G --cmin 4 --cmax 10 --mset LG

best_model=$(grep "Model of substitution:" ${{tree_prefix}}.iqtree | cut -f2 -d':' | sed 's/ //')
loglik=$(grep "Log-likelihood of the tree:" ${{tree_prefix}}.iqtree | cut -f2 -d':' | cut -f2 -d' ')

mv ${{tree_prefix}}.treefile {output.tree}

echo -e "{wildcards.i}\\t$best_model\\t$loglik\\t$(cat {output.tree})" > {output.treeline}

mv ${{tree_prefix}}.log {log}
rm -f ${{tree_prefix}}.iqtree ${{tree_prefix}}.model.gz ${{tree_prefix}}.splits.nex ${{tree_prefix}}.contree ${{tree_prefix}}.ckp.gz
'''

##### FOLDTREE #####

rule foldseek_allvall_tree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}.ids"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_allvall.txt"
    params: config['structure_dir']
    log: outdir+"/log/foldseek/{seed}_{i}_{mode}.log"
    benchmark: outdir+"/benchmarks/foldseek/{seed}_{i}_{mode}.txt"
    shell:'''
structdir=$(dirname {input})/structs_{wildcards.mode}
mkdir -p $structdir

for id in $(cat {input}); do
zcat {params}*/high_cif/${{id}}-model_v4.cif.gz > $structdir/${{id}}.cif
done

foldseek easy-search $structdir $structdir {output} $TMPDIR/{wildcards.i} \
--format-output 'query,target,fident,lddt,alntmscore' --exhaustive-search -e inf > {log}

rm -r $structdir
'''

rule foldseek_distmat:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_allvall.txt"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di_fident.txt"
        # outdir+"/seeds/{seed}/{i}/{i}_alntmscore.txt",
        # outdir+"/seeds/{seed}/{i}/{i}_lddt.txt"
    script: "../../software/foldtree/foldseekres2distmat_simple.py"
#     shell:'''
# python ./scripts/foldseek2distmat.py -i {input} -o {output}
# '''

rule foldtree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_fident.txt"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_FT.nwk"
    wildcard_constraints:
        alphabet="3Di"
    benchmark: outdir+"/benchmarks/foldtree/{seed}_{i}_{mode}_{alphabet}_FT.txt"
    shell:'''
quicktree -i m {input} | paste -s -d '' > {output}
'''

rule quicktree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_QT.nwk"
    wildcard_constraints:
        alphabet="aa"
    params: config["distboot"]
    log: outdir+"/log/quicktree/{seed}_{i}_{mode}_{alphabet}_QT.log"
    benchmark: outdir+"/benchmarks/quicktree/{seed}_{i}_{mode}_{alphabet}_QT.txt"
    shell:'''
esl-reformat stockholm {input} | quicktree -boot {params} -in a -out t /dev/stdin | paste -s -d '' > {output} 2> {log}
'''

rule RangerDTL:
    input:
        sptree=config['species_tree'],
        genetree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.nwk",
        taxidmap=rules.make_taxidmap_sp.output
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}_DTL.txt"
    shell:'''
echo -e "$(cat {input.sptree})\\n$(nw_rename {input.genetree} {input.taxidmap})" | \
Ranger-DTL.linux -i /dev/stdin  -o /dev/stdout -T 2000 -q -s | grep reco | \
grep -E "[0-9]+" -o | paste -s -d'\\t' > {output}
'''

rule root_tree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{algorithm}_{model}.nwk"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{algorithm}_{model}.nwk.rooted"
    log: outdir+"/log/mad/{seed}_{i}_{mode}_{algorithm}_{model}.log"
    shell: '''
mad {input} > {log}
sed -i \'2,$d\' {output}
'''

# rule foldtree:
#     input: outdir+"/seeds/{seed}/{i}/{i}_{mode}.ids"
#     output:
#         distmat=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_FT.txt",
#         tree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_FT.nwk"
#     wildcard_constraints:
#         alphabet="3Di"
#     params: config['structure_dir']
#     log: outdir+"/log/ft/{seed}_{i}_{mode}_{alphabet}.log"
#     benchmark: outdir+"/benchmarks/ft/{seed}_{i}_{mode}_{alphabet}.txt"
#     shell: '''
# indir=$(dirname {input})
# mkdir -p $indir/structs_{wildcards.mode}

# for id in $(cat {input}); do
# zcat {params}/*/high_cif/${{id}}-model_v4.cif.gz > $indir/structs_{wildcards.mode}/${{id}}.cif
# done

# python ./software/foldtree/foldtree.py -i $indir/structs_{wildcards.mode} -o $indir/{wildcards.i}_{wildcards.mode} \
# -t $TMPDIR/{wildcards.i}_{wildcards.mode}_ft --outtree {output.tree} -c $indir/{wildcards.i}_{wildcards.mode}_core \
# --corecut --correction --kernel fident > {log}

# rm -r $indir/structs_{wildcards.mode}
# rm -r $indir/{wildcards.i}_{wildcards.mode}_core
# rm {output.distmat}_fastme_stat.txt {output.distmat}.tmp
# rm $indir/{wildcards.i}_{wildcards.mode}_allvall.tsv
# rm $indir/{wildcards.i}_{wildcards.mode}_core_allvall.tsv
# '''
# rm -r $TMPDIR/{wildcards.i}_{wildcards.mode}_ft
