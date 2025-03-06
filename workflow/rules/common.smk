rule get_seqs_aa:
    input: 
        ids=outdir+"/seeds/{seed}/{i}/{i}_{mode}.ids",
        fa=homodir+"/db/all_seqs.fa"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_aa.seqs",
    conda: "../envs/sp_utils.yaml"
    shell:'''
seqkit grep -f {input.ids} {input.fa} > {output}
'''

rule aln_aa:
    input: rules.get_seqs_aa.output
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_aa.alg"
    log: outdir+"/log/mafft/{seed}_{i}_{mode}_aa.log"
    benchmark: outdir+"/benchmarks/mafft/{seed}_{i}_{mode}_aa.txt"
    threads: 4
    conda: "../envs/sp_tree.yaml"
    shell:'''
mafft --auto --thread {threads} {input} > {output} 2> {log}
'''

# If foldseek, retrieve sequences from "translated version"
rule get_seqs_3Di:
    input:
        ids=outdir+"/seeds/{seed}/{i}/{i}_{mode}.ids",
        fa=homodir+"/db/all_seqs_fsdb_ss.fa"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.seqs"
    conda: "../envs/sp_utils.yaml"
    shell:'''
seqkit grep -f {input.ids} {input.fa} > {output}
'''

rule mask_seqs_3Di:
    input: rules.get_seqs_3Di.output
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.masked.seqs"
    params: 
        structdir=config['data_dir']+"structures/",
        min_lddt=config['min_lddt']
    conda: "../envs/sp_python.yaml"
    script: "../scripts/mask_structures.py"

rule struct_tmp_db:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}.ids"
    output: temp(directory(outdir+"/seeds/{seed}/{i}/structs_{i}_{mode}"))
    params: config['data_dir']+"structures/"
    shell:'''
mkdir -p {output}

for id in $(cat {input}); do
    if [ -f {params}*/high_cif/${{id}}-model_v4.cif.gz ]; then
        zcat {params}*/high_cif/${{id}}-model_v4.cif.gz > {output}/${{id}}.cif
    fi
done
'''

rule foldmason:
    input: 
        db=rules.struct_tmp_db.output
    output: 
        fa=outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.alg"
    log: outdir+"/log/foldmason/{seed}_{i}_{mode}_3Di.log"
    benchmark: outdir+"/benchmarks/foldmason/{seed}_{i}_{mode}_3Di.txt"
    threads: 4
    conda: "../envs/sp_tree.yaml"
    shell: '''
output={output}.tmp
foldmason easy-msa --threads {threads} {input} $output $TMPDIR/{wildcards.i}_{wildcards.mode}_fm > {log}
rm ${{output}}.nw ${{output}}_aa.fa
mv ${{output}}_3di.fa {output}
sed -i 's/.cif//g' {output}
'''
# you need to do a tmp file because foldmason automaticlly computes an aa file.

# rule aln_3Di:
#     input: rules.mask_seqs_3Di.output
#     output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.alg"
#     params: 
#         submat=config['subst_matrix'][']
#     log: outdir+"/log/mafft/{seed}_{i}_{mode}_3Di.log"
#     benchmark: outdir+"/benchmarks/mafft/{seed}_{i}_{mode}_3Di.txt"
#     threads: 4
#     conda: "../envs/sp_tree.yaml"
#     shell:'''
# mafft --auto --thread {threads} --aamatrix {params.submat} {input} > {output} 2> {log}
# '''

rule trim_aln:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    conda: "../envs/sp_tree.yaml"
    shell: '''
trimal -in {input} -out /dev/stdout -gappyout | seqtk seq -A | \
awk '!/^[X-]+$/' | seqtk seq -L 1 -l 60 > {output}
'''
# trimal -in {input} -out {output} -cons {trimal_cons} -gt {trimal_gt}

rule concat_aln:
    input:
        aa=outdir+"/seeds/{seed}/{i}/{i}_{mode}_aa.alg.clean",
        msa_3di=outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.alg.clean"
    output:
        outdir+"/seeds/{seed}/{i}/{i}_{mode}_comb.alg.clean"
    shell:"""
seqkit concat {input.aa} {input.msa_3di} > {output}
"""

rule iqtree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    output: 
        tree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.nwk",
        # treeline=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.treeline",
        # ufboot=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.ufboot",
        iqtree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.iqtree"
    wildcard_constraints:
        model="GTR|LG|3Di|LLM|AF"
    params: 
        threedi_submat=config['subst_matrix']['3di'],
        LLM_submat=config['subst_matrix']['LLM'],
        AF_submat=config['subst_matrix']['AF'],
        ufboot=config['UF_boot']
    log: outdir+"/log/iqtree/{seed}_{i}_{mode}_{alphabet}_{model}.log"
    benchmark: outdir+"/benchmarks/iqtree/{seed}_{i}_{mode}_{alphabet}_{model}.txt"
    threads: 4
    conda: "../envs/sp_tree.yaml"
    shell: '''
tree_prefix=$(echo {output.tree} | sed 's/.nwk//')
model={wildcards.model}

if [ "$model" == "GTR" ]; then
    model="GTR20"
elif [ "$model" == "3Di" ]; then
    model="3DI -mdef {params.threedi_submat}"
elif [ "$model" == "LLM" ]; then
    model={params.LLM_submat}
elif [ "$model" == "AF" ]; then
    model={params.AF_submat}
fi

iqtree2 -s {input} --prefix $tree_prefix -B {params.ufboot} -T {threads}  --quiet \
--mem 4G --cmin 4 --cmax 10 --mset $model

mv ${{tree_prefix}}.treefile {output.tree}
mv ${{tree_prefix}}.log {log}
rm -f ${{tree_prefix}}.model.gz ${{tree_prefix}}.splits.nex ${{tree_prefix}}.contree ${{tree_prefix}}.ckp.gz
'''
# --boot-trees

rule get_part: 
    input: 
        LG=outdir+"/seeds/{seed}/{i}/{i}_{mode}_aa_LG.iqtree",
        ThreeDI=outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di_3Di.iqtree"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_comb.alg.clean.part"
    shell: """
model_LG=$(grep "^Model of substitution" {input.LG} | cut -f2 -d ':')
model_3Di=$(grep "^Model of substitution" {input.ThreeDI} | cut -f2 -d ':')
length_LG=$(grep "^Input data" {input.LG} | cut -f6 -d ' ')
length_3Di=$(grep "^Input data" {input.ThreeDI} | cut -f6 -d ' ')

end=$(expr $length_LG + $length_3Di)
start_3Di=$(expr $length_LG + 1)

echo "$model_LG, {wildcards.i}_{wildcards.mode}_aa = 1-$length_LG" > {output}
echo "$model_3Di, {wildcards.i}_{wildcards.mode}_3Di = $start_3Di-$end\n" >> {output}
"""
# should look like this: 
# LG+R4, p1_COG2207_common_3Di = 1-220
# 3DI+I+G4, p2_COG2207_common_aa = 221-442

rule iqtree_partitioned:
    input: 
        fa=rules.concat_aln.output,
        part=rules.get_part.output
    output: 
        tree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_comb_part.nwk",
        iqtree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_comb_part.iqtree"
    wildcard_constraints:
        model="GTR|LG|3Di|LLM|AF"
    params: 
        threedi_submat=config['subst_matrix']['3di'],
        LLM_submat=config['subst_matrix']['LLM'],
        AF_submat=config['subst_matrix']['AF'],
        ufboot=config['UF_boot']
    log: outdir+"/log/iqtree/{seed}_{i}_{mode}_comb_part.log"
    benchmark: outdir+"/benchmarks/iqtree/{seed}_{i}_{mode}_comb_part.txt"
    threads: 4
    conda: "../envs/sp_tree.yaml"
    shell: """
tree_prefix=$(echo {output.tree} | sed 's/.nwk//')

iqtree2 -s {input.fa} -p {input.part} --prefix $tree_prefix -mdef {params.threedi_submat} -B {params.ufboot} -T {threads} --quiet

mv ${{tree_prefix}}.treefile {output.tree}
mv ${{tree_prefix}}.log {log}
rm -f ${{tree_prefix}}.model.gz ${{tree_prefix}}.splits.nex ${{tree_prefix}}.contree ${{tree_prefix}}.ckp.gz
"""


##### FOLDTREE #####

rule foldmason_report:
    input: 
        fa=outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.alg",
        db=homodir+"/db/all_seqs_fsdb"
    output: 
        fa=temp(outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.alg.tmp"),
        html=outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.html"
    conda: "../envs/sp_tree.yaml"
    shell: '''
cat <(seqkit grep -p "AF-{wildcards.i}-F1" {input.fa}) <(seqkit grep -v -p "AF-{wildcards.i}-F1" {input.fa}) |\
seqkit replace -p $ -r -model_v4.cif > {output.fa}
foldmason msa2lddtreport {input.db} {output.fa} {output.html}
'''

rule foldseek_allvall_tree:
    input: rules.struct_tmp_db.output
    output: temp(outdir+"/seeds/{seed}/{i}/{i}_{mode}_allvall.txt")
    log: outdir+"/log/foldseek/{seed}_{i}_{mode}.log"
    benchmark: outdir+"/benchmarks/foldseek/{seed}_{i}_{mode}.txt"
    conda: "../envs/sp_homology.yaml"
    shell:'''
foldseek easy-search {input} {input} {output} $TMPDIR/{wildcards.i} \
--format-output 'query,target,fident,lddt,alntmscore' --exhaustive-search -e inf \
--alignment-type 2 > {log}
'''

rule foldseek_distmat:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_allvall.txt"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di_fident.txt"
    conda: "../envs/sp_python.yaml"
    script: "../scripts/foldtree/foldseekres2distmat_simple.py"

rule foldtree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_fident.txt"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_FT.nwk"
    wildcard_constraints:
        alphabet="3Di"
    benchmark: outdir+"/benchmarks/foldtree/{seed}_{i}_{mode}_{alphabet}_FT.txt"
    conda: "../envs/sp_tree.yaml"
    shell:'''
quicktree -i m {input} | paste -s -d '' > {output}
'''

rule foldtree_py:
    input: rules.struct_tmp_db.output
    output:
        distmat=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_FTPY.txt",
        tree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_FTPY.nwk"
    wildcard_constraints:
        alphabet="3Di"
    params: config['data_dir']+"structures/"
    log: outdir+"/log/ft/{seed}_{i}_{mode}_{alphabet}_FTPY.log"
    benchmark: outdir+"/benchmarks/ft/{seed}_{i}_{mode}_{alphabet}_FTPY.txt"
    conda: "../envs/sp_python.yaml"
    shell: '''
indir=$(dirname {input})

python workflow/scripts/foldtree/foldtree.py -i {input} -o $indir/{wildcards.i}_{wildcards.mode} \
-t $TMPDIR/{wildcards.i}_{wildcards.mode}_ft --outtree {output.tree} -c $indir/{wildcards.i}_{wildcards.mode}_core \
--corecut --correction --kernel fident > {log}

rm -r $indir/structs_{wildcards.mode}
rm -r $indir/{wildcards.i}_{wildcards.mode}_core
rm {output.distmat}_fastme_stat.txt {output.distmat}.tmp
rm $indir/{wildcards.i}_{wildcards.mode}_allvall.tsv
rm $indir/{wildcards.i}_{wildcards.mode}_core_allvall.tsv
'''
# rm -r $TMPDIR/{wildcards.i}_{wildcards.mode}_ft

rule convert_phylip:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_aa.alg.clean"
    output: temp(outdir+"/seeds/{seed}/{i}/{i}_{mode}_aa.alg.clean.phy")
    conda: "../envs/sp_python.yaml"
    script: "../scripts/fasta2phylip.py"

rule fastme:
    input: rules.convert_phylip.output
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_aa_FM.nwk"
    params: config["distboot"]
    log: outdir+"/log/fastme/{seed}_{i}_{mode}_aa_FM.log"
    benchmark: outdir+"/benchmarks/fastme/{seed}_{i}_{mode}_aa_FM.txt"
    threads: 4
    conda: "../envs/sp_tree.yaml"
    shell:'''
fastme -q -p -T {threads} -b {params} -i {input} -o {output} > {log}
rm {input}_fastme*
'''


# rule quicktree:
#     input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
#     output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_QT.nwk"
#     wildcard_constraints:
#         alphabet="aa"
#     params: config["distboot"]
#     log: outdir+"/log/quicktree/{seed}_{i}_{mode}_{alphabet}_QT.log"
#     benchmark: outdir+"/benchmarks/quicktree/{seed}_{i}_{mode}_{alphabet}_QT.txt"
#     conda: "../envs/sp_tree.yaml"
#     shell:'''
# esl-reformat stockholm {input} | quicktree -boot {params} -in a -out t /dev/stdin | paste -s -d '' > {output} 2> {log}
# '''


# rule RangerDTL:
#     input:
#         sptree=config['species_tree'],
#         genetree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.nwk",
#         taxidmap=rules.make_taxidmap_sp.output
#     output: 
#         full=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}_ranger.txt"
#     conda: "../envs/sp_tree.yaml"
#     shell:'''
# echo -e "$(cat {input.sptree})\\n$(nw_rename {input.genetree} {input.taxidmap})" | \
# Ranger-DTL.linux -i /dev/stdin -o {output.full} -T 2000 -q
# '''


rule Notung:
    input: 
        genetree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.nwk",
        sptree=config['species_tree'],
        taxidmap=homodir+"/db/taxidmap_sps"
    output: 
        map_ids=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.map",
        genetree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}_rnm.nwk",
        reco=outdir+"/reco/notung/{seed}/{i}_{mode}_{alphabet}_{model}_rnm.nwk.rooting.0.parsable.txt"
    log: outdir+"/log/notung/{seed}_{i}_{mode}_{alphabet}_{model}_reco.log"
    conda: "../envs/reco.yaml"
    shell: '''
nw_labels {input.genetree} -I | \
csvtk join -H -t -f 1 {input.taxidmap} - | \
awk '{{print $1"\\t"$1"_"$2}}' > {output.map_ids}

nw_rename {input.genetree} {output.map_ids} > {output.genetree}

notung --root --maxtrees 1 -g {output.genetree} -s {input.sptree} \
--speciestag postfix --parsable --outputdir $(dirname {output.reco}) > {log}
'''
