import pandas as pd
import os

def get_valid_subdirs(parent_dir):
    valid_subdirs = []
    
    # List all items in the directory
    for subdir in os.listdir(parent_dir):
        subdir_path = os.path.join(parent_dir, subdir)
        
        # Check if it's a directory and iqtree has been computed
        if os.path.isdir(subdir_path) and os.path.exists(os.path.join(subdir_path, 'treescores_sequences_iq.muscle.json')):
            valid_subdirs.append(subdir)
    
    return valid_subdirs

cwd = os.getcwd()  # This should be inside snakemake workflow but I can't find it

folders = get_valid_subdirs(cwd)
print(folders)

if len(folders) == 0:
    exit("No folder has trees")

# Load fixed config parameters
configfile: workflow.basedir+'/../config/params_ortho_benchmark.yaml'

mattypes = ['fident', 'alntmscore', 'lddt' ]
metrics = ["apro", "disco", "notung"]

rule all:
    input:
        "reconciliation/TCS.txt",
        expand("reconciliation/{metric}_{mat}.txt", metric=metrics, mat = mattypes + ['iq'])
        # expand('reconciliation/plots/comparison.pdf', db=config['dbs'])


rule get_map:
    input: expand("{oma}/sequence_dataset.csv", oma=folders)
    output: 
        ids="reconciliation/ids.txt",
        sps="reconciliation/sps.txt"
    shell: """
cat {input} | awk -F',' '{{print $3"\\t"$2}}' | grep -v Entry > {output.ids}
awk '{{split($2, arr, "_"); print $2" "arr[2]}}' {output.ids} > {output.sps}
"""


rule rename_iq:
    input:
        tree="{oma}/sequences.aln.muscle.fst.treefile",
        ids=rules.get_map.output.ids
    output: "reconciliation/trees/tmp/{oma}_iq.nwk"
    shell: """
gotree rename -i {input.tree} -m {input.ids} -o {output}
"""


rule rename_ft:
    input:
        tree="{oma}/{mat}_1_exp_struct_tree.PP.nwk",
        ids=rules.get_map.output.ids
    output: "reconciliation/trees/tmp/{oma}_{mat}.nwk"
    shell: """
gotree rename -i {input.tree} -m {input.ids} -o {output}
"""

ruleorder: rename_iq > rename_ft > concat

rule concat:
    input: 
        trees=expand("reconciliation/trees/tmp/{oma}_{{mat}}.nwk", oma=folders),
        sps=rules.get_map.output.sps
    output: 
        ids="reconciliation/trees/{mat}.nwk",
        sps="reconciliation/trees/{mat}_sps.nwk"
    conda: workflow.basedir+'/envs/reco.yaml'
    shell: """
cat {input.trees} > {output.ids} 
nw_rename {input.trees} {input.sps} > {output.sps} 
"""


rule subset_sptree:
    input:
        sptree=config["sptree"],
        sps=rules.get_map.output.sps
    output: "reconciliation/sptree.nwk"
    conda: workflow.basedir+'/envs/reco.yaml'
    shell: """
nw_prune -v {input.sptree} $(cut -f2 {input.sps} | sort -u | tr '\\n' ' ') > {output}
"""


rule astral_pro:
    input:
        sptree=rules.subset_sptree.output,
        sps=rules.get_map.output.sps,
        gt=rules.concat.output.ids
    output: "reconciliation/apro_{mat}.txt"
    conda: workflow.basedir+'/envs/reco.yaml'
    shell: """
astral-pro -c {input.sptree} -a {input.sps} -u 2 -i {input.gt} -o {output} -C
"""


rule disco:
    input:
        gt=rules.concat.output.sps
    output: "reconciliation/disco_{mat}.txt"
    conda: workflow.basedir+'/envs/reco.yaml'
    shell: """
disco.py -i {input.gt} -o {output}
"""


rule notung:
    input:
        sptree=rules.subset_sptree.output,
        gt="reconciliation/trees/tmp/{oma}_{mat}.nwk"
    output: "reconciliation/notung/single/{oma}_{mat}.nwk.rooting.0.parsable.txt"
    conda: workflow.basedir+'/envs/reco.yaml'
    shell: """
notung --root --maxtrees 1 -g {input.gt} -s {input.sptree} --speciestag postfix --parsable --outputdir $(dirname {output})
"""


rule merge_notung:
    input: expand("reconciliation/notung/single/{oma}_{{mat}}.nwk.rooting.0.parsable.txt", oma=folders)
    output: "reconciliation/notung_{mat}.txt"
    shell: """
for file in {input}; do 
    echo -e "$(basename $file | cut -f1 -d'_')\\t$(awk 'NR==1' $file | cut -f2,5)"
done > {output}
"""


rule TCS_score:
    input: 
        expand("{oma}/{mat}_1_exp_treescores_struct_tree.json", oma=folders, mat=mattypes),
        expand("{oma}/treescores_sequences_iq.muscle.json", oma=folders)
    output: "reconciliation/TCS.txt"
    shell: """
cat {input} | jq -r 'to_entries[] | "\(.key | split("/")[0])	\
    \(.key | split("/")[1] | split("_")[0])\\t\(.value.score)"' | \
    sed 's/sequences.aln.muscle.fst.treefile.rooted/iqtree/' > {output}    
"""
