rule make_blastdb:
    input:
        fa=expand(config['data_dir']+'fastas/{code}.fa', code=codes),
        maps=expand(outdir+'/db/{code}.taxid', code=codes)
    output:
        fa=outdir+"/db/all_seqs.fa",
        blast=outdir+"/db/all_seqs_blastdb.pdb",
        mapid=outdir+"/db/taxidmap"
    conda: "../envs/sp_homology.yaml"
    shell:'''
cat {input.maps} > {output.mapid}
cat {input.fa} > {output.fa}
out_file=$(echo {output.blast} | rev | cut -f 2- -d '.' | rev)
makeblastdb -in {output.fa} -out $out_file -parse_seqids -taxid_map {output.mapid} -dbtype prot
'''

rule make_taxidmap_sp:
    input:
        taxid=rules.make_blastdb.output.mapid,
        table=config['data_dir']+'meta/'+config["homology_dataset"]+'_uniprot_genomes.tsv'
    output: outdir+"/db/taxidmap_sps"
    conda: "../envs/sp_utils.yaml"
    shell:'''
awk 'NR>1' {input.table} | cut -f2,3 | csvtk join -H -t -f"2;1" {input.taxid} - | \
awk '{{print $1"\\t"$3}}' | shuf > {output}
'''
# I shuf because in very few instances nw_rename fails!!! Magic stuff

# rule make_taxidmap_ale:
#     input: 
#         taxid=rules.make_blastdb.output.mapid,
#         groups=config['taxons_file']
#     output: outdir+"/db/taxidmap_ale"
#     shell:'''
# awk 'NR>1' {input.groups} | cut -f2,11 | csvtk join -H -t -f"2;1" {input.taxid} - | \
# awk '{{split($1,a,"-"); print $1"\\t"$3"_"a[2]}}' > {output}
# '''

rule make_blastdb_single:
    input:
        fa=rules.make_fastas.output,
        mapid=rules.make_blastdb.output.mapid
    output:
        outdir+"/db/single_dbs/{code}.pdb"
    conda: "../envs/sp_homology.yaml"
    shell:'''
out_file=$(echo {output} | rev | cut -f 2- -d '.' | rev)
makeblastdb -in {input.fa} -out $out_file -parse_seqids -taxid_map {input.mapid} -dbtype prot
'''

rule blast:
    input:
        q=config['data_dir']+'fastas/{seed}.fa',
        db=rules.make_blastdb.output.blast,
    output: outdir+"/{seed}_blast.tsv"
    params: config['target_seqs']
    benchmark: outdir+"/benchmarks/homology/{seed}_blast.txt"
    threads: 24
    conda: "../envs/sp_homology.yaml"
    shell:'''
dbfile=$(echo {input.db} | rev | cut -f 2- -d '.' | rev)
blastp -query {input.q} -db $dbfile -out {output} -max_hsps 1 -max_target_seqs {params} \
-outfmt "6 std qcovs qcovhsp qlen slen staxids" -num_threads {threads}
'''

rule filter_blast:
    input: rules.blast.output
    output: outdir+"/{seed}_blast_filtered.tsv"
    params: 
        eval_both=config['eval_both'],
        coverage=config['coverage'],
        max_seqs=config['max_seqs']
    shell: """
awk '$11<{params.eval_both} && $4/$15*100>{params.coverage} && $4/$16*100>{params.coverage}' {input} > {output}    
"""

rule blast_allvall:
    input:
        q=config['data_dir']+'fastas/{seed}.fa',
        t=config['data_dir']+'fastas/{code}.fa',
        q_db=outdir+"/db/single_dbs/{seed}.pdb",
        t_db=outdir+"/db/single_dbs/{code}.pdb",
    output: 
        q_t=outdir+"/allvall/{seed}_{code}_blast.tsv",
        t_q=outdir+"/allvall/{code}_{seed}_blast.tsv"
    params: config['max_seqs_brh']
    benchmark: outdir+"/benchmarks/homology/{code}_{seed}_blast.txt"
    threads: 4
    conda: "../envs/sp_homology.yaml"
    shell:'''
dbfile1=$(echo {input.t_db} | rev | cut -f 2- -d '.' | rev)
blastp -query {input.q} -db $dbfile1 -out {output.q_t} \
-outfmt "6 std qcovs qcovhsp qlen slen staxids" -num_threads {threads} -max_hsps 1 -max_target_seqs {params}

dbfile2=$(echo {input.q_db} | rev | cut -f 2- -d '.' | rev)
blastp -query {input.t} -db $dbfile2 -out {output.t_q} \
-outfmt "6 std qcovs qcovhsp qlen slen staxids" -num_threads {threads} -max_hsps 1 -max_target_seqs {params}
'''

rule blast_brh:
    input: 
        a=expand(outdir+"/allvall/{seed}_{code}_blast.tsv", seed=config['seed'], code=codes),
        b=expand(outdir+"/allvall/{code}_{seed}_blast.tsv", seed=config['seed'], code=codes)
    output: outdir+"/{seed}_blast_brh.tsv"
    params: config['eval_brh']
    conda: "../envs/sp_python.yaml"
    script: "../scripts/get_BRH.py"


rule gene_map:
    input:
        table=config['data_dir']+'meta/'+config["homology_dataset"]+'_uniprot_genomes.tsv',
        taxidmap=rules.make_blastdb.output.mapid
    output:
        outdir+"/db/gene_species.map"
    conda: "../envs/sp_utils.yaml"
    shell:'''
csvtk join -H -t -f 2 -L {input.taxidmap} {input.table} | cut -f1,4 > {output}
'''