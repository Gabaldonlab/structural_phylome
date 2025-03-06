import pandas as pd

# Load fixed config parameters
configfile: 'config/params_ortho_benchmark.yaml'

input_table = pd.read_csv(config['uniprot_df'], sep='\t')
input_table.columns = ['uniprot', 'taxid', 'mnemo']
input_dict = input_table.set_index('uniprot').T.to_dict()

codes = list(input_table['uniprot'])
# aln_type = [2]

rule all:
    input:
        expand('results/ortho_bench/{db}/accuracy.pdf', db=config['dbs'])


rule get_ids:
    input:
        idmap = 'results/homology/hsap_euka/db/taxidmap',
        meta = 'results/data/meta/hsap_euka_uniprot_genomes.tsv',
        ids=expand('results/data/ids/{code}_ids.tsv', code=codes)
    output: 'results/ortho_bench/{db}/db/ids.tsv'
    params: min_seqs = config["min_seqs"]
    conda: 'envs/sp_R.yaml'
    script: 'scripts/get_ids.R'


rule get_fasta:
    input: 
        ids=rules.get_ids.output,
        fa=expand('results/data/fastas/{code}.fa', code=codes)
    output: 'results/ortho_bench/{db}/db/db.fa'
    conda: 'envs/sp_utils.yaml'
    shell: '''
seqkit grep -f {input.ids} {input.fa} > {output}    
'''

rule get_tsv:
    input: 
        ids=rules.get_ids.output,
        folder='results/data/structures/'
    output: 'results/ortho_bench/{db}/db/struct_files.tsv'
    shell: '''
find {input.folder}/*/high_cif -name '*cif.gz' | grep -f {input.ids} > {output}
'''

rule make_foldseekdb:
    input:
        tsv=rules.get_tsv.output
    output: 'results/ortho_bench/{db}/db/all_seqs_fsdb'
    conda: 'envs/sp_homology.yaml'
    shell:'''
foldseek createdb {input.tsv} {output}
'''

rule make_blastdb:
    input:
        fa=rules.get_fasta.output
    output: 'results/ortho_bench/{db}/db/all_seqs_blastdb.pdb'
    conda: 'envs/sp_homology.yaml'
    shell:'''
out_file=$(echo {output} | rev | cut -f 2- -d '.' | rev)
makeblastdb -in {input.fa} -out $out_file -dbtype prot
'''

rule blast:
    input:
        q=rules.get_fasta.output,
        db=rules.make_blastdb.output,
    output: 'results/ortho_bench/{db}/hits_blast.tsv'
    params: config['target_seqs']
    benchmark: 'results/ortho_bench/{db}/benchmarks/blast.txt'
    threads: 28
    conda: 'envs/sp_homology.yaml'
    shell:'''
dbfile=$(echo {input.db} | rev | cut -f 2- -d '.' | rev)
blastp -query {input.q} -db $dbfile -out {output} -max_hsps 1 \
-outfmt '6 std qcovs qcovhsp qlen slen' -num_threads {threads}  -max_target_seqs {params}
'''


rule foldseek:
    input:
        db=rules.make_foldseekdb.output
    output: 'results/ortho_bench/{db}/hits_fs.tsv'
    params: config['target_seqs']
    benchmark: 'results/ortho_bench/{db}/benchmarks/fs.txt'
    threads: 28
    conda: 'envs/sp_homology.yaml'
    shell:'''
foldseek easy-search {input.db} {input.db} {output} $TMPDIR/{wildcards.db}_fs --threads {threads} --max-seqs {params} \
--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,rmsd,prob,qcov,tcov 

sed -i 's/-model_v4.cif//g' {output} 
'''
# --alignment-type {wildcards.aln}

rule plot_results:
    input:
        fa=rules.get_fasta.output,
        blast=rules.blast.output,
        fs=rules.foldseek.output,
        ids=expand('results/data/ids/{code}_ids.tsv', code=codes)
    params: 
        min_seqs=config['min_seqs'],
        evalue=config['evalue']
    output: 'results/ortho_bench/{db}/accuracy.pdf'
    conda: 'envs/sp_R.yaml'
    script: 'scripts/analyze_orthobench.R'