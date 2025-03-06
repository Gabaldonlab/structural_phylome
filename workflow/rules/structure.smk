rule make_foldseekdb:
    input:
        all_struct=expand(config['data_dir']+'structures/{code}/high_cif', code=codes),
    output:
        outdir+"/db/all_seqs_fsdb"
    conda: "../envs/sp_homology.yaml"
    shell:'''
structsdir=$(dirname {output})/structs
mkdir -p $structsdir

for id in {input}; do
bn=$(dirname $id | rev | cut -f1 -d'/' | rev)

if [ -L $structsdir/$bn ]; then
    unlink $structsdir/$bn
fi
ln -s -f -r $id $structsdir/$bn
done

foldseek createdb $structsdir {output}
'''

rule make_foldseekdb_seq:
    input: rules.make_foldseekdb.output
    output: outdir+"/db/all_seqs_fsdb_ss.fa"
    conda: "../envs/sp_homology.yaml"
    shell:'''
foldseek lndb {input}_h {input}_ss_h
foldseek convert2fasta {input}_ss {output}
sed -i 's/-model_v4.cif//g' {output}
'''

rule make_foldseekdb_single:
    input: config['data_dir']+'structures/{code}/high_cif',
    output: outdir+"/db/single_dbs/{code}_fsdb"
    threads: 8
    conda: "../envs/sp_homology.yaml"
    shell:'''
foldseek createdb --threads {threads} {input} {output}
'''

rule foldseek:
    input:
        q=config['data_dir']+'structures/{seed}/high_cif',
        db=rules.make_foldseekdb.output
    output: 
        outdir+"/{seed}_fs.tsv"
    params: config['target_seqs']
    log: outdir+"/log/homology/{seed}_fs.log"
    benchmark: outdir+"/benchmarks/homology/{seed}_fs.txt"
    threads: 24
    conda: "../envs/sp_homology.yaml"
    shell:'''
foldseek easy-search {input.q} {input.db} {output} $TMPDIR/{wildcards.seed} --threads {threads} --max-seqs {params} \
--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,rmsd,prob,qcov,tcov > {log}
sed -i 's/-model_v4.cif//g' {output}
'''
# foldseek results are sorted by bitScore * sqrt(alnlddt * alntmscore) https://github.com/steineggerlab/foldseek/issues/85

rule filter_foldseek:
    input: rules.foldseek.output
    output: outdir+"/{seed}_fs_filtered.tsv"
    params: 
        eval_both=config['eval_both'],
        coverage=config['coverage'],
        max_seqs=config['max_seqs']
    shell: """
awk '$11<{params.eval_both} && $17*100>{params.coverage} && $18*100>{params.coverage}' {input} > {output}
"""

rule foldseek_allvall:
    input:
        q=config['data_dir']+'structures/{seed}/high_cif',
        t=config['data_dir']+'structures/{code}/high_cif'
    output:
        q_t=outdir+"/allvall/{seed}_{code}_fs.tsv",
        t_q=outdir+"/allvall/{code}_{seed}_fs.tsv"
    params: config['max_seqs_brh']
    log: outdir+"/log/homology/{code}_{seed}_fs.log"
    benchmark: outdir+"/benchmarks/homology/{code}_{seed}_fs.txt"
    threads: 4
    conda: "../envs/sp_homology.yaml"
    shell:'''
foldseek easy-search {input.q} {input.t} {output.q_t} $TMPDIR/{wildcards.seed}_{wildcards.code} \
--threads {threads} --max-seqs {params} \
--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,rmsd,prob,qcov,tcov > {log}
sed -i 's/-model_v4.cif//g' {output.q_t}

foldseek easy-search {input.t} {input.q} {output.t_q} $TMPDIR/{wildcards.code}_{wildcards.seed} \
--threads {threads} --max-seqs {params} \
--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,rmsd,prob,qcov,tcov > {log}
sed -i 's/-model_v4.cif//g' {output.t_q}
'''

rule foldseek_brh:
    input: 
        a=expand(outdir+"/allvall/{seed}_{code}_fs.tsv", seed=config['seed'], code=codes),
        b=expand(outdir+"/allvall/{code}_{seed}_fs.tsv", seed=config['seed'], code=codes)
    output:
        outdir+"/{seed}_fs_brh.tsv"
    params: config['eval_brh']
    conda: "../envs/sp_python.yaml"
    script: "../scripts/get_BRH.py"


