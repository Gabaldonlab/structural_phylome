outdir=config['outdir']+'homology/'+config['homology_dataset']

input_table = pd.read_csv(config['taxids'], sep='\t')
input_table.columns = ['uniprot', 'taxid', 'mnemo']
input_dict = input_table.set_index('uniprot').T.to_dict()

codes = list(input_table['uniprot'])


rule get_nums:
    input: 
        seqs=expand(config['data_dir']+'fastas/{code}.fa', code=codes),
        structs=expand(config['data_dir']+'structures/{code}/high_cif', code=codes)
    output:
        num=outdir+'/db/db_num.tsv',
        seqs=temp(outdir+'/db/db_num_seqs.tsv'),
        structs=temp(outdir+'/db/db_num_structs.tsv')
    conda: "../envs/sp_utils.yaml"
    shell:'''
find {input.structs} -name "*cif.gz" | rev | cut -d/ -f3 | rev | \
sort | uniq -c | awk -v OFS='\\t' '{{$1=$1;print $2,$1}}' > {output.structs}
echo {input.seqs} | tr ' ' '\\n' | seqkit stats -T -b -X - | sed 's/.fa//g' | \
cut -f1,4 | awk 'NR>1' > {output.seqs}

csvtk join -H -t -f1 {output.structs} {output.seqs} > {output.num}
'''

rule make_fastas:
    input: config['data_dir']+'structures/{code}/high_cif'
    output: config['data_dir']+'fastas/{code}.fa'
    params: config['min_len']
    threads: 8
    conda: "../envs/sp_python.yaml"
    shell:'''
python workflow/scripts/cif2fasta.py -i {input} -o {output} -c {threads} -l {params}
'''

rule make_taxidmap:
    input: rules.make_fastas.output
    output: temp(outdir+'/db/{code}.taxid')
    params:
        taxid=lambda wcs: str(input_dict[wcs.code]['taxid'])
    shell:'''
grep ">" {input} | sed 's/>//' | awk '{{print $0"\\t"{params.taxid}}}' > {output}
'''


