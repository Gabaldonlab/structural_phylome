import glob
import pandas as pd

# User-specified dataset-specific config file
dataset_config_file = config.get("dataset_config_file", None)

# Load fixed config parameters
configfile: "config/params.yaml"

# Check if the dataset-specific config file is provided
if dataset_config_file is not None:
    # Load dataset-specific parameters
    configfile: dataset_config_file

data_dir=config['data_dir']
# dataset specific info
seed=config['seed']
files = ['cif', 'pae', 'confidence']

input_table = pd.read_csv(config['taxids'], sep='\t')
input_table.columns = ['uniprot', 'taxid', 'mnemo']
input_dict = input_table.set_index('uniprot').T.to_dict()

codes = list(input_table['uniprot'])

rule all:
    input:
        # expand('data/ids/{code}.txt', code=codes),
        # expand(data_dir+'{code}/cif', code=codes),
        # expand(data_dir+'{code}/pae', code=codes),
        # expand(data_dir+'{code}/confidence', code=codes),
        expand(data_dir+'structures/{code}/low_cif', code=codes),
        expand(data_dir+'structures/{code}/high_cif', code=codes),
        expand(data_dir+'gffs/{code}.gff', code=codes),
        expand(data_dir+'cath/{code}_cath.tsv', code=codes),
        expand(data_dir+'ids/{code}_ids.tsv', code=codes),
        expand(data_dir+'ids/{code}_3ds.tsv', code=codes),
        data_dir+'meta/'+config["homology_dataset"]+'_uniprot_genomes.tsv'


rule download_pdbs:
    output:
        temp(directory(data_dir+'structures/raw_files/{code}'))
    params:
        taxid=lambda wcs: str(input_dict[wcs.code]['taxid'])
    shell:'''
currdir=$PWD
mkdir -p {output}
cd {output}
gsutil -m cp -r gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-{params.taxid}-* .
cd $currdir
'''

rule download_uniprot_ids:
    output:
        ids=data_dir+'structures/{code}/{code}.txt'
        # rev_ids='data/ids/{code}_rev.txt'
    params:
        taxid=lambda wcs: str(input_dict[wcs.code]['taxid'])
    shell: '''
set +o pipefail;

kingdom=$(echo {params.taxid} | taxonkit reformat -I 1 -f "{{k}}" | cut -f2)
url=ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$kingdom/{wildcards.code}/{wildcards.code}_{params.taxid}.fasta.gz

if curl --head --silent --fail $url 2> /dev/null; then
    echo "Download from reference proteome"
    wget -O - $url | zgrep ">" | cut -f2 -d'|' > {output}
else
    echo "Download from fasta query"
    curl "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28proteome%3A{wildcards.code}%29%29" --output - | \
    zgrep ">" | cut -f2 -d'|' > {output}
fi
'''

rule untar_pdbs:
    input: 
        ids=rules.download_uniprot_ids.output,
        files=rules.download_pdbs.output
    output:
        cif=temp(directory(data_dir+'structures/{code}/cif')),
        pae=directory(data_dir+'structures/{code}/pae'),
        conf=directory(data_dir+'structures/{code}/confidence')
    shell:'''
mkdir -p {output.cif}
mkdir -p {output.pae}
mkdir -p {output.conf}

for tarfile in {input.files}/*tar; do

    ids=$(echo $tarfile | sed 's/.tar//')
    tar --list --file=$tarfile | awk -F- 'NR == FNR {{ a[$0]; next }} $2 in a' {input.ids} /dev/stdin > ${{ids}}.txt
    if [[ $(wc -l <${{ids}}.txt) -ge 1 ]]; then
        grep "model_v4.cif.gz" ${{ids}}.txt > ${{ids}}_cif.txt
        grep "predicted_aligned_error_v4.json.gz" ${{ids}}.txt > ${{ids}}_pae.txt
        grep "confidence_v4.json.gz" ${{ids}}.txt > ${{ids}}_conf.txt
        tar xf $tarfile -C {output.cif} -T ${{ids}}_cif.txt
        tar xf $tarfile -C {output.pae} -T ${{ids}}_pae.txt
        tar xf $tarfile -C {output.conf} -T ${{ids}}_conf.txt
    fi
done
'''


rule move_lowconf:
    input: 
        conf=rules.untar_pdbs.output.conf,
        cif=rules.untar_pdbs.output.cif
    output: 
        stats=data_dir+'structures/{code}/{code}_mean_plddt.tsv',
        low_dir=directory(data_dir+'structures/{code}/low_cif'),
        high_dir=directory(data_dir+'structures/{code}/high_cif')
    params:
        config['low_confidence']
    shell:'''
mkdir -p {output.low_dir}
mkdir -p {output.high_dir}

> {output.stats}
for struct in {input.conf}/*json.gz; do
    protein=$(basename $struct | cut -f2 -d'-')
    avg_lddt=$(zcat $struct | jq '.confidenceScore |  add/length*1000 | round/1000')
    echo -e "$protein\\t$avg_lddt" >> {output.stats}

    if (( $(echo "$avg_lddt < {params}" | bc -l) )); then
        cp {input.cif}/AF-${{protein}}-F1-model_v4.cif.gz {output.low_dir}
    else 
        cp {input.cif}/AF-${{protein}}-F1-model_v4.cif.gz {output.high_dir}
    fi
done
'''

rule download_up_meta:
    input: sps=config["taxids"]
    output: data_dir+'meta/'+config["homology_dataset"]+'_uniprot_genomes.tsv'
    shell: '''
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/STATS -O - | \
tail -n +16 | sed 's/ //' | csvtk join -t -f "Tax_ID,Proteome_ID" {input.sps} - | cut -f1,2,3,7,9 | \
taxonkit reformat -I 2 -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" | awk 'NR>1' | \
csvtk add-header -t -n Proteome_ID,Tax_ID,mnemo,Assembly,Species_name,Kingdom,Phylum,Class,Order,Family,Genus,Species > {output}
'''


rule get_gff:
    output: data_dir+'gffs/{code}.gff'
    shell: '''
wget "https://rest.uniprot.org/uniprotkb/stream?format=gff&query=%28%28proteome%3A{wildcards.code}%29%29" \
-O /dev/stdout | awk 'NF' > {output}  
'''

rule get_CATH:
    output: data_dir+'cath/{code}_cath.tsv'
    shell: '''
wget "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cxref_pfam%2Cxref_gene3d&format=tsv&query=%28%28proteome%3A{wildcards.code}%29%29" -O {output}  
'''

rule get_phygeno:
    output: data_dir+'ids/{code}_ids.tsv'
    shell: '''
wget "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cxref_genetree%2Cxref_hogenom%2Cxref_inparanoid%2Cxref_phylomedb%2Cxref_orthodb%2Cxref_oma%2Cxref_treefam%2Cxref_eggnog&format=tsv&query=%28%28proteome%3A{wildcards.code}%29%29" \
-O {output}  
'''

rule get_3Ds:
    output: data_dir+'ids/{code}_3ds.tsv'
    shell: '''
wget "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cxref_pdb%2Cxref_emdb%2Cxref_bmrb%2Cxref_alphafolddb%2Cxref_pcddb%2Cxref_pdbsum%2Cxref_smr%2Cxref_sasbdb&format=tsv&query=%28proteome%3A{wildcards.code}%29" \
-O {output}  
'''

