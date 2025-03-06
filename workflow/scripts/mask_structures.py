#!/usr/bin/env python
from Bio import SeqIO
import pandas as pd
import glob


if __name__ == '__main__':
    min_lddt = snakemake.params[1]

    out_seqs = []

    for record in SeqIO.parse(snakemake.input[0], "fasta"):
        gene_id = record.name
        json_conf = glob.glob(snakemake.params[0]+'/*/confidence/'+gene_id+'-confidence_v4.json.gz')[0]

        confidence = pd.read_json(json_conf)
        low_confidence_sites = list(confidence[confidence['confidenceScore']<min_lddt]['residueNumber'])

        for site in low_confidence_sites:
            record.seq = record.seq[:site-1] + '-' + record.seq[site:]

        out_seqs.append(record)

    SeqIO.write(out_seqs, snakemake.output[0], "fasta")
