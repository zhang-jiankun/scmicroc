#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
from biomart import BiomartServer

# filter_coding_genes = True

# Downloaded data from biomart
server = BiomartServer("http://uswest.ensembl.org/biomart")
hsapiens = server.datasets['hsapiens_gene_ensembl']
response = hsapiens.search({'attributes':['ensembl_gene_id', 'gene_biotype', 'ensembl_transcript_id']})
hsapiens_gene_ensembl = []
for line in response.iter_lines():
    line = line.decode('utf-8')
    hsapiens_gene_ensembl.append(line.split('\t'))
ensembl_genes = pd.DataFrame(hsapiens_gene_ensembl, columns=['ensembl_gene_id', 'gene_biotype', 'ensembl_transcript_id'])

chroms = [ 'chr'+str(i) for i in np.arange(1,23)]
chroms = np.append(chroms, 'chrX')

# read data downloaded from UCSC table browser
knownGenes = pd.read_csv(
    '~/data/hsp_genome/GENCODEv32_hg38_knownGene.gz',
    compression='gzip',
    sep='\t'
)

knownGenes.columns = map(lambda x: x.replace('hg38.knownGene.', ''), knownGenes.columns)
knownGenes.columns = map(lambda x: x.replace('hg38.kgXref.', ''), knownGenes.columns)
knownGenes['length'] = knownGenes['txEnd'] - knownGenes['txStart']
knownGenes['#name'] = knownGenes['#name'].str.slice(0, 15)
# Each gene is represented by its longest transcript isoform.
knownGenes = knownGenes.sort_values(['geneSymbol', 'length']).drop_duplicates(['geneSymbol'], keep='last')
knownEnsemblGenes = pd.merge(knownGenes, ensembl_genes, left_on=['#name'], right_on=['ensembl_transcript_id'], how='inner').drop(columns=['ensembl_transcript_id'])
knownEnsemblGenes.to_csv('./ensembl_genes.tsv',index=False, sep='\t')