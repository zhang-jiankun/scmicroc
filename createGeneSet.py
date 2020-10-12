#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd

# Downloaded from UCSC table browser
knownGenes = pd.read_csv(
    '~/data/hsp_genome/GENCODEv32_hg38_knownGene.gz',
    compression='gzip',
    sep='\t'
)

knownGenes.columns = map(lambda x: x.replace('hg38.knownGene.', ''), knownGenes)
knownGenes['length'] = knownGenes['End'] - knownGenes['Start']
# Each gene is represented by its longest transcript isoform.
knownGenes = knownGenes.sort_values(['geneSymbol', 'length']).drop_duplicates(['geneSymbol'], keep='last')
