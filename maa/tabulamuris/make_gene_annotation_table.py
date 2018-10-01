# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/10/18
content:    Make gene annotation table from gene names and biomart export.
'''
import os
import numpy as np
import pandas as pd


if __name__ == '__main__':

    with open('../../data/tabula_muris/FACS/gene_names.csv', 'rt') as f:
        gene_names = f.read().strip('\n').split(',')

    biomart = pd.read_csv(
        '../../data/tabula_muris/FACS/tmp/mart_export.tsv',
        sep='\t',
        index_col=0,
        )
    biomart[biomart.index.name] = biomart.index
    print('Biomart\'s uinque is not really unique, make it unique and restrict to gene list')
    biomartu = []
    genes = set()
    n_genes = 0
    for _, row in biomart.iterrows():
        gene = row['Gene name']
        if gene in genes:
            continue
        if gene not in gene_names:
            continue
        biomartu.append(row)
        genes.add(gene)
        n_genes += 1
        print(n_genes, gene)
    biomartu = pd.DataFrame(biomartu)
    biomartu.set_index('Gene name', drop=False, inplace=True)
    for gene in gene_names:
        if gene not in biomartu.index:
            biomartu.loc[gene] = ''
    biomartu['Gene name'] = biomartu.index
    # Reformat dtypes
    if 'Gene % GC content' in biomartu.columns:
        biomartu['Gene % GC content'] = [float(x) if x else np.nan for x in biomartu['Gene % GC content']]
    biomartu = biomartu.loc[gene_names]
    biomartu.to_csv(
        '../../data/tabula_muris/FACS/metadata_FACS_features.tsv',
        sep='\t',
        index=False,
        )
