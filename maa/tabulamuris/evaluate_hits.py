# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/10/18
content:    Build a model to evaluate the cell type enrichment with combinations
            of hits (from KS or other sources).
'''
import os
import sys
import argparse
import numpy as np
import pandas as pd
import xarray as xr
from Bio import AlignIO

import matplotlib.pyplot as plt
import seaborn as sns

from maa.tabulamuris.base import DatasetTM, tissues, load_mouse_human_homology, load_hits


def select_model(models):
    # Sort models according to some decent criteria
    # TODO: refine these
    models = models.copy()
    models['PPV+TPR'] = models['PPV'] + models['TPR']
    models_sorted = models.sort_values(['PPV+TPR', 'PPV', 'TPR'], ascending=False)
    del models_sorted['PPV+TPR']
    del models['PPV+TPR']

    # Take first pick with decent immunofluorescence antibodies
    ab = ('Supported', 'Enhanced', 'Approved')
    models_sorted['select'] = False
    for _, model in models_sorted.iterrows():
        if (model['HPA_Rel_IF1'] in ab) and (model['HPA_Rel_IF2'] in ab):
            models_sorted.loc[_, 'select'] = True
    models_select = models_sorted.loc[models_sorted['select']]
    del models_select['select']
    if len(models_select) == 0:
        raise ValueError('No model has decent antibodies')
    return models_select


def find_protein_refseqs(gene, mho):
    # Here the manually curated exceptions, mouse gene name to human HomoloGeneID
    exceptions = {'Bst2': 48277}

    tmp_mouse = mho.loc[(gene, 'mouse')]
    hgi = exceptions.get(gene, tmp_mouse['HomoloGeneID'])
    tmp_human = mho.query('(HomoloGeneID == @hgi) & (Organism == "human")')
    if tmp_human.shape[0] < 1:
        raise ValueError('Homolog not found for mouse gene {:}'.format(gene))
    tmp_human = tmp_human.iloc[0]
    refseqs_mouse = tmp_mouse['ProteinRefSeqIDs'].split(',')
    refseqs_human = tmp_human['ProteinRefSeqIDs'].split(',')
    return {
        'mouse': refseqs_mouse,
        'human': refseqs_human,
        }


def download_protein_sequences(refseqs):
    from Bio import Entrez, SeqIO
    seqs = []
    refseqs_all = refseqs['mouse'] + refseqs['human']

    # Request Entrez
    Entrez.email = 'fabio.zanini@fastmail.fm'
    handle = Entrez.efetch(db="protein", id=refseqs_all, rettype="gb", retmode="text")
    seqs = list(SeqIO.parse(handle, 'genbank'))

    # Assign common species
    for seq in seqs:
        if seq.name in refseqs['mouse']:
            seq.annotations['OrganismCommon'] = 'mouse'
        elif seq.name in refseqs['human']:
            seq.annotations['OrganismCommon'] = 'human'
        else:
            raise ValueError('Where is this refseq coming from?')

    # Get rid of predicted stuff
    seqs_not_predicted = []
    for seq in seqs:
        if not seq.description.startswith('PREDICTED'):
            seqs_not_predicted.append(seq)

    # Check that we have at least one per species
    n_seqs = {'mouse': 0, 'human': 0}
    for seq in seqs_not_predicted:
        if seq.annotations['OrganismCommon'] == 'mouse':
            n_seqs['mouse'] += 1
        elif seq.annotations['OrganismCommon'] == 'human':
            n_seqs['human'] += 1
    if (n_seqs['mouse'] < 1) or (n_seqs['human'] < 1):
        raise ValueError('Missing a species (non predicted only)')
    return seqs_not_predicted


def align_multiple_sequences(seqs):
    import subprocess as sp
    from Bio import SeqIO, AlignIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet.IUPAC import protein
    from Bio.Align import MultipleSeqAlignment
    fn_tmp = '/tmp/tmp.fasta'
    fn_tmp2 = '/tmp/tmp_aligned.fasta'
    seqs_muscle = []
    for i, seq in enumerate(seqs):
        name = 'seq{:}'.format(i)
        seqm = SeqRecord(seq.seq, id=name, name=name, description='')
        seqs_muscle.append(seqm)
    SeqIO.write(seqs_muscle, fn_tmp, 'fasta')
    sp.run(
        'muscle -in {:} -out {:} -diags'.format(fn_tmp, fn_tmp2),
        shell=True,
        )
    ali_muscle = AlignIO.read(fn_tmp2, 'fasta')
    ali = []
    for seq_muscle in ali_muscle:
        i = int(seq_muscle.id[3:])
        rec = seqs[i]
        seq = SeqRecord(
            seq_muscle.seq,
            id=rec.annotations['OrganismCommon']+'-'+rec.id,
            name=rec.name,
            description=rec.description,
            )
        seq.annotations = rec.annotations
        ali.append(seq)
    ali = MultipleSeqAlignment(ali)
    os.remove(fn_tmp)
    os.remove(fn_tmp2)
    return ali


def build_doublepositive(ds, ct, hits, tissue):
    '''Build model to enrich cell types by double positive gates.

    Args:
        ds (DatasetTM): the normalized dataset to analyze
        ct (str): the cell type to enrich
        hits (DataFrame): the table of double positive hits

    Returns:
        model of the gates
    '''
    def test_enrichment(ds, ct, gene1, gene2, th1, th2, df):
        # NOTE: this is an AND gate, i.e. double positives (square)
        df['prediction'] = (ds.counts.loc[gene1] >= th1) & (ds.counts.loc[gene2] >= th2)
        res = df.groupby(['identity', 'prediction']).count()['tmp'].unstack().fillna(0)
        if True not in res.columns:
            res[True] = 0
        if False not in res.columns:
            res[False] = 0
        return res

    # Restrict to the hits
    ds = ds.query_features_by_name(hits.index)

    # Test all combinations
    thresholds = np.logspace(0, 5, 15)

    df = ds.samplesheet[['cell_ontology_class']] == ct
    df.rename(columns={'cell_ontology_class': 'identity'}, inplace=True)
    df['tmp'] = 1
    models = []
    for i1, gene1 in enumerate(hits.index):
        for i2, gene2 in enumerate(hits.index[:i1]):
            for th1 in thresholds:
                for th2 in thresholds:
                    res = test_enrichment(ds, ct, gene1, gene2, th1, th2, df)
                    mod = {
                        'gene1': gene1,
                        'gene2': gene2,
                        'threshold1': th1,
                        'threshold2': th2,
                        'prevalence': 1.0 * res.loc[True].sum() / res.values.sum(),
                        'TPR': 1.0 * res.loc[True, True] / res.loc[True].sum(),
                        'FPR': 1.0 * res.loc[False, True] / res.loc[False].sum(),
                        'PPV': 1.0 * res.loc[True, True] / res.loc[:, True].sum(),
                        'HPA_Rel_IH1': hits.loc[gene1, 'HPA_Rel_IH'],
                        'HPA_Rel_IH2': hits.loc[gene2, 'HPA_Rel_IH'],
                        'HPA_Rel_IF1': hits.loc[gene1, 'HPA_Rel_IF'],
                        'HPA_Rel_IF2': hits.loc[gene2, 'HPA_Rel_IF'],
                        'cell_type': ct,
                        'tissue': tissue,
                    }
                    # Enrichment: PPV over prevalence
                    mod['enrichment'] = mod['PPV'] / mod['prevalence']
                    models.append(mod)
            models = pd.DataFrame(models)
            return models


if __name__ == '__main__':

    pa = argparse.ArgumentParser(
        description='Test loading various tabula muris tissues',
        )
    pa.add_argument(
        '--tissue',
        required=True,
        choices=tissues,
        help='Tissue of origin',
        )
    pa.add_argument(
        '--hit-method',
        choices=['KS'],
        default='KS',
        help='Method used to find the hits',
        )
    args = pa.parse_args()

    mho = load_mouse_human_homology()

    print('Load data for {:}'.format(args.tissue))
    ds = DatasetTM(
        counts_table=args.tissue,
        samplesheet='alltissues',
        featuresheet='alltissues',
        )

    print('Parse Human Protein Atlas')
    hpa = ds.parse_hpa()

    print('QC cells')
    ds.query_samples_by_metadata(
        '(is_annotated == True) and (coverage >= 50000)',
        inplace=True)

    print('Normalize')
    ds.counts.normalize(inplace=True)

    print('Load hits')
    hits_alltypes = load_hits(args.tissue, method=args.hit_method)
    cell_types = np.unique(hits_alltypes['cell_type'])

    print('Look for double positive gates, they are the most robust')
    for ct in cell_types:
        print(ct)
        hits = hits_alltypes.query('(cell_type == @ct) & (diff > 0)')

        print('Build models')
        models = build_doublepositive(ds, ct, hits, args.tissue)

        print('Select best model')
        models_select = select_model(models)

        for _, model in models_select.iterrows():
            try:
                print('Align the proteins and check how close they are')
                refseqs1 = find_protein_refseqs(model['gene1'], mho)
                refseqs2 = find_protein_refseqs(model['gene2'], mho)
                ali1 = align_multiple_sequences(download_protein_sequences(refseqs1))
                ali2 = align_multiple_sequences(download_protein_sequences(refseqs2))
                break
            except ValueError:
                continue
        else:
            raise ValueError('No suitable model found for {:}, {:}'.format(
                args.tissue,
                ct,
                ))

        print('Write output folder')
        fdn = '../../data/tabula_muris/FACS/models/double_positive_{:}_{:}/'.format(
            args.tissue,
            ct.replace(' ', '_'))
        os.makedirs(fdn, exist_ok=True)
        AlignIO.write(ali1, fdn+'alignment_{:}.fasta'.format(model['gene1']), 'fasta')
        AlignIO.write(ali2, fdn+'alignment_{:}.fasta'.format(model['gene2']), 'fasta')
        model_out = model.rename(columns={
            'threshold1': 'threshold gene 1 (counts per million)',
            'threshold2': 'threshold gene 2 (counts per million)',
            })
        model_out.to_csv(fdn+'model.tsv', sep='\t')

    print('Merge all in a summary file')
    fn = '../../data/tabula_muris/FACS/models/{:}_summary.txt'.format(
        args.tissue)
    with open(fn, 'wt') as f:
        f.write('Summary for antibody search for tissue: {:}\n'.format(args.tissue))
        f.write('This type of models is a double positive AND gate, so it classifies a cell if it stains for both genes.\n')
        f.write('Legend: PPV = positive predictive value, TPR = true positive rate, FPR = false positive rate, HPA_Rel_IF1 = reliability of Human Protein Atlas antibodies for immunocytochemistry of gene 1, HPA_Rel_IF2 = same for gene 2, HPA_Rel_IH1 = same but for immunohistochemistry (tissues), HPA_Rel_IH2 = same but for immunohistochemistry (tissues) and gene 2.\n\n')

        f.write('Models:\n')
        is_first = True
        for ct in cell_types:
            fdn = '../../data/tabula_muris/FACS/models/double_positive_{:}_{:}/'.format(
                args.tissue,
                ct.replace(' ', '_'))
            model = pd.read_csv(fdn+'model.tsv', sep='\t', squeeze=True, header=None, index_col=0)
            # Change order for clarity
            cols_first = ['cell_type', 'gene1', 'gene2', 'PPV', 'FPR', 'TPR']
            cols_new = cols_first + [x for x in model.index if x not in (cols_first + ['tissue'])]
            model = model[cols_new]
            # Nans in the HPA don't format well
            for key in ('IF1', 'IF2', 'IH1', 'IH2'):
                if not isinstance(model['HPA_Rel_'+key], str) and np.isnan(model['HPA_Rel_'+key]):
                    model['HPA_Rel_'+key] = 'not available'

            # Approximate the fractions
            for key in ('TPR', 'FPR', 'PPV'):
                model[key] = '{:.2f}'.format(np.round(float(model[key]), 2))

            # Format enrichment and prevalence
            model['enrichment'] = '{:.1f}'.format(float(model['enrichment']))
            model['prevalence'] = '{:.2%}'.format(float(model['prevalence']))

            # Format thresholds
            model['threshold1'] = '{:.1f}'.format(float(model['threshold1']))
            model['threshold2'] = '{:.1f}'.format(float(model['threshold2']))

            # First columns can get quite long
            fmt = '{:'+str(len(max(cell_types, key=len)) + 1)+'s}'
            if is_first:
                f.write('\t'.join([fmt.format('cell type'), '{:12s}'.format('gene 1'), '{:12s}'.format('gene 2')] + [str(x) for x in model.index[3:]])+'\n')
                is_first = False
            f.write('\t'.join([fmt.format(model['cell_type']), '{:12s}'.format(model['gene1']), '{:12s}'.format(model['gene2'])] + [str(x) for x in model.iloc[3:]])+'\n')
        f.write('\nAlignments:\n')
        for ct in cell_types:
            fdn = '../../data/tabula_muris/FACS/models/double_positive_{:}_{:}/'.format(
                args.tissue,
                ct.replace(' ', '_'))
            model = pd.read_csv(fdn+'model.tsv', sep='\t', squeeze=True, header=None, index_col=0)
            model.index.name = 'column'
            model.name = 'model'
            ali1 = AlignIO.read(fdn+'alignment_{:}.fasta'.format(model['gene1']), 'fasta')
            ali2 = AlignIO.read(fdn+'alignment_{:}.fasta'.format(model['gene2']), 'fasta')
            f.write(model['gene1']+':\n')
            for seq in ali1:
                f.write('{:}\t{:}\n'.format(str(seq.seq), seq.id))
            f.write('\n')


