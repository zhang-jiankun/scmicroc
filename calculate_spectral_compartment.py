#!/usr/env/bin python


"""
Author: Yueying He (GaoYQ lab)
"""

import sys
import copy, time
import click
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

import straw
import pybedtools as pb

def Laplacian(adjacency):
    """calculate L=D^-0.5 * (A+I) * D^-0.5"""
    degree = np.array(adjacency.sum(1))
    d_hat = sp.diags(np.power(degree, -0.5).flatten())
    return d_hat.dot(adjacency).dot(d_hat).tocoo()

@click.command()
@click.argument('hic', type=str, required=True)
@click.argument('res', type=int, required=True)
@click.argument('chromsize', type=str, required=True)
@click.argument('track', type=str, required=True)
def calculate_spectral_compartment(hic, res, chromsize, track):
    
    chrom_size = dict(pd.read_csv(chromsize, sep='\t', header=None).values)
    chromosomes = list(chrom_size.keys())
    n_chroms = len(chrom_size.keys())

    # chromosome binning using pybedtools
    # w: res, s: res
    chrom_bin = pb.BedTool().window_maker(
        g=chromsize, w=res, s=res
    ).to_dataframe().sort_values(by=['chrom', 'start', 'end'])

    matrix_size = chrom_bin.shape[0]
    # add offset of chromosome 1
    matrix_offset = np.insert(np.cumsum(chrom_bin['chrom'].value_counts()[chromosomes].values)[:-1], 0, 0)
    chrom_bin_offset = dict(zip(chromosomes, matrix_offset))

    chrom_sign = np.zeros(matrix_size)

    k = 0
    for i in chromosomes:
        k += 1
        offset = chrom_bin_offset[i]
        chrom_sign[offset:] = k

    track_sign = chrom_sign.copy()
    track_pos = np.loadtxt(track, dtype=str, usecols=[1,2,3])

    k = 0
    for i in chromosomes:
        k += 1
        track_pos_chr = track_pos[track_pos[:,0]==i][:,[1,2]].astype(int).mean(1)
        track_counts = np.histogram(track_pos_chr, range = (0,np.sum(chrom_sign==k)*res), bins = np.sum(chrom_sign==k))[0]
        track_sign[chrom_sign == k] = track_counts

    hic_matrix_all = sp.coo_matrix(([0],([0],[0])),shape = (matrix_size, matrix_size),dtype = "float32")

    for chr1 in range(n_chroms):
        for chr2 in np.arange(chr1, n_chroms):
            results = straw.straw('NONE', hic, chromosomes[chr1].replace('chr', ''), chromosomes[chr2].replace('chr', ''), 'BP', res)
            row = np.floor(np.array(results[0])/res) + chrom_bin_offset[chromosomes[chr1]]
            col = np.floor(np.array(results[1])/res) + chrom_bin_offset[chromosomes[chr2]]
            hic_matrix = sp.coo_matrix((np.array(results[2]),(row,col)),shape = (matrix_size, matrix_size),dtype = "float32")
            hic_matrix_all = hic_matrix_all + hic_matrix 

    hic_matrix_full = hic_matrix_all + hic_matrix_all.transpose()
    # set upper limits
    q_99 = np.quantile(hic_matrix_full.data, 0.99)
    hic_matrix_full[hic_matrix_full > q_99] = q_99
    s = np.array(hic_matrix_full.sum(1)).reshape(-1)
    hic_matrix_full = hic_matrix_full[s>0]
    hic_matrix_full = hic_matrix_full[:, s>0]

    Lap = Laplacian(hic_matrix_full)
    vals, vecs = eigsh(Lap, k = 100)
    vecs = np.real_if_close(vecs)
    vals = np.real_if_close(vals)

    o = np.argsort(1 - vals)
    vals = 1 - vals[o]
    vecs = vecs[:,o]

    lda = LinearDiscriminantAnalysis().fit_transform(vecs,(track_sign>0)[s>0])
    lda = 1/(1 + np.exp(-lda)) - 0.5
    comp_index = np.repeat(np.nan, s.shape[0])
    comp_index[s > 0] = lda.reshape(-1)

    comp_index_df = chrom_bin.copy()
    comp_index_df['comp_index'] = comp_index
    comp_index_df.to_csv(sys.stdout, index=False, sep='\t')

if __name__ == '__main__':
    calculate_spectral_compartment()
