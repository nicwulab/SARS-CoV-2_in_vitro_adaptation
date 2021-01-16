#!/usr/bin/env python

import sys
import glob
import os
import pandas as pd
import seaborn as sns
from matplotlib import use as mpl_use
mpl_use('agg')
import matplotlib.pyplot as plt
import numpy as np
plt.rc('axes', labelsize=15)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)


def read_bed(bed):
    return pd.read_csv(bed, names = ['chrom','start','end','read_coverage'], sep='\t') \
        .assign(samplename = os.path.basename(os.path.dirname(bed))) \
        .assign(log_cov = lambda d: d.read_coverage.transform(np.log))


def main(mosdepth_bed_file: str, figure: str):
    '''
    plotting coverage plot from a mosdepth bed file

    :param  mosdepth_bed_file: bed filr generated from mosdepth
    :param figure: output figure name
    '''
    df = read_bed(mosdepth_bed_file)

    #make figure
    p = sns.FacetGrid(data = df, col_wrap = 1, col = 'samplename')
    p.map(plt.plot, 'end', 'log_cov')
    p.set_titles(col_template = '{col_name}')
    sns.despine()
    p.set(xlabel = 'Read coverage (log)', ylabel = 'Genomic position')
    p.savefig(figre, bbox_inches='tight')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit('[usage] python %s <mosdepth bed file> <figurename> ' %sys.argv[0])
    main(sus/agrv[1])

