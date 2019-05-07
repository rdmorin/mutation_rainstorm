'''
Rainstorm calculation - cohort-wide variation of rainfall plots based on a MAF file from many cancer genomes
Based on the R version created by Ryan Morin and Aixiang Jang, 2017

Author: Matthew Nguyen, 2019
'''

import argparse as ap
import pyranges as pr
import pandas as pd
import multiprocessing
import os
import logging

from Bio import AlignIO

logger = logging.getLogger()

if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Rainstorm\n' +
                               'Copyright (C) 2019 Ryan Morin, Aixiang Jang, Matthew Nguyen',
                               formatter_class=ap.RawTextHelpFormatter)

    parser.add_argument('-ll', '--loglevel', type=str, default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='Set the logging level')
    parser.add_argument('-m', '--input_maf', type=str, metavar='INPUT_MAF',
                        help='MAF file containing mutation calls from many patient genomes')
    parser.add_argument('-o', '--output_base_name', type=str, metavar='OUTPUT_BASE_NAME',
                        help='Specify a base file name prefix for all outputs')
    parser.add_argument('-c', '--cpu_num', type=int, metavar='CPU_NUM',
                        help='Set to number of CPUs you would like to use to perform calculation in '
                             'parallel (consumes lots of RAM)')
    parser.add_argument('-g', '--genome_fai', type=str, metavar='GENOME_FAI', required=True,
                        help='Provide the corresponding fasta index for the genome you use. Must '
                             'match the chromosome naming style used in your MAF!')
    parser.add_argument('-p', '--plot', action='store_true', default=False,
                        help='Produce rainstorm plot for each chromosome')
    parser.add_argument('-M', '--max_mut', type=int, metavar='MAX_MUT',
                        help='Genomes skipped if their total mutation load exceeds this value')
    parser.add_argument('-k', '--off_by', type=int, metavar='OFF_BY',
                        help='Take mean of the distance to the k closest mutations to determine '
                             'rainstorm distance value')
    parser.add_argument('-b', '--calc_background', action='store_true', default=False,
                        help='If you have done this once for a cohort, you can reload the result in '
                             'future runs by using this parameter')

    param = parser.parse_args()

    logging.basicConfig(level=param.loglevel,
                        format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s: %(message)s',
                        datefmt='%I:%M:%S %p')

    if not os.path.exists(param.genome_fai):
        logger.error("Missing genome")
        exit(1)

    genomeDetails = pd.read_csv(param.genome_fai, sep='\t', header=None)
    print(genomeDetails)
    chrlengths = genomeDetails[1]
    chrlengths.index = genomeDetails[0]



