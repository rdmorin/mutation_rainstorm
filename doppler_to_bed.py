'''
Doppler to BED - Convert Doppler output to BED files

Author: Matthew Nguyen, 2020
'''

import argparse as ap
import pandas as pd
import logging

logger = logging.getLogger()

if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Doppler to BED\n' +
                               'Copyright (C) 2020 Matthew Nguyen',
                               formatter_class=ap.RawTextHelpFormatter)

    parser.add_argument('-ll', '--loglevel', type=str, default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='Set the logging level')
    parser.add_argument('input', metavar='DOPPLER_OUTPUT', help='Doppler file to filter')
    parser.add_argument('-o', '--output', metavar='FILENAME', help='Output file prefix or file name')

    param = parser.parse_args()

    logging.basicConfig(level=param.loglevel,
                        format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s: %(message)s',
                        datefmt='%I:%M:%S %p')

    # conditionally determine whether to add .bed extension
    if "bed" not in str(param.output):
        out_file = str(param.output) + '.BED'
    else:
        out_file = str(param.output)

    in_doppler = pd.read_csv(param.input, '\t')
    out_bed = pd.DataFrame(in_doppler[['chromosome', 'leftPosition', 'rightPosition', 'mutPerKb']])

    if "chr" not in out_bed['chromosome'][1]:
        out_bed['chromosome'] = 'chr' + out_bed['chromosome'].astype(str)

    out_bed.to_csv(out_file, sep='\t', index=False, header=False)
    