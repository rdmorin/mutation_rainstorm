'''
Rainstorm calculation - cohort-wide variation of rainfall plots based on a MAF file from many cancer genomes
Based on the R version created by Ryan Morin and Aixiang Jang, 2017

Author: Matthew Nguyen, 2019
'''

import pdb
import argparse as ap
import pyranges as pr
import pandas as pd
import numpy as np
import multiprocessing as mp
import os
import logging
import pymaf
import skmisc.loess as loess
import time
import math

pd.options.mode.chained_assignment = None

logger = logging.getLogger()

"""
Correct local mutation rate

Uses a loess model that fits a smoothed curve to the mutation rate across the chromosome
"""
def correctLocalMutrate(positions, distval, model, logged_mutrate):
    # predrate =  predict(model,newdata=data.frame(starts=positions))
    predrate = model.predict(positions)
    adjusted = np.log(distval + predrate.values + logged_mutrate)
    return adjusted


"""
Obtain a value for each mutation that is later scaled for local mutation rate

Compares pairs of genomes and is called when performing a one-to-all comparison
"""
def getMutDists(pos1, pos2):
    # Merge together both lists of positions
    if len(pos1) == 0 or len(pos2) == 0:
        return [np.nan] * len(pos1)

    pos2 = np.append(pos2, 1000000000)  # To ensure last p1 position always gets a diff value
    these = pd.DataFrame({'names': (['p1'] * len(pos1)) + (['p2'] * len(pos2)), 'mut': np.concatenate((pos1, pos2))})

    # Sort on position
    sorted = these.sort_values('mut')
    diffs = sorted.iloc[:-1]
    diffs.loc[:, 'mut'] = np.diff(sorted['mut'])
    # Determine adjacencies in the same genome (to mask out as NA)
    sorted['names_shifted'] = sorted['names'].shift(-1)
    adjacents = sorted.index[sorted['names'] == sorted['names_shifted']].tolist()
    diffs.loc[adjacents, 'mut'] = np.nan
    # keep only the positions with names indicating they derive from a position in p1
    pos1diffs = diffs.loc[diffs['names'] == 'p1', 'mut'].values
    return pos1diffs


"""
Helper function for each row of distsort 
"""
def offby_mutations(x, offby=3, use_mean=True):
    if use_mean:
        return np.mean(x[0:offby])
    else:
        # Original approach is to just return the kth value instead of mean from 1:k
        return x[offby - 1]


"""
Calls getMutDists() on all cases for a single index case (ID)
"""
def getMinDistByGenome(maf, id, chrom, IDs, start, end, offby=3, use_mean=True):
    # Extract mutations in region for this genome and compute the N-closest minimum distance to each variant among all
    # genomes (default N, 2). Self is ignored, nearest genome is ignored.

    thesemut = maf.nonSyn_df.loc[(maf.nonSyn_df['Tumor_Sample_Barcode'] == id) & (maf.nonSyn_df['Chromosome'] == chrom)
                                 & (maf.nonSyn_df['Start_Position'] >= start)
                                 & (maf.nonSyn_df['End_Position'] < end)]['Start_Position'].values
    thesemut = np.sort(thesemut)
    IDs.remove(id)

    all_dists = dict()
    for case in IDs:
        thosemut = maf.nonSyn_df.loc[(maf.nonSyn_df['Tumor_Sample_Barcode'] == case)
                                     & (maf.nonSyn_df['Chromosome'] == chrom)
                                     & (maf.nonSyn_df['Start_Position'] >= start)
                                     & (maf.nonSyn_df['End_Position'] < end)]['Start_Position'].values
        thosemut = np.sort(thosemut)
        all_dists[case] = getMutDists(thesemut, thosemut)

    distmat = pd.DataFrame(np.vstack([i for i in all_dists.values()]))

    # before removing any cases where every value is NA,
    # the indexes in alldists correspond to the mutation positions in thesemut
    # need to mark all-NA positions for removal and removal of corresponding position in thesemut
    # this command removes any patient that contributed only NAs to the matrix

    # Remove positions with only NA
    allna_pos = distmat.columns[distmat.isna().all()].tolist()
    if len(allna_pos) > 0:
        distmat = distmat.drop(columns=allna_pos)
        thesemut = np.delete(thesemut, allna_pos)

    # Remove patients with only NA
    allna_pat = pd.isnull(distmat).all(1).to_numpy().nonzero()[0]
    if len(allna_pat) > 0:
        distmat = distmat.drop(index=allna_pat)

    distsort = np.sort(distmat.values.transpose())
    # print(table(unlist(lapply(distsort,length))==0))
    keepdist = np.apply_along_axis(offby_mutations, 1, distsort, offby=offby, use_mean=use_mean)
    IDs.append(id)
    return pd.DataFrame({'position': thesemut, 'mindist': keepdist})


def runByCaseSmooth_multiprocess(case, maf, model, variants, IDs, chrom, start, end, nathres=0.3, offby=3):
    output = runByCaseSmooth(case, maf, model, variants, IDs, chrom, start, end, nathres, offby)
    return case, output


def runByCaseSmooth(case, maf, model, variants, IDs, chrom, start, end, nathres=0.3, offby=3):
    start_time = time.time()
    stored_all = {'mutdiff': [], 'position': [], 'mutrate': [], 'mutrate_noadj': [], 'patient': []}
    use_mean = True

    these = getMinDistByGenome(maf, case, chrom, IDs.tolist(), start, end, offby=offby, use_mean=use_mean)

    if these.shape[0] == 0:
        logger.info("Skip due to lack of mutations on chromosome {0}".format(case))
        return stored_all

    density_table = these.isna()
    density_table_true = density_table['position'].sum() + density_table['mindist'].sum()
    density_table_false = (len(density_table) - density_table['position'].sum()) + \
                          (len(density_table) - density_table['mindist'].sum())

    if density_table_false > 0 and density_table_true > 0:
        if density_table_false / density_table_true > nathres:
            logger.info("Skip due to high NA count {0}".format(case))
            return stored_all

    tot = these.shape[0]
    genometot = maf.variant_count.loc[maf.variant_count['Tumor_Sample_Barcode'] == case, 'Variants']
    # genometot = maf.full @ variants.per.sample[Tumor_Sample_Barcode == case, Variants]
    ltot = math.log(genometot / 280000000)
    # ltot = log(genometot / 280000000)
    logger.debug("Shifting by {0}, {1}".format(genometot, ltot))
    napos = these['mindist'].index[these['mindist'].apply(np.isnan)].tolist()
    these_keep = these.drop(index=napos)
    nonatot = these_keep.shape[0]

    # Should we get rid of all NA values first? Seems reasonable since they're being counted here in the denominator?
    # Though they are in fact mutations, so maybe not...
    scaled = these_keep['mindist'].apply(lambda x: math.log(x+1) + math.log(genometot) - math.log(280000000))
    # Add one to get rid of the annoying -Inf issue. These are definitely things that need to be retained.
    scaled -= np.median(scaled)
    localadj = correctLocalMutrate(these_keep['position'], these_keep['mindist']+1, model, ltot)
    scaled_localadj = localadj - np.nanmedian(localadj)
    stored_all['mutdiff'] = these_keep['mindist'].tolist()
    stored_all['position'] = these_keep['position'].tolist()
    stored_all['mutrate'] = scaled_localadj.tolist()
    stored_all['mutrate_noadj'] = scaled.tolist()
    stored_all['patient'] = [case] * len(these_keep['position'])
    done_time = time.time()
    overall_time = done_time - start_time
    logger.info("{0} took {1} seconds".format(case, overall_time))
    return stored_all


def viewMeans(bins, cvg):
    # mean_per_bin = []
    mut_index = cvg.values.nonzero()[0]
    # mut_loc = deque()
    mut_loc = []
    loc = 0

    for ind, run in enumerate(cvg.runs):
        loc += run
        if ind in mut_index:
            mut_loc.append(loc)

    # for _,bin in bins.iterrows():
    #     count = 0
    #     while len(mut_loc) > 0 and bin['Start'] <= mut_loc[0] < bin['End']:
    #         count += 1
    #         mut_loc.popleft()
    #     mean = float(count) / (float(bin['End']) - float(bin['Start']))
    #     mean_per_bin.append(mean)

    bins['binned_score'] = 0.0

    for mut in mut_loc:
        bins.loc[(bins['Start'] <= mut) & (bins['End'] > mut), 'binned_score'] += 1

    mask = (bins['binned_score'] != 0.0)
    bins_valid = bins[mask]

    bins.loc[mask, 'binned_score'] /= (bins_valid['End'] - bins_valid['Start'])

    # return mean_per_bin
    return bins['binned_score'].values.tolist()


def binnedAverage(bins, cvg):
    # bins_per_chrom = {chrom: chrom_df for chrom, chrom_df in bins.df.groupby('Chromosome')}
    # means_list = []

    # for chrom in bins_per_chrom.keys():
    #     means_list += viewMeans(bins_per_chrom[chrom], cvg[chrom])
    # bins.binned_score = means_list
    bins.binned_score = viewMeans(bins.df, cvg)
    return bins


def plotRainstorm(points, name):
    # ggplot(points, aes(x=position, y=mutrate, colour=patient), size=1) + geom_point(
    #     alpha=0.2) + theme_classic() + theme(legend.position = "none") + ylim(NA, 0)
    # ggsave(file=name, width=7, height=4)
    return


if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Rainstorm\n' +
                               'Copyright (C) 2019 Ryan Morin, Aixiang Jang, Matthew Nguyen',
                               formatter_class=ap.RawTextHelpFormatter)

    parser.add_argument('-ll', '--loglevel', type=str, default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='Set the logging level')
    parser.add_argument('-i', '--input_maf', type=str, metavar='INPUT_MAF',
                        help='MAF file containing mutation calls from many patient genomes')
    parser.add_argument('-nc', '--nonCoding', action='store_true', default=False,
                        help='Limit to non-coding range only')
    parser.add_argument('-o', '--output_base_name', type=str, metavar='OUTPUT_BASE_NAME',
                        help='Specify a base file name prefix for all outputs')
    parser.add_argument('-c', '--cpu_num', type=int, metavar='CPU_NUM', default=1,
                        help='Set to number of CPUs you would like to use to perform calculation in '
                             'parallel (consumes lots of RAM)')
    parser.add_argument('-g', '--genome_fai', type=str, metavar='GENOME_FAI', default="hg19.ensembl.fa.fai",
                        help='Provide the corresponding fasta index for the genome you use. Must '
                             'match the chromosome naming style used in your MAF!')
    parser.add_argument('-p', '--plot', action='store_true', default=True,
                        help='Produce rainstorm plot for each chromosome')
    parser.add_argument('-M', '--max_mut', type=int, metavar='MAX_MUT', default=50000,
                        help='Genomes skipped if their total mutation load exceeds this value')
    parser.add_argument('-k', '--off_by', type=int, metavar='OFF_BY', default=4,
                        help='Take mean of the distance to the k closest mutations to determine '
                             'rainstorm distance value')
    parser.add_argument('-b', '--calc_background', action='store_true', default=False,
                        help='If you have done this once for a cohort, you can reload the result in '
                             'future runs by using this parameter')
    parser.add_argument('-na', '--nathresh', type=float, metavar='NA_THRESH', default=0.3,
                        help='Threshold for of NA to skip a patient')

    param = parser.parse_args()

    logging.basicConfig(level=param.loglevel,
                        format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s: %(message)s',
                        datefmt='%I:%M:%S %p')

    if not os.path.exists(param.genome_fai):
        logger.error("Missing genome")
        exit(1)

    genomeDetails = pd.read_csv(param.genome_fai, sep='\t', header=None)
    print(genomeDetails)
    chrlengths = genomeDetails[[0, 1]]
    chrlengths.insert(1, 'Start', 0)
    chrlengths.columns = ['Chromosome', 'Start', 'End']
    chrlengths_pr = pr.PyRanges(chrlengths)

    maf = pymaf.MAF(param.input_maf)

    # Get IDs of cases passing the max mutation criteria
    IDs = maf.variant_count[maf.variant_count.Variants < param.max_mut]['Tumor_Sample_Barcode']

    # Choice for both full and non-coding range
    if param.nonCoding:
        variants = {"3'Flank", "IGR", "Intron", "3'UTR", "5'Flank", "5'UTR", "Targeted_Region", "RNA"}
    else:
        variants = maf.codingVars

    goodchrom = chrlengths['Chromosome'].values

    cols = {'Chromosome': }
    snvs_df = maf.nonSyn_df[['Chromosome', 'Start_Position', 'End_Position', 'Tumor_Sample_Barcode']]
    snvs_df.columns = ['Chromosome', 'Start', 'End', 'Tumor_Sample_Barcode']
    snvs_df['End'] += 1

    if not param.calc_background:
        bins_chr = chrlengths_pr.window(100000, tile=True)

        bincounts_all = []
        binstarts_all = []
        bincounts_chrom = []
        binstops_all = []
        binlength = 100000

        for chrom in goodchrom:
            logger.info("Calculating {0}".format(chrom))
            patient = IDs[0]
            snvs_df_subset = snvs_df.loc[snvs_df['Tumor_Sample_Barcode'] == IDs[0]]
            snvs_subset = pr.PyRanges(snvs_df_subset)

            cvg = snvs_subset.coverage()
            npat_tot = snvs_df_subset.shape[0]

            tile = binnedAverage(bins_chr[chrom], cvg[chrom])
            ntile = len(tile[chrom].df['binned_score'])
            testmat = np.empty((ntile, len(IDs)))

            # r <- runsum(cvg, 1) --> doesn't actually do anything?

            for num in range(len(IDs)):
                patient = IDs[num]
                snvs_df_subset = snvs_df.loc[snvs_df['Tumor_Sample_Barcode'] == IDs[0]]

                cvg = snvs_subset.coverage()
                npat_tot = snvs_df_subset.shape[0]

                tile = binnedAverage(bins_chr[chrom], cvg[chrom])
                a = tile[chrom].df['binned_score']
                testmat[:, num] = a[:ntile]

                # r <- runsum(cvg, 1) --> doesn't actually do anything?

            means = np.mean(testmat, axis=1)
            means *= 0.000000001

            # note, natural log scale (harmonize to log10 as per the distance/rainfall approach?)
            logmeans = np.log(means * binlength)
            bincounts_all += logmeans.tolist()
            binstarts_all += tile.df['Start'].tolist()
            binstops_all += tile.df['End'].tolist()
            bincounts_chrom += [chrom for i in range(len(logmeans))]

        logger.info("Done calculating background correction")

        # Load this next file in as a data frame to skip the steps leading up to this line
        cols = {'chrom': bincounts_chrom, 'starts': binstarts_all, 'ends': binstops_all, 'counts': bincounts_all}
        all_df = pd.DataFrame(cols, index=None)
        all_df.to_csv(param.output_base_name+'_background_100k_binned_density.tsv', sep='\t',
                      index=False)
    else:
        all_df = pd.read_csv(param.output_base_name+'_background_100k_binned_density.tsv', sep='\t')

    n = len(IDs) + 22
    '''
    cols=colors()[c(23:n)]
    names(cols) = use.cases
    '''

    for chrom in goodchrom:
        logger.info("Running calculation for {0}".format(str(chrom)))
        start = 1
        end = maf.nonSyn_df.loc[maf.nonSyn_df['Chromosome'] == chrom]['Start_Position'].max()
        data = all_df.loc[(all_df['chrom'] == chrom) & (all_df['counts'] != -np.inf)]
        if data.empty:
            logger.warning("No bins with mutations on chromosome {0}".format(chrom))
            continue
        # model= loess(counts~starts,data=alldf[alldf[,1]== chrom & alldf[,"counts"]!=-Inf,],span=0.01,surface='direct')
        success = False
        # span = 0.01
        span = 0.1
        while not success:
            try:
                # todo: currently loess loops infinitely for span that is too small
                model = loess.loess(data['starts'], data['counts'], span=span, surface='direct')
                model.fit()
                success = True
            except:
                span += 0.02

        if param.cpu_num > 1:
            pool = mp.Pool(processes=param.cpu_num)
            outputs = [pool.apply(runByCaseSmooth_multiprocess, args=(case, maf, model, variants, IDs, chrom, start,
                                                                      end, param.nathresh, param.off_by))
                       for case in IDs]
            all_data_all_patients = dict(outputs)
            # all_data_all_patients = {case: runByCaseSmooth(case, maf, model, variants, IDs, chrom, start, end,
            #                                                param.nathres, param.off_by)
            #                          for case in IDs}  # parallelize this
        else:
            start_time = time.time()
            case_times = {}
            all_data_all_patients = {}
            lu = len(IDs)
            j = 1

            for case in IDs[50:100]:
                all_data_all_patients[case] = runByCaseSmooth(case, maf, model, variants, IDs, chrom, start, end,
                                                              param.nathresh, offby=param.off_by)

                end_time = time.time()
                duration = end_time - start_time
                logger.info("Time for {0} was {1}".format(case, duration))
                logger.info("Done {0} of {1}".format(j, lu))
                j += 1
                case_times[case] = duration

            meantime = np.asarray(list(case_times.values())).mean()
            logger.info("Average time per query genome comparison: {0}".format(meantime))

        # Convert to lists with like elements combined and grouped by patient
        patients = []
        positions = []
        mutrate = []
        unadj = []
        mutdiff = []

        for patient in all_data_all_patients.keys():
            n = len(unadj)
            logger.info("Patient {0}\n------------".format(n))
            unadj += all_data_all_patients[patient]['mutrate_noadj']
            patients += all_data_all_patients[patient]['patient']
            positions += all_data_all_patients[patient]['position']
            mutrate += all_data_all_patients[patient]['mutrate']
            mutdiff += all_data_all_patients[patient]['mutdiff']

        # Ready points for ggplot rendering
        all_counted = pd.DataFrame({'mutrate': mutrate, 'unadj': unadj, 'position': positions, 'patient': patients,
                                    'mutdiff': mutdiff})
        filename = "{0}_rainstorm_k_{1}_mean_{2}.tsv".format(param.output_base_name, param.off_by, chrom)
        all_counted.to_csv(filename, sep='\t', )
        # write.table(allcounted,file=filen,sep="\t",quote=F);
        # plotRainstorm(allcounted,gsub(".tsv",".pdf",filen));
