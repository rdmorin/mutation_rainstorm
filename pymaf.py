"""
Helper functions to work with MAF files. Based on maftools for R.

Author: Matthew Nguyen, 2019
"""

import pandas as pd


"""
Class storing summarizations of MAF files.
"""
class MAF():

    def __init__(self, maf):
        self.nonSyn_df, self.syn_df = self.read_maf(maf)
        self.variants_per_sample = self.variants_per_sample(self.nonSyn_df)

    """
    Reads MAF file, and returns nonsynonymous variants in a dataframe.
    """
    def read_maf(self, maf):
        maf_df = pd.read_csv(maf, sep='\t')

        # Remove synonymous variants
        nonSilentVars = {"Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation",
                         "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"}
        nonSyn_df = maf_df[maf_df['Variant_Classification'].isin(nonSilentVars)]
        syn_df = maf_df[~maf_df['Variant_Classification'].isin(nonSilentVars)]

        return nonSyn_df, syn_df


    """
    Returns dataframe with number of variants per sample.
    """
    def variants_per_sample(self, maf_df):
        variant_count_df = maf_df.groupby('Tumor_Sample_Barcode').size().sort_values(ascending=False)
        return variant_count_df
