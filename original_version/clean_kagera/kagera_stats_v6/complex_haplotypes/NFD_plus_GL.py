'''
Answers the question of how many individuals have the N86, F184, D1246 alleles
of MDR1 (both as homozygotes and as any presence). Similarly answers the
question of how many individuals have the DHFR 164L allele and the DHPS 581G
allele (both of which seemed to co-localize with each other in geographic data).

Approach: for each sample, checks AA tables for presence of 'N' AND 'F' AND 'D'-
resets all values in corresponding sample to '0' to the minimum count associated
with each allele (simple presence) and resets all values in corresponding sample
to '0' if any value is not equal to the coverage (homozygous).
'''

import pandas as pd
import numpy as np

def reformat_columns(sample_names, header_df, NFD_any, NFD_pure, super_any, super_pure, output_path):
    result = pd.DataFrame({
        'MDR1_NFD_any': NFD_any,
        'MDR1_NFD_pure': NFD_pure,
        'DHFR164L_DHPS581G_any': super_any,
        'DHFR164L_DHPS581G_pure': super_pure
    }, index=sample_names)
    result = result.reset_index().rename(columns={'index': 'Sample ID'})
    new_header = pd.MultiIndex.from_arrays(header_df.values)
    result.columns = new_header
    result.to_csv(output_path, index=False)


def process_haplotypes(coverage_path, alternate_path, output_path):
    # Read the tables (skipping the first 6 rows of metadata/headers to get to numeric data)
    cov_raw = pd.read_csv(coverage_path, header=None)
    alt_raw = pd.read_csv(alternate_path, header=None)

    # Extract mutation names from row 2 and sample names from row 6 on, column 0
    mutation_names = cov_raw.iloc[2].values
    sample_names = cov_raw.iloc[6:, 0].values

    # Extract numeric data
    cov_data = cov_raw.iloc[6:, 1:].astype(float)
    alt_data = alt_raw.iloc[6:, 1:].astype(float)

    # mutation names will skip the first column (which contains sample names)
    cov_data.columns = mutation_names[1:]
    alt_data.columns = mutation_names[1:]
    cov_data.index = sample_names
    alt_data.index = sample_names

    def get_stats(mut_list):
        '''
        gets the coverage and alternate counts associated with a mutation. If a
        single reference has multiple mutations (e.g. a tr-allelic site),
        returns the maximum coverage among the mutation columns and the sum of
        alternate alleles.
		'''
        cols = [column for column in cov_data.columns if column in mut_list]
        if not cols:
            return pd.Series(0.0, index=cov_data.index), pd.Series(0.0, index=cov_data.index)
        return cov_data[cols].max(axis=1), alt_data[cols].sum(axis=1)

    # 1. MDR1 NFD Logic (N86, 184F, D1246)
    cov_86, alt_86 = get_stats(['mdr1-Asn86Phe', 'mdr1-Asn86Tyr'])
    N_cnt = (cov_86 - alt_86).clip(lower=0)
    
    cov_184, alt_184 = get_stats(['mdr1-Tyr184Phe'])
    F_cnt = alt_184
    
    cov_1246, alt_1246 = get_stats(['mdr1-Asp1246Tyr'])
    D_cnt = (cov_1246 - alt_1246).clip(lower=0)

    # Any NFD
    nfd_any = np.where((N_cnt > 0) & (F_cnt > 0) & (D_cnt > 0), 
                       np.minimum(np.minimum(N_cnt, F_cnt), D_cnt), 0.0)
    #NFD coverage
    nfd_cov = np.where((cov_86 > 0) & (cov_184 > 0) & (cov_1246 > 0), 
                       np.minimum(np.minimum(cov_86, cov_184), cov_1246), 0.0)

    # Pure NFD
    nfd_pure = np.where((cov_86 > 0) & (alt_86 == 0) & (cov_184 > 0) & 
                        (alt_184 == cov_184) & (cov_1246 > 0) & (alt_1246 == 0),
                        np.minimum(np.minimum(cov_86, cov_184), cov_1246), 0.0)

    # 2. DHFR 164L + DHPS 581G Logic
    cov_164, alt_164 = get_stats(['dhfr-ts-Ile164Leu'])
    cov_581, alt_581 = get_stats(['dhps-Ala581Gly'])

    # Any 164L/581G
    super_any = np.where((alt_164 > 0) & (alt_581 > 0), 
                         np.minimum(alt_164, alt_581), 0.0)
    # Pure 164L/581G
    super_pure = np.where((cov_164 > 0) & (alt_164 == cov_164) & 
                          (cov_581 > 0) & (alt_581 == cov_581),
                          np.minimum(cov_164, cov_581), 0.0)

    super_cov = np.where((cov_164 > 0) & (cov_581 > 0), 
                         np.minimum(cov_164, cov_581), 0.0)

    # Create and save result
    header_df = pd.read_csv('header_file.csv', header=None)
    reformat_columns(sample_names, header_df, nfd_any, nfd_pure, super_any, super_pure, output_path+'_alternate.csv')
    reformat_columns(sample_names, header_df, nfd_cov, nfd_cov, super_cov, super_cov, output_path+'_coverage.csv')

# Usage
twenty_one_cov='../cleaned_filtered_AA_tables/2021_AA_tables/coverage_AA_table.csv'
twenty_one_alt='../cleaned_filtered_AA_tables/2021_AA_tables/alternate_AA_table.csv'
twenty_two_cov='../cleaned_filtered_AA_tables/2022_AA_tables/coverage_AA_table.csv'
twenty_two_alt='../cleaned_filtered_AA_tables/2022_AA_tables/alternate_AA_table.csv'
twenty_three_cov='../cleaned_filtered_AA_tables/2023_AA_tables/coverage_AA_table.csv'
twenty_three_alt='../cleaned_filtered_AA_tables/2023_AA_tables/alternate_AA_table.csv'

process_haplotypes(twenty_one_cov, twenty_one_alt, '2021_haplotype_counts_AA')
process_haplotypes(twenty_two_cov, twenty_two_alt, '2022_haplotype_counts_AA')
process_haplotypes(twenty_three_cov, twenty_three_alt, '2023_haplotype_counts_AA')