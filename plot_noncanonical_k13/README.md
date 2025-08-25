# plot_noncanonical_k13

The goal of this script is to generate data underlying Figure 2 - namely, which
samples have at least one noncanonical K13 propeller domain mutation.

## Inputs and Outputs
The snakemake program plot_noncanonical_k13.smk executes the python script
plot_noncanonical_k13.py which takes canonical_AA_tables (copied from
clean_kagera/kagera_stats_v6/cleaned_filtered_AA_tables) and generates
noncanonical_AA_tables as output.

## Filtering criteria/algorithm
The script iterates through the AA tables and looks for k13 mutations in the
propeller domain (defined as position 350 or greater) and eliminates any that
belong to known k13 mutations, which include: 'k13-Ala675Val', 'k13-Arg622Ile',
'k13-Cys580Tyr', 'k13-Pro574Leu', 'k13-Val568Gly', 'k13-Arg561His',
'k13-Pro553Leu', 'k13-Ile543Thr', 'k13-Arg539Thr', 'k13-Gly538Val',
'k13-Asn537Ile', 'k13-Asn537Asp', 'k13-Pro527His', 'k13-Arg515Lys',
'k13-Tyr493His', 'k13-Ala481Val', 'k13-Met476Ile', 'k13-Cys469Phe',
'k13-Cys469Tyr', 'k13-Asn458Tyr', 'k13-Gly449Ala', 'k13-Phe446Ile',
'k13-Pro441Leu', 'k13-Ala578Ser'

Among these mutations, only those that pass filtering thresholds of coverage 10
and alternate count 3 are kept. If a sample has at least one mutation that
passes these filters, the mutation with the highest alternate count is chosen to
represent the sample, and if a sample has mutations that pass the coverage
filter but not the alternate count filter, the first mutation that passes the
coverage filters is arbitrarily selected, and the alternate count is set to 0.
Finally, if a sample fails both filters, coverage and alternate counts are both
set to 0. The output AA tables are suitable for downstream parsing with the
AA_table_visualization and plot_static_maps pipelines.