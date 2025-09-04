'''
It feels a little silly to have a snakemake program that runs a single script,
but it is nice to see inputs and outputs clearly labeled and to have automatic
output folder creation, as well as the automatic and parallelized analysis of
all 3 years and automatic determination of whether the script needs to be rerun
or not.
'''

rule all:
	input:
		coverage_output=expand('noncanonical_AA_tables/{year}_AA_tables/coverage_AA_table.csv', year=['2021', '2022', '2023'])

rule plot_noncanonical_k13:
	input:
		coverage_input='../clean_kagera_v2/AA_tables/{year}_AA_tables/coverage_AA_table.csv',
		alternate_input='../clean_kagera_v2/AA_tables/{year}_AA_tables/alternate_AA_table.csv',
		metadata='../clean_kagera_v2/kagera_rerun_cleaned/cleaned_metadata/{year}_cleaned_kagera_metadata.csv'
	params:
		coverage_threshold=10,
		alternate_threshold=3,
		propeller_start=350
	output:
		coverage_final='noncanonical_AA_tables/{year}_AA_tables/coverage_AA_table.csv',
		alternate_final='noncanonical_AA_tables/{year}_AA_tables/alternate_AA_table.csv',
		coverage_propeller='k13_novel_propeller_AA_tables/{year}_AA_tables/coverage_AA_table.csv',
		alternate_propeller='k13_novel_propeller_AA_tables/{year}_AA_tables/alternate_AA_table.csv',
		distributions='noncanonical_AA_tables/{year}_AA_tables/{year}_distributions.tsv'
	script:
		'plot_noncanonical_k13.py'
