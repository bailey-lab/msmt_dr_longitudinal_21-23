'''
It feels a little silly to have a snakemake program that runs a single script
but it is nice to see inputs and outputs clearly labeled and to have automatic
output folder creation.
'''

rule all:
	input:
		coverage_output=expand('noncanonical_AA_tables/{year}_AA_tables/coverage_AA_table.csv', year=['2021', '2022', '2023'])

rule plot_noncanonical_k13:
	input:
		coverage_input='canonical_AA_tables/{year}_AA_tables/coverage_AA_table.csv',
		alternate_input='canonical_AA_tables/{year}_AA_tables/alternate_AA_table.csv'
	params:
		coverage_threshold=10,
		alternate_threshold=3,
		propeller_start=350
	output:
		coverage_output='noncanonical_AA_tables/{year}_AA_tables/coverage_AA_table.csv',
		alternate_output='noncanonical_AA_tables/{year}_AA_tables/alternate_AA_table.csv'
	script:
		'plot_noncanonical_k13.py'
