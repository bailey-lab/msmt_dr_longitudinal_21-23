'''
calculates confidence intervals for mutation prevalences
'''

input_folder='../AA_table_visualization'
output_folder='confidence_interval_outputs'

rule all:
	input:
		confidence_intervals=expand(output_folder+'/{year}_{threshold}_CIs/{year}_{threshold}_confidence_intervals.tsv', year=['2021', '2022', '2023'], threshold=['10_3', '3_1'])

rule calculate_confidence_intervals:
	input:
		prevalence_table=input_folder+'/{year}_{threshold}_output/prevalence_summary.tsv'
	output:
		confidence_intervals=output_folder+'/{year}_{threshold}_CIs/{year}_{threshold}_confidence_intervals.tsv'
	conda:
		'statsmodels.yaml'
	script:
		'calculate_confidence_intervals.py'