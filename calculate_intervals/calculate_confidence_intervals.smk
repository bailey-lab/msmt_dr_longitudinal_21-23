'''
calculates confidence intervals for mutation prevalences and then reformats
validated and candidate k13 mutations into a table.
'''

input_folder='../AA_table_visualization'
output_folder='confidence_interval_outputs'

rule all:
	input:
		validated_key=output_folder+'/validated_non-k13_prevalences.tsv',
		validated_k13=output_folder+'/validated_k13_prevalences.tsv'

rule calculate_confidence_intervals:
	'''
	calculates confidence intervals and sends outputs into a table.
	'''
	input:
		prevalence_table=input_folder+'/{year}_{threshold}_output/prevalence_summary.tsv'
	output:
		confidence_intervals=output_folder+'/{year}_{threshold}_CIs/{year}_{threshold}_confidence_intervals.tsv'
	conda:
		'statsmodels.yaml'
	script:
		'scripts/calculate_confidence_intervals.py'

rule validated_prevalences:
	'''
	reformats confidence intervals/prevalences associated with candidate and
	validated k13 mutations and aggregates these across years into a single
	table.
	'''
	input:
		all_intervals=expand(output_folder+'/{year}_{threshold}_CIs/{year}_{threshold}_confidence_intervals.tsv', year=['2021', '2022', '2023'], threshold=['10_3', '3_1'])
	params:
		interval_folder=output_folder
	output:
		validated_prevalences=output_folder+'/validated_k13_prevalences.tsv'
	script:
		'scripts/generate_validated_prevalences_v2.py'

rule key_DR_prevalences:
	'''
	reformats confidence intervals/prevalences associated with validated non-k13
	mutations and aggregates these across years into a single table.
	'''
	input:
		all_intervals=expand(output_folder+'/{year}_{threshold}_CIs/{year}_{threshold}_confidence_intervals.tsv', year=['2021', '2022', '2023'], threshold=['10_3', '3_1'])
	params:
		interval_folder=output_folder
	output:
		validated_prevalences=output_folder+'/validated_non-k13_prevalences.tsv'
	script:
		'scripts/generate_key_DR_prevalences_v2.py'