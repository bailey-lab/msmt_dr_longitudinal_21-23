'''
calculates confidence intervals for mutation prevalences and then reformats
validated and candidate k13 mutations into a table.
'''

output_folder='confidence_interval_outputs_exact_95p_avgs'

complex_haps=['mdr1-NFD1ANY', 'mdr1-NFD1PUR', 'sulfa-SUP1ANY', 'sulfa-SUP1PUR']


rule all:
	input:
		validated_prevalences=output_folder+'/complex_haplotype_prevalences_3_1_v2.tsv'

rule calculate_confidence_intervals:
	'''
	calculates confidence intervals and sends outputs into a table.
	'''
	input:
		prevalence_table='{year}_{threshold}_output/prevalence_summary.tsv'
	output:
		confidence_intervals=output_folder+'/{year}_{threshold}_CIs/{year}_{threshold}_confidence_intervals.tsv'
	conda:
		'statsmodels.yaml'
	script:
		'scripts/calculate_confidence_intervals_v2.py'


rule validated_prevalences:
	'''
	reformats confidence intervals/prevalences associated with candidate and
	validated k13 mutations and aggregates these across years into a single
	table.
	'''
	input:
		all_intervals=expand(output_folder+'/{year}_{threshold}_CIs/{year}_{threshold}_confidence_intervals.tsv', year=['2021', '2022', '2023'], threshold=['3_1'])
	params:
		interval_folder=output_folder,
		key_muts=complex_haps,
		threshold='3_1'
	output:
		validated_prevalences=output_folder+'/complex_haplotype_prevalences_3_1_v2.tsv'
	script:
		'scripts/generate_key_DR_prevalences_v2.py'