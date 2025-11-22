'''
calculates confidence intervals for mutation prevalences and then reformats
validated and candidate k13 mutations into a table.
'''

input_folder='../AA_table_visualization'
output_folder='confidence_interval_outputs_exact_95p_avgs'

key_non_k13=['dhfr-ts-Ala16Val', 'dhfr-ts-Asn51Ile', 'dhfr-ts-Cys59Arg',
'dhfr-ts-Ile164Leu', 'dhfr-ts-Ser108Asn', 'dhfr-ts-Ser108Thr',
'mdr1-Asn1042Asp', 'mdr1-Asn86Tyr', 'mdr1-Asp1246Tyr', 'mdr1-Ser1034Cys',
'mdr1-Tyr184Phe', 'crt-Val73Leu', 'crt-Ala220Ser', 'crt-Arg371Ile',
'crt-Asn75Glu', 'crt-Cys101Phe', 'crt-Cys72Ser', 'crt-Gln271Glu',
'crt-Gly353Val', 'crt-His97Leu', 'crt-His97Tyr', 'crt-Ile218Phe',
'crt-Ile356Thr', 'crt-Lys76Thr', 'crt-Met343Leu', 'crt-Met74Ile',
'crt-Phe145Ile', 'crt-Thr93Ser', 'dhps-Ala437Gly', 'dhps-Ala581Gly',
'dhps-Ala613Ser', 'dhps-Ala613Thr', 'dhps-Ile431Val', 'dhps-Lys540Glu',
'dhps-Ser436Ala', 'dhps-Ser436Phe', 'crt-Asn326Ser']

key_k13=['k13-Pro441Leu', 'k13-Phe446Ile', 'k13-Gly449Ala',
'k13-Asn458Tyr', 'k13-Cys469Phe', 'k13-Cys469Tyr', 'k13-Met476Ile',
'k13-Ala481Val', 'k13-Tyr493His', 'k13-Arg515Lys', 'k13-Pro527His',
'k13-Asn537Asp', 'k13-Asn537Ile', 'k13-Gly538Val', 'k13-Arg539Thr',
'k13-Ile543Thr', 'k13-Pro553Leu', 'k13-Arg561His', 'k13-Val568Gly',
'k13-Pro574Leu', 'k13-Cys580Tyr', 'k13-Arg622Ile', 'k13-Ala675Val']

rule all:
	input:
		validated_key=output_folder+'/validated_non-k13_prevalences_v2.tsv',
		validated_k13=output_folder+'/validated_k13_prevalences_v2.tsv',
		validated_k13_noCI=output_folder+'/validated_k13_prevalences_noCI.tsv',
		validated_noCI=output_folder+'/validated_non-k13_prevalences_noCI.tsv'


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
#		'scripts/calculate_confidence_intervals.py'
		'scripts/calculate_confidence_intervals_v2.py'


rule validated_prevalences:
	'''
	reformats confidence intervals/prevalences associated with candidate and
	validated k13 mutations and aggregates these across years into a single
	table.
	'''
	input:
		all_intervals=expand(output_folder+'/{year}_{threshold}_CIs/{year}_{threshold}_confidence_intervals.tsv', year=['2021', '2022', '2023'], threshold=['10_3', '3_1'])
	params:
		interval_folder=output_folder,
		key_muts=key_k13
	output:
		validated_prevalences=output_folder+'/validated_k13_prevalences_v2.tsv'
	script:
		'scripts/generate_key_DR_prevalences_v2.py'

rule key_DR_prevalences:
	'''
	reformats confidence intervals/prevalences associated with validated non-k13
	mutations and aggregates these across years into a single table.
	'''
	input:
		all_intervals=expand(output_folder+'/{year}_{threshold}_CIs/{year}_{threshold}_confidence_intervals.tsv', year=['2021', '2022', '2023'], threshold=['10_3', '3_1'])
	params:
		interval_folder=output_folder,
		key_muts=key_non_k13
	output:
		validated_prevalences=output_folder+'/validated_non-k13_prevalences_v2.tsv'
	script:
		'scripts/generate_key_DR_prevalences_v2.py'

rule validated_no_CI:
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
		validated_prevalences=output_folder+'/validated_k13_prevalences_noCI.tsv'
	script:
		'scripts/generate_validated_prevalences_noCI.py'

rule key_DR_no_CI:
	'''
	reformats confidence intervals/prevalences associated with validated non-k13
	mutations and aggregates these across years into a single table.
	'''
	input:
		all_intervals=expand(output_folder+'/{year}_{threshold}_CIs/{year}_{threshold}_confidence_intervals.tsv', year=['2021', '2022', '2023'], threshold=['10_3', '3_1'])
	params:
		interval_folder=output_folder
	output:
		validated_prevalences=output_folder+'/validated_non-k13_prevalences_noCI.tsv'
	script:
		'scripts/generate_key_DR_prevalences_noCI.py'