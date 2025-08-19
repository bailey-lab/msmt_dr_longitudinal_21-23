'''
This pipeline somewhat automates district plotting but has many hard-coded
features specific to Kagera region of Tanzania. Among the useful automations:
1. Automatically detects the mutants from prevalence_summary tables (generated
during variant_graphing step of miptools). 
2. converts mutants to percentages and sample sizes
3. graphs every mutant (in parallel).
4. Automatically filters out regions that are outside boundaries of geographic
regions.

Among the hard-coded features:
1. shape file
2. geographic region
3. neighboring country labels
'''
min_long=30.3
max_long=32
min_lat=-3.4
max_lat=-0.89


def get_muts(prevalence_file):
	columns=open(prevalence_file).readline().strip().split('\t')
	muts=[mut.replace('-', '_') for mut in columns[3:]]
	return muts

year='2021'
threshold='3_1'

prevalence_file='/home/alfred/msmt_re_analysis_with_cleaned_metadata/v3_08-11-25_official_ms_github/msmt_dr_longitudinal_21-23/AA_table_visualization/'+year+'_'+threshold+'_output/prevalence_summary.tsv'
mutations=get_muts(prevalence_file)
print('mutations are', mutations)

rule all:
	input:
		final_files=expand('kagera_graphs_'+year+'_'+threshold+'/png_plots/kagera_{mutation}_'+year+'_'+threshold+'.png', mutation=mutations)

rule convert_prevalence_table:
	input:
		prevalence_table=prevalence_file
	params:
		longitude_min=min_long,
		longitude_max=max_long,
		latitude_min=min_lat,
		latitude_max=max_lat
	output:
		R_muts='R_tables/muts_district_'+year+'.csv'
	script:
		'python_scripts/prevalence_to_R.py'

rule graph_kagera:
	input:
		input_csv='R_tables/muts_district_'+year+'.csv'
	params:
		longitude_min=min_long,
		longitude_max=max_long,
		latitude_min=min_lat,
		latitude_max=max_lat
	output:
		output_png='kagera_graphs_'+year+'_'+threshold+'/png_plots/kagera_{mutation}_'+year+'_'+threshold+'.png',
		output_svg='kagera_graphs_'+year+'_'+threshold+'/svg_plots/kagera_{mutation}_'+year+'_'+threshold+'.svg'
	conda:
		"envs/R_environment.yml"
	script:
		'R_scripts/kagera_only_all_muts.R'