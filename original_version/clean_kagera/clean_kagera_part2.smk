'''
uses manual corrections made at the end of part 1 (with the help of lookup_
candidates.py and the odd_fastq_files and odd_AA_tables outputs of step 1) to
fix the 2021, 2022, and 2023 metadata, AA tables, and fastq file entries. For
'''

output_folder='kagera_stats_v6'

rule all:
	input:
		filtered_coverage=expand(output_folder+'/cleaned_filtered_AA_tables/{year}_AA_tables/coverage_AA_table.csv', year=['2021', '2022', '2023'])

rule clean_kagera:
	input:
		metadata_dict_path=output_folder+'/sample_dicts/metadata_dict_{year}.pkl',
		aa_dict_path=output_folder+'/sample_dicts/AA_dict_{year}.pkl',
		fastq_dict_path=output_folder+'/sample_dicts/fastq_dict_{year}.pkl',
		metadata_replacements=output_folder+'/manual_corrections/metadata_corrections_{year}.tsv',
		aa_replacements=output_folder+'/manual_corrections/AA_corrections_{year}.tsv',
		fastq_replacements=output_folder+'/manual_corrections/fastq_corrections_{year}.tsv',
		original_metadata=output_folder+'/kagera_metadata.csv',
		original_coverage='AA_tables/{year}_AA_tables/coverage_AA_table.csv',
		original_reference='AA_tables/{year}_AA_tables/reference_AA_table.csv',
		original_alternate='AA_tables/{year}_AA_tables/alternate_AA_table.csv',
		original_fastq='fastq_filenames/{year}_fastq_files.txt'
	params:
		year='{year}'
	output:
		cleaned_metadata=output_folder+'/cleaned_metadata/{year}_cleaned_kagera_metadata.csv',
		cleaned_aa_coverage=output_folder+'/cleaned_AA_tables/{year}_AA_tables/coverage_AA_table.csv',
		cleaned_aa_reference=output_folder+'/cleaned_AA_tables/{year}_AA_tables/reference_AA_table.csv',
		cleaned_aa_alternate=output_folder+'/cleaned_AA_tables/{year}_AA_tables/alternate_AA_table.csv',
		cleaned_fastq_script=output_folder+'/fastq_filenames/{year}_cleaned_fastq_commands.sh',
		cleaned_sample_sheet=output_folder+'/sample_sheets/{year}_sample_sheet.tsv'
	script:
		'scripts/clean_kagera.py'

rule filter_tables:
	input:
		cleaned_coverage=output_folder+'/cleaned_AA_tables/{year}_AA_tables/coverage_AA_table.csv',
		cleaned_reference=output_folder+'/cleaned_AA_tables/{year}_AA_tables/reference_AA_table.csv',
		cleaned_alternate=output_folder+'/cleaned_AA_tables/{year}_AA_tables/alternate_AA_table.csv',
	params:
		target_genes=['dhfr-ts', 'crt', 'mdr1', 'dhps', 'k13', 'cytb'],
		year='{year}'
	output:
		filtered_coverage=output_folder+'/cleaned_filtered_AA_tables/{year}_AA_tables/coverage_AA_table.csv',
		filtered_reference=output_folder+'/cleaned_filtered_AA_tables/{year}_AA_tables/reference_AA_table.csv',
		filtered_alternate=output_folder+'/cleaned_filtered_AA_tables/{year}_AA_tables/alternate_AA_table.csv',
	script:
		'scripts/filter_AA_tables.py'