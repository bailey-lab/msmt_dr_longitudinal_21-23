'''
These scripts are needed because this study focuses on the Kagera region,
whereas our original metadata includes additional regions that are still being
actively sequenced. The originally analyzed amino acid tables, metadata, and
fastq files included some of these additional regions but we don't yet have
enough sequencing data to write them up yet. This pipeline removes samples that
are not within the Kagera region. Additionally, the sample nomenclature changed
somewhat between years, and between metadata, amino acid count tables, and fastq
files. This pipeline cleans up and standardizes this nomenclature. Only the
cleaned outputs of the pipeline are included with the github repository.

This pipeline is of course not generalizable - the cleaning programs are built
around the observed idiosyncrasies of the underlying data.

Part 1 of this pipeline has the main goal of flagging AA sample names that
aren't in the metadata and any AA sample names that are present in the metadata
but missing from the fastq samples.

Part 2 uses the flagged samples from part 1 plus a table of manual corrections
to fix the names of samples in the metadata, fastq, and AA table files.
'''
output_folder='kagera_stats_v6'
rule all:
	input:
		odd_AA=expand(output_folder+'/odd_AA_tables/extra_{year}_AA_samples.txt', year=['2021', '2022', '2023'])

rule filter_kagera_metadata:
	'''
	filters metadata to only include Kagera samples
	'''
	input:
		all_metadata='/home/alfred/msmt_re_analysis_with_cleaned_metadata/v3_08-11-25_official_ms_github/msmt_dr_longitudinal_21-23/metadata/merged_corrected_renamed_samples_metadata.csv'
	output:
		kagera_metadata=output_folder+'/kagera_metadata.csv'
	script:
		'scripts/filter_kagera_metadata.py'

rule compare_file_lists:
	'''
	finds which samples from the metadata file are also present in AA tables
	and fastq files. Flags any AA table entries that are missing from metadata,
	and any metadata entries that are present in AA tables but missing from
	fastq files. The metadata dict is overwritten somewhat redundantly but the
	computational cost is small.
	'''
	input:
		fastq_path='fastq_filenames/{year}_fastq_files.txt',
		AA_path='AA_tables/{year}_AA_tables/coverage_AA_table.csv',
		kagera_metadata=output_folder+'/kagera_metadata.csv'
	params:
		year='{year}',
	output:
		metadata_dict_path=output_folder+'/sample_dicts/metadata_dict_{year}.pkl',
		aa_dict_path=output_folder+'/sample_dicts/AA_dict_{year}.pkl',
		fastq_dict_path=output_folder+'/sample_dicts/fastq_dict_{year}.pkl',
		odd_AA=output_folder+'/odd_AA_tables/extra_{year}_AA_samples.txt',
		odd_fastq=output_folder+'/odd_fastq_files/missing_{year}_fastq_samples.txt',
		stats=output_folder+'/sample_comparison_stats/{year}_sample_comparison_stats.txt'
	script:
		'scripts/compare_file_lists.py'

