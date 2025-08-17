'''
Because our initial analyses contained both Kagera samples and samples beyond
Kagera that are not yet fully sequenced, this script removes these additional
samples from the AA tables.
'''
import pickle

metadata_path=snakemake.input.kagera_metadata
AA_path=snakemake.input.AA_path
fastq_path=snakemake.input.fastq_path
year=snakemake.params.year
aa_dict_path=snakemake.output.aa_dict_path
fastq_dict_path=snakemake.output.fastq_dict_path
metadata_dict_path=snakemake.output.metadata_dict_path

odd_AA=open(snakemake.output.odd_AA, 'w')
odd_fastq=open(snakemake.output.odd_fastq, 'w')
stats_file=open(snakemake.output.stats, 'w')

def get_metadata_samples(metadata_path, year):
	sample_dict={}
	for line_number, line in enumerate(open(metadata_path)):
		line=line.strip().split(',')
		if line_number>0 and line[1]==year:
			original_sample=line[0]
			new_sample=original_sample.replace('-', '').replace(' ', '').replace('2021MSMT1', '')
			sample_dict[new_sample]=original_sample
	return sample_dict

#print('valid samples are', valid_samples)
def get_AA_samples(AA_path):
	sample_dict={}
	for line in open(AA_path):
		original_sample=line.strip().split(',')[0]
		new_sample=original_sample.replace('-', '').replace(' ', '').replace('2023MSMT1', '').replace('2021MSMT1', '')
		good_sample=False
		if new_sample.startswith('KG'):
			good_sample=True
		elif new_sample.startswith('KIT') and new_sample[3:].isdigit():
			good_sample=True
		elif new_sample.startswith('KTR') and new_sample[3:].isdigit():
			good_sample=True
		elif new_sample.startswith('NKR') and new_sample[3:].isdigit():
			good_sample=True
		elif new_sample.startswith('RUB') and new_sample[3:].isdigit():
			good_sample=True
		elif new_sample.startswith('RUK') and new_sample[3:].isdigit():
			good_sample=True
		if good_sample:
			sample_dict[new_sample]=original_sample
	return sample_dict

def get_fastq_samples(fastq_path):
	sample_dict={}
	for line in open(fastq_path):
		original_sample=line.strip().split('/')[-1].replace('_R2_001.fastq.gz', '').replace('_R1_001.fastq.gz', '')
		new_sample=original_sample.replace('-', '').replace(' ', '').replace('2023MSMT1', '').replace('2022MSMT1', '').replace('2021MSMT1', '')
		if '2022' in year and new_sample[-3:].isalpha(): #2022 samples got an extra 3 letters added at the end for some reason
			new_sample=new_sample[:-3]
		sample_dict[new_sample]=original_sample
	return sample_dict

metadata_dict=get_metadata_samples(metadata_path, year)
fastq_dict=get_fastq_samples(fastq_path)
AA_dict=get_AA_samples(AA_path)

metadata_samples=set(metadata_dict.keys())
fastq_samples=set(fastq_dict.keys())
AA_samples=set(AA_dict.keys())
stats_file.write(f'\n\nyear is {year}\n')
stats_file.write(f'AA sample size is {len(AA_samples)}\n')
stats_file.write(f'AA with metadata is {len(AA_samples&metadata_samples)}\n')
stats_file.write(f'AA with metadata and fastq is {len(AA_samples&metadata_samples&fastq_samples)}\n')
stats_file.write(f'fastq with metadata is {len(fastq_samples&metadata_samples)}\n')
stats_file.write(f'extra fastq is {len(fastq_samples-metadata_samples)}\n')
pickle.dump(AA_dict, open(aa_dict_path, 'wb'))
pickle.dump(fastq_dict, open(fastq_dict_path, 'wb'))
pickle.dump(metadata_dict, open(metadata_dict_path, 'wb'))
for extra_AA in AA_samples-metadata_samples:
	odd_AA.write(extra_AA+'\n')
for missing_fastq in (AA_samples&metadata_samples)-fastq_samples:
	odd_fastq.write(missing_fastq+'\n')