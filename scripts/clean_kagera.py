'''
iterates through the no whitespace sample names in the 'replacement' tables, and
makes replacements or removals in metadata, AA tables, and fastq dictionaries.
Remaining keys that are not found in both metadata and fastq are removed. Any AA
dict keys that are not in metadata and fastq dicts are also removed. Remaining
keys (present in metadata and fastq or present in all three dicts) are renamed -
7 letters followed by 3 numbers are trimmed to this length only. 3 letters
followed by all digits are also kept. Flags any samples that don't follow this
convention for additional parsing. Dashes added between first 4 letters, next 3
letters, and last 3 numbers (and between first 3 letters and remaining digits).
Year-MSMT-1 appended to all sample names at the ends (because some sample names
are redundant across years). For every 'new' key renamed in this way, looks up
the original sample name from the original metadata, aa, or fastq file (values
of dictionaries). Iterates through original metadata, AA tables, and fastq
paths, and replaces any found value in original file with new key name.
'''
import pickle

metadata_dict=snakemake.input.metadata_dict_path
aa_dict=snakemake.input.aa_dict_path
fastq_dict=snakemake.input.fastq_dict_path
aa_replacements=snakemake.input.aa_replacements
metadata_replacements=snakemake.input.metadata_replacements
fastq_replacements=snakemake.input.fastq_replacements
original_metadata=snakemake.input.original_metadata
original_coverage=snakemake.input.original_coverage
original_reference=snakemake.input.original_reference
original_alternate=snakemake.input.original_alternate
original_fastq=snakemake.input.original_fastq

year=snakemake.params.year

cleaned_metadata=snakemake.output.cleaned_metadata
cleaned_coverage=snakemake.output.cleaned_aa_coverage
cleaned_reference=snakemake.output.cleaned_aa_reference
cleaned_alternate=snakemake.output.cleaned_aa_alternate
cleaned_fastq_script=snakemake.output.cleaned_fastq_script

metadata_dict=pickle.load(open(metadata_dict, 'rb'))
aa_dict=pickle.load(open(aa_dict, 'rb'))
fastq_dict=pickle.load(open(fastq_dict, 'rb'))

def replace_keys(sample_dict, replacement_file, dict_type):
	'''
	replaces (or removes) any keys that are marked in the replacement file for
	replacement or removal.
	'''
	for line in open(replacement_file):
		original, replacement=line.strip().split('\t')
#		print('year is', year, 'type is', dict_type)
#		print('original is', original, 'replacement is', replacement)
		if original in sample_dict and 'remove' not in replacement:
			sample_dict[replacement]=sample_dict.pop(original)
		elif original in sample_dict and 'remove' in replacement:
			sample_dict.pop(original)
		elif original not in sample_dict:
			print('Error', year, dict_type, original, replacement)
	return sample_dict

def remove_nonoverlapping(metadata_dict, aa_dict, fastq_dict):
	'''
	removes any metadata dict keys that don't overlap with the fastq dict, and
	any fastq keys that don't overlap with metadata, and any aa_dict keys that
	are not in the intersection of fastq and metadata.
	'''
	metadata_samples=set(metadata_dict.keys())
	aa_samples=set(aa_dict.keys())
	fastq_samples=set(fastq_dict.keys())
	overlappers=metadata_samples&fastq_samples
	#print('year is', year)
	#print('metadata count is', len(metadata_samples))
	#print('fastq count is', len(fastq_samples))
	#print('aa count is', len(aa_samples))
	#print('overlap count is', len(overlappers))
	metadata_dict={overlapper:metadata_dict[overlapper] for overlapper in overlappers}
	aa_dict={overlapper:aa_dict[overlapper] for overlapper in overlappers if overlapper in aa_dict}
	fastq_dict={overlapper:fastq_dict[overlapper] for overlapper in overlappers}
	return metadata_dict, aa_dict, fastq_dict

def reformat_names(sample_dict, dict_type):
	'''
	reformats names from AAAAAAADDDxxx to AAAA-AAA-DDD-year-MSMT-1 or from
	AAADDDDDDD to AAA-DDDDDDD-year-MSMT-1, where A is an alphabetical character
	and D is a digit character and xxx are additional characters that we wish to
	trim off. Checks for any names that don't fit this convention and outputs
	them for further examination.
	'''
	samples=list(sample_dict.keys())
	for sample in samples:
		new_name=''
		if sample[:7].isalpha and sample[7:10].isdigit():
			new_name=f'{sample[:4]}-{sample[4:7]}-{sample[7:10]}-{year}-MSMT-1'
		elif sample[:3].isalpha and sample[3:].isdigit():
			new_name=f'{sample[:3]}-{sample[3:]}-{year}-MSMT-1'
		else:
			print(sample, year, dict_type, 'flagged for further cleaning')
		if new_name:
			sample_dict[new_name]=sample_dict.pop(sample)
	return sample_dict

def output_new_csv(original_dict, original_csv, new_csv, dict_type, header_count):
	'''
	copies entries matching the renamed and pared down input dictionary to a
	new cleaned metadata (or AA) file.
	'''
	new_csv=open(new_csv, 'w')
	flipped_dict={original_dict[key]:key for key in original_dict}
	count=0
	for line_number, line in enumerate(open(original_csv)):
		line=line.strip().split(',')
		if line[0] in flipped_dict:
			line[0]=flipped_dict[line[0]]
			new_csv.write(','.join(line)+'\n')
			count+=1		
		elif line_number<header_count:
			new_csv.write(','.join(line)+'\n')
	if count!=len(original_dict):
		print(year, dict_type, 'missing some expected entries')

def output_new_fastq(fastq_dict, original_fastq, cleaned_fastq_script):
	'''
	Because the original fastq data is stored on multiple separate computers,
	the output of this function is a shell script that can be used to rename and
	copy fastq files of interest.
	'''
	flipped_dict={fastq_dict[key]:key for key in fastq_dict}
	cleaned_fastq_script=open(cleaned_fastq_script, 'w')
	count=0
	#print('flipped is', flipped_dict)
	for line in open(original_fastq):
		line=line.strip()
		original_sample=line.split('/')[-1]
		if line.endswith('_R1_001.fastq.gz'):
			suffix='_R1_001.fastq.gz'
		elif line.endswith('_R2_001.fastq.gz'):
			suffix='_R2_001.fastq.gz'
		sample_lookup=original_sample.replace('_R1_001.fastq.gz', '').replace('_R2_001.fastq.gz', '')
		#print('sample lookup is', sample_lookup)
		if sample_lookup in flipped_dict:
			cleaned_fastq_script.write(f'cp {line} {flipped_dict[sample_lookup]}{suffix}\n')
			count+=1
	if count!=len(fastq_dict)*2:
		print('unexpected fastq count', year, count, len(fastq_dict))

metadata_dict=replace_keys(metadata_dict, metadata_replacements, 'metadata')
aa_dict=replace_keys(aa_dict, aa_replacements, 'aa')
fastq_dict=replace_keys(fastq_dict, fastq_replacements, 'fastq')
metadata_dict, aa_dict, fastq_dict=remove_nonoverlapping(metadata_dict, aa_dict, fastq_dict)
metadata_dict=reformat_names(metadata_dict, 'metadata')
aa_dict=reformat_names(aa_dict, 'aa')
fastq_dict=reformat_names(fastq_dict, 'fastq')
output_new_csv(metadata_dict, original_metadata, cleaned_metadata, 'metadata', 1)
output_new_csv(aa_dict, original_coverage, cleaned_coverage, 'coverage', 6)
output_new_csv(aa_dict, original_reference, cleaned_reference, 'reference', 6)
output_new_csv(aa_dict, original_alternate, cleaned_alternate, 'alternate', 6)
output_new_fastq(fastq_dict, original_fastq, cleaned_fastq_script)