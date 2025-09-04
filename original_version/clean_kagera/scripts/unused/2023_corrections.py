'''
The original metadata contained several misformatted 2023 sample names. I've
included the correction script here (but not the original misformatted input
files) so that people can see how these corrections were performed.

This script is still buggy - many of the sample names are poorly formatted in
both the AA tables and the metadata and this current script doesn't deal with
this problem.
'''

#output_file=open('2023_automated_corrections_v3.tsv', 'w')
#output_file.write('original_name\tfastq_name\n')

metadata_file='../merged_corrected_renamed_samples_metadata_kagera.csv'
corrected_metadata_file='../merged_corrected_renamed_samples_metadata_kagera_2023_corrections.csv'
AA_file='../AA_tables/2023_AA_tables/coverage_AA_table.csv'

def reformat_sample(sample):
	new_sample=sample.replace('-', '').replace(' ', '').replace('2023MSMT1', '')
	first_part, second_part, third_part=new_sample[:4], new_sample[4:7], new_sample[7:-3]
	if first_part.isalpha() and second_part.isalpha and third_part.isdigit():
		pass
		#print('valid sample', new_sample)
	elif new_sample[:3].isalpha() and new_sample[3:].isdigit():
		pass
		#print('valid sample', new_sample)
	else:
		print('odd sample', new_sample)
	return new_sample

def get_metadata_samples(metadata_file):
	metadata_samples=[]
	for metadata_line in open(metadata_file):
		metadata_line=metadata_line.strip().split(',')
		if metadata_line[1]=='2023':
			new_sample=reformat_sample(metadata_line[0])
			metadata_samples.append(new_sample)
	return metadata_samples

def get_AA_samples(AA_file):
	AA_samples=[]
	for line_number, line in enumerate(open(AA_file)):
		line=line.strip().split(',')
		if line_number>5:
			new_sample=reformat_sample(line[0])
			AA_samples.append(new_sample)
	return AA_samples

metadata_samples=set(get_metadata_samples(metadata_file))
AA_samples=set(get_AA_samples(AA_file))
print('metadata is', len(metadata_samples))
print('AA is', len(AA_samples))
print('shared is', len(metadata_samples&AA_samples))
print('extra AA is', len(AA_samples-metadata_samples))
print('extra metadata is', len(metadata_samples-AA_samples))
for sample in AA_samples-metadata_samples:
	if sample[:2]=='KG':
		print('extra AA sample is', sample)

#for sample in metadata_samples-AA_samples:
#	if 'BUT276' in sample:
#		print('potential BUT276 metadata sample is', sample)
'''
seen_hits=set([])
good_dict={}
duplicated_set=set([])
for missing_sample in metadata_samples:
	first_four, second_three, last_numbers=False,False,False
	testing=False
	if len(missing_sample.split('-'))>2:
		first_part=missing_sample.split('-')[0]
		second_part=missing_sample.split('-')[1]
		third_part=missing_sample.split('-')[2]
		if len(first_part)==4 and first_part.isalpha():
			first_four=first_part
		else:
			print(missing_sample, 'has invalid first part')
		if len(second_part)==3 and second_part.isalpha() and third_part.isdigit():
			second_three=second_part
			last_numbers=third_part
		elif len(second_part)==6 and second_part[:3].isalpha() and second_part[3:].isdigit():
			second_three=second_part[:3]
			last_numbers=second_part[3:]
		else:
			print(missing_sample, 'has invalid second part')
		if testing:
			print('parts are', first_four, second_three, last_numbers)
		if first_four and second_three and last_numbers:
			hits=[] #these are possible matchers from the AA tables
			for AA_sample in AA_samples:
				if AA_sample.startswith(first_four) and second_three in AA_sample[4:] and last_numbers in AA_sample[7:]:
					hits.append(AA_sample)
			if len(hits)==0:
				print(missing_sample, 'not found in AA samples')
			elif len(hits)==2:
				print(missing_sample, 'found the following hits in AA samples:', hits)
			else:
				good_dict.setdefault(hits[0], []).append(missing_sample)
				if hits[0] in seen_hits:
					duplicated_set.add(hits[0])
				seen_hits.add(hits[0])
	else:
		print(missing_sample, 'has fewer than 3 parts')
'''
#for hit in good_dict:
	#The conditionals are commented out here because if a metadata sample maps
	#to two 'missing' samples, I would like to arbitrarily analyze one of them
	#rather than tossing both.
	#comment the 'if' and 'else' statements back in to ignore duplicated entries
	#if hit not in duplicated_set and len(good_dict[hit])==1:
#	output_file.write(hit+'\t'+good_dict[hit][0]+'\n')
#	else:
#		print('duplicates:', good_dict[hit])
