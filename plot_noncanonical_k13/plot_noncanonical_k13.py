'''
outputs an AA table with a single noncanonical k13 column. Uses the following
steps:
1. iterates through a single line of original AA coverage table and notes k13
column numbers that are not A578S or within blacklist and where 'coverage'
surpasses threshold.
2. iterates through same line of original AA alternate table and chooses the
column from step 1 that has the highest 'alternate' count (or sets this to '0'
if these parameters aren't met).
3. Outputs the coverage and alternate value of this mutation as the
representative noncanonical mutation from this sample.

The 'blacklist' of known K13 mutations and the output AA table header are
hardcoded.
'''

coverage_input=snakemake.input.coverage_input
alternate_input=snakemake.input.alternate_input
metadata_path=snakemake.input.metadata
coverage_threshold=snakemake.params.coverage_threshold
alternate_threshold=snakemake.params.alternate_threshold
propeller_start=snakemake.params.propeller_start
coverage_final=snakemake.output.coverage_final
alternate_final=snakemake.output.alternate_final
coverage_propeller=snakemake.output.coverage_propeller
alternate_propeller=snakemake.output.alternate_propeller
distribution_file=snakemake.output.distributions

coverage_list=[line.strip().split(',') for line in open(coverage_input)]
alternate_list=[line.strip().split(',') for line in open(alternate_input)]

blacklist=['k13-Ala675Val', 'k13-Arg622Ile', 'k13-Cys580Tyr', 'k13-Pro574Leu',
'k13-Val568Gly', 'k13-Arg561His', 'k13-Pro553Leu', 'k13-Ile543Thr',
'k13-Arg539Thr', 'k13-Gly538Val', 'k13-Asn537Ile', 'k13-Asn537Asp',
'k13-Pro527His', 'k13-Arg515Lys', 'k13-Tyr493His', 'k13-Ala481Val',
'k13-Met476Ile', 'k13-Cys469Phe', 'k13-Cys469Tyr', 'k13-Asn458Tyr',
'k13-Gly449Ala', 'k13-Phe446Ile', 'k13-Pro441Leu', 'k13-Ala578Ser']
labels=['Gene ID', 'Gene', 'Mutation Name', 'ExonicFunc', 'AA Change',
'Targeted']
output_header=['PF3D7_1343700', 'k13', 'k13-all350all', 'missense_variant',
'all350all', 'No']


def filter_missense(header, value_list):
	'''
	filters AA table to only include missense mutations
	'''
	filtered_columns=[0]
	filtered_list=[]
	filtered_header=[]
	for line in value_list:
		filtered_line=[line[0]]
		for value_number, value in enumerate(line):
			mutation_type, mut=header[3][value_number], header[2][value_number]
			if mutation_type=='missense_variant':
				filtered_line.append(value)
				filtered_columns.append(value_number)
		filtered_list.append(filtered_line)
	filtered_columns=set(filtered_columns)
	for line in header:
		filtered_line=[]
		for value_number, value in enumerate(line):
			if value_number in filtered_columns:
				filtered_line.append(value)
		filtered_header.append(filtered_line)
	return filtered_list, filtered_header

def get_covered_mutations(blacklist, coverage_threshold, header, line):
	'''
	only returns mutations in the propeller domain of k13 that are not in a
	blacklist of known k13 mutations and have coverage above threshold.
	'''
	covered_columns=[]
	for column_number, column in enumerate(line):
		if column_number>0 and int(float(column))>=coverage_threshold:
			mutation=header[2][column_number]
			gene='-'.join(mutation.split('-')[:-1])
			aa_change=mutation.split('-')[-1]
			pos=int(aa_change[3:-3])
			if mutation not in blacklist and pos>=propeller_start and gene=='k13':
				covered_columns.append(column_number)
#	if line[0]=='KGKA-BKG-222-2021-MSMT-1':
#		print('in cov function')
#		print('covered column numbers are', covered_columns)
#		print('covered column values are', [line[column] for column in covered_columns])
	return covered_columns

def get_alternate_mutations(alternate_threshold, coverage_columns, line):
	'''
	finds the mutation among the 'covered' mutations that has the highest
	alternate count
	'''
	max_alt_count=0
	highest_alt=-1
	for column_number in coverage_columns:
		alt_count=int(float(line[column_number]))
		if alt_count>max_alt_count and alt_count>=alternate_threshold:
			highest_alt=column_number
#	if line[0]=='KGKA-BKG-222-2021-MSMT-1':
#		print('in alt function')
#		print('covered column numbers are', coverage_columns)
#		print('alternate column values are', [line[column] for column in coverage_columns])
	return highest_alt

def parse_tables(coverage_values, alternate_values, missense_header):
	results=[]
	all_covered=[0]
	for line_number, cov_line in enumerate(coverage_values):
		alt_line=alternate_values[line_number]
		covered_columns=get_covered_mutations(blacklist, coverage_threshold, missense_header, cov_line)
		all_covered.extend(covered_columns)
		alt_column=get_alternate_mutations(alternate_threshold, covered_columns, alt_line)
		best_cov, best_alt=0, 0
		sample=cov_line[0]
		if alt_column>0:
			best_alt=alt_line[alt_column]
			best_cov=cov_line[alt_column]
		elif covered_columns:
			print('covered columns is', covered_columns)
			best_cov=cov_line[covered_columns[0]]
		print('sample is', sample, 'best cov is', best_cov, 'best alt is', best_alt)
		results.append([sample, best_cov, best_alt])
	return results, set(all_covered)

def write_values(all_covered, value_list, output_file):
	'''
	A helper function for write_covered - iterates through a list of values and
	prints out any columns that are in all_covered
	'''
	for line in value_list:
		output_line=[]
		for column_number in all_covered:
			output_line.append(line[column_number])
		output_file.write(','.join(output_line)+'\n')

def write_covered(all_covered, missense_header, value_list, output_path):
	'''
	outputs AA tables of all covered mutations in the propeller domain of k13
	'''
	output_file=open(output_path, 'w')
	write_values(all_covered, missense_header, output_file)
	write_values(all_covered, value_list, output_file)

def special_sort(mut_list):
	sorting_list=[]
	for mut_district in mut_list:
		mut=mut_district.split('_')[0]
		gene='-'.join(mut.split('-')[:-1])
		pos=int(mut.split('-')[-1][3:-3])
		sorting_list.append([gene, pos, mut_district])
	return [item[-1] for item in sorted(sorting_list)]

def make_district_dict(metadata_path):
	'''
	looks up the district associated with each sample
	'''
	district_dict={}
	for line_number, line in enumerate(open(metadata_path)):
		if line_number>0:
			line=line.strip().split(',')
			district_dict[line[0]]=line[4]
	return district_dict


def get_distributions(all_covered, missense_header, coverage_list, alternate_list, output_path):
	'''
	gets the number of samples that have each type of mutation
	'''
	output_file=open(output_path, 'w')
	count_dict={}
	total_count=0
	district_dict=make_district_dict(metadata_path)
	coverage_dict={}
	for line_number, line in enumerate(missense_header):
		if line_number==2:
			muts=line
			print('muts are', muts)
	for line_number, cov_line in enumerate(coverage_list):
		district=district_dict[cov_line[0]]
		alt_line=alternate_list[line_number]
		for column_number in all_covered[1:]:
			mut_district=muts[column_number]+'_'+district
			cov_value=int(float(cov_line[column_number]))
			alt_value=int(float(alt_line[column_number]))
			if cov_value>=coverage_threshold and alt_value>=alternate_threshold:
				count_dict[mut_district]=count_dict.setdefault(mut_district, 0)+1
			elif cov_value>=coverage_threshold:
				coverage_dict[district]=coverage_dict.setdefault(district, 0)+1
	for mut in special_sort(list(count_dict.keys())):
		output_file.write('\t'.join(mut.split('_'))+'\t'+str(count_dict[mut])+'\n')
	for district in coverage_dict:
		output_file.write(district+' total\t'+str(coverage_dict[district])+'\n')

def write_final(output_path, results, index_number):
	output_file=open(output_path, 'w')
	for header_number, header_line in enumerate(output_header):
		output_file.write(labels[header_number]+','+header_line+'\n')
	for result_line in results:
		output_file.write(f'{result_line[0]},{result_line[index_number]}\n')


def all_covered_sort(all_covered, missense_header):
	'''
	sorts columns by amino acid position and returns the sorted list. Ignores
	unparseable first column.
	'''
	mut_list=[]
	for column_number in all_covered:
		if column_number!=0:
			mut=missense_header[2][column_number]
			gene='-'.join(mut.split('-')[:-1])
			pos=int(mut.split('-')[-1][3:-3])
			mut_list.append([gene, pos, column_number])
	print('mut list is', mut_list)
	return [0]+[item[-1] for item in sorted(mut_list)]


coverage_values=coverage_list[6:]
alternate_values=alternate_list[6:]
header=coverage_list[:6]
coverage_values, missense_header=filter_missense(header, coverage_values)
alternate_values, junk=filter_missense(header, alternate_values)
results, all_covered=parse_tables(coverage_values, alternate_values, missense_header)
all_covered=all_covered_sort(all_covered, missense_header)
write_covered(all_covered, missense_header, coverage_values, coverage_propeller)
write_covered(all_covered, missense_header, alternate_values, alternate_propeller)
write_final(coverage_final, results, 1)
write_final(alternate_final, results, 2)
get_distributions(all_covered, missense_header, coverage_values, alternate_values, distribution_file)