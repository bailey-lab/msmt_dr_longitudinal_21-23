'''
filters all metadata to retrieve only lines corresponding to samples from the
kagera region
'''

original_metadata=snakemake.input.all_metadata
output_file=open(snakemake.output.kagera_metadata, 'w')

for line_number, line in enumerate(open(original_metadata)):
	split_line=line.strip().split(',')
	if line_number>0:
		if split_line[2]=='Kagera':
			output_file.write(line)
	else:
		output_file.write(line)