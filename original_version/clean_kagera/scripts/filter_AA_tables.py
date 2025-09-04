'''
Removes excess columns from the AA tables to shrink them down to a manageable
size. In this case, we're keeping only mutations belonging to dhfr, crt, mdr1,
dhps, k13, and cytb (pulled in from snakemake params)
'''

year=snakemake.params.year
target_genes=set(snakemake.params.target_genes)
cleaned_coverage=snakemake.input.cleaned_coverage
cleaned_reference=snakemake.input.cleaned_reference
cleaned_alternate=snakemake.input.cleaned_alternate
filtered_coverage=snakemake.output.filtered_coverage
filtered_reference=snakemake.output.filtered_reference
filtered_alternate=snakemake.output.filtered_alternate

def filter_genes(input_table, output_table):
	good_columns=[0]
	output_file=open(output_table, 'w')
	for line_number, line in enumerate(open(input_table)):
		line=line.strip().split(',')
		if line_number==1:
			for column_number, gene in enumerate(line):
				if gene in target_genes:
					good_columns.append(column_number)
			break
	for line_number, line in enumerate(open(input_table)):
		input_line=line.strip().split(',')
		output_line=[]
		for column_number, column in enumerate(input_line):
			if column_number in good_columns:
				output_line.append(column)
		output_file.write(','.join(output_line)+'\n')

filter_genes(cleaned_coverage, filtered_coverage)
filter_genes(cleaned_reference, filtered_reference)
filter_genes(cleaned_alternate, filtered_alternate)