'''
Removes excess columns from the AA tables to shrink them down to a manageable
size. In this case, we're keeping only mutations belonging to dhfr, crt, mdr1,
dhps, k13, and cytb
'''

year=snakemake.params.year
cleaned_coverage=snakemake.input.cleaned_coverage
cleaned_reference=snakemake.input.cleaned_reference
cleaned_alternate=snakemake.input.cleaned_alternate
genes=snakemake.params.genes

def filter_genes(input_table, output_table):
	for line_number, line in enumerate(open(input_table)):
		line=line.strip().split()
		if line_number==x:
			for column_number, column in enumerate(line)