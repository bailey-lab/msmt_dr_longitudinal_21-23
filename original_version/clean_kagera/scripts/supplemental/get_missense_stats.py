'''
Current AA tables are too big to upload to github. This is a simple script to
count the number of simple missense mutations in the AA tables vs. the number of
more complex mutations (to see if limiting to missense mutations would be enough
to upload these).
'''
gene_dict={}
coverage_sums={}
for line_number, line in enumerate(open('/home/alfred/msmt_re_analysis_with_cleaned_metadata/v3_08-11-25_official_ms_github/msmt_dr_longitudinal_21-23/AA_tables/2021_AA_tables/alternate_AA_table.csv')):
	line=line.strip().split(',')
	if line_number==1:
		genes=line
	if line_number==3:
		muts=line
	elif line_number>5:
		counts=line
		for count_number, count in enumerate(counts):
			if count_number>0:
				coverage_sums.setdefault(count_number, 0)
				coverage_sums[count_number]+=int(float(count))

count=0
for mut in muts:
	if mut=='missense_variant':
		count+=1
for gene in genes:
	gene_dict[gene]=gene_dict.setdefault(gene, 0)+1
print('missense is', count)
print('total is', len(muts))
print('genes are', gene_dict)

covered_muts=0
for column in coverage_sums:
	if coverage_sums[column]>0:
		covered_muts+=1
print('mutations with non-zero counts are', covered_muts)