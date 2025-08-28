'''
examines only validated or candidate k13 mutations and reformats them as a table
with three columns to match the 3 prevalence values in table 1 of the study.
'''

output_path=snakemake.output.validated_prevalences
CI_folder=snakemake.params.interval_folder
candidate_validated=['k13-Pro441Leu', 'k13-Phe446Ile', 'k13-Gly449Ala',
'k13-Asn458Tyr', 'k13-Cys469Phe', 'k13-Cys469Tyr', 'k13-Met476Ile',
'k13-Ala481Val', 'k13-Tyr493His', 'k13-Arg515Lys', 'k13-Pro527His',
'k13-Asn537Ile', 'k13-Asn537Asp', 'k13-Gly538Val', 'k13-Arg539Thr',
'k13-Ile543Thr', 'k13-Pro553Leu', 'k13-Arg561His', 'k13-Val568Gly',
'k13-Pro574Leu', 'k13-Cys580Tyr', 'k13-Arg622Ile', 'k13-Ala675Val']

def special_sort(mut_list):
	sorting_list=[]
	for mut in mut_list:
		gene='-'.join(mut.split('-')[:-1])
		pos=int(mut.split('-')[-1][3:-3])
		sorting_list.append([gene, pos, mut])
	return [item[-1] for item in sorted(sorting_list)]

def populate_dict(CI_folder):
	prev_dict, all_muts={}, set([])
	for year in ['2021', '2022', '2023']:
		prevalence_table=CI_folder+f'/{year}_3_1_CIs/{year}_3_1_confidence_intervals.tsv'
		for line_number, line in enumerate(open(prevalence_table)):
			line=line.strip().split('\t')
			if line_number==0:
				header=line
			else:
				for mut_column, value in enumerate(line):
					mut_name=header[mut_column]
					if mut_name in candidate_validated:
						all_muts.add(mut_name)
						district=line[0]
						if district!='overall':
							prev_dict.setdefault(district, {})
							prev_dict[district].setdefault(year, {})
							prev_dict[district][year].setdefault(mut_name, ['','',''])
							index=line_number%3-1
							prev_dict[district][year][mut_name][index]=value
	all_muts=special_sort(all_muts)
	return prev_dict, all_muts

def print_table(prev_dict, all_muts, output_path):
	output_file=open(output_path, 'w')
	header_line=['District', 'year']
	for mut in all_muts:
		header_line.extend([mut]*3)
	output_file.write('\t'.join(header_line)+'\n')
	for district in prev_dict:
		for year in prev_dict[district]:
			output_line=[district, year]
			for mut in all_muts:
				if mut in prev_dict[district][year]:
					output_line.extend(prev_dict[district][year][mut])
				else:
					output_line.extend(['0.0/unk', '0.0%', '0.0%-0.0%'])
			output_file.write('\t'.join(output_line)+'\n')


prev_dict, all_muts=populate_dict(CI_folder)
print('prev dict is', prev_dict)
print_table(prev_dict, all_muts, output_path)