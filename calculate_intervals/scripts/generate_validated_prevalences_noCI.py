'''
examines only validated or candidate k13 mutations and reformats them as a table
with three columns to match the 3 prevalence values in table 1 of the study.

This version strips out confidence intervals (for use with main tables)
'''

output_path=snakemake.output.validated_prevalences
CI_folder=snakemake.params.interval_folder
candidate_validated=['k13-Pro441Leu', 'k13-Phe446Ile', 'k13-Gly449Ala',
'k13-Asn458Tyr', 'k13-Cys469Phe', 'k13-Cys469Tyr', 'k13-Met476Ile',
'k13-Ala481Val', 'k13-Tyr493His', 'k13-Arg515Lys', 'k13-Pro527His',
'k13-Asn537Asp', 'k13-Asn537Ile', 'k13-Gly538Val', 'k13-Arg539Thr',
'k13-Ile543Thr', 'k13-Pro553Leu', 'k13-Arg561His', 'k13-Val568Gly',
'k13-Pro574Leu', 'k13-Cys580Tyr', 'k13-Arg622Ile', 'k13-Ala675Val']

all_districts=['Biharamulo', 'Bukoba DC', 'Karagwe', 'Kyerwa', 'Misenyi', 'Muleba', 'Ngara']

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
							prev_dict.setdefault(mut_name, {})
							prev_dict[mut_name].setdefault(year, {})
							prev_dict[mut_name][year].setdefault(district, ['','',''])
							index=line_number%3-1
							prev_dict[mut_name][year][district][index]=value
	all_muts=special_sort(all_muts)
	return prev_dict, all_muts

def print_table(prev_dict, all_muts, output_path):
	#output_line.extend(['0.0/unk', '0.0%', '0.0%-0.0%'])
	output_file=open(output_path, 'w')
	header_line=['', '']
	for district in all_districts:
		header_line.extend([district])
	output_file.write('\t'.join(header_line)+'\n')
	header_line=['mutation', 'year']
	header_line.extend(['N (prev)']*len(all_districts))
	output_file.write('\t'.join(header_line)+'\n')
	for mut_name in all_muts:
		for year in prev_dict[mut_name]:
			output_line=[mut_name, year]
			for district in all_districts:
				if district in prev_dict[mut_name][year]:
					frac, prev, CI=prev_dict[mut_name][year][district]
					frac_prev=f'{frac} ({prev})'
					#prev_CI=f'{prev} ({CI})'
					output_line.append(frac_prev)
				else:
					output_line.append('-')
			output_file.write('\t'.join(output_line)+'\n')

prev_dict, all_muts=populate_dict(CI_folder)
print('prev dict is', prev_dict)
print_table(prev_dict, all_muts, output_path)