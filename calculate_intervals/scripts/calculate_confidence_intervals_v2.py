'''
converts a prevalence table (of the type output by miptools) into confidence
intervals (of the type used by MSMT longitudinal papers).

This version uses a wilson confidence interval instead of a normal confidence
interval.
'''
prevalence_file=snakemake.input.prevalence_table
output_table=open(snakemake.output.confidence_intervals, 'w')

def format_line(location, latitude, longitude, mut_values, outuput_table):
	from statsmodels.stats.proportion import proportion_confint
	output_lines=[[location, latitude, longitude],
	[location, latitude, longitude],
	[location, latitude, longitude]]
	for mut_value in mut_values:
		alt_count, cov_count=map(int, mut_value.split(' ')[1][1:-1].split('/'))
		if cov_count>0:
			lower_bound, upper_bound=proportion_confint(count=alt_count, nobs=cov_count, method='beta')
			lower_bound, upper_bound=round(lower_bound*100, 1), round(upper_bound*100, 1)
			output_lines[0].append(f'{alt_count}/{cov_count}')
			output_lines[1].append(f'{round(alt_count/cov_count*100, 1)}')
			output_lines[2].append(f'{lower_bound}-{upper_bound}')
		else:
			output_lines[0].append(f'no coverage')
			output_lines[1].append(f'no coverage')
			output_lines[2].append(f'no coverage')
	for output_line in output_lines:
		output_table.write('\t'.join(output_line)+'\n')

for line_number, line in enumerate(open(prevalence_file)):
	if line_number==0:
		output_table.write(line)
	else:
		line=line.strip().split('\t')
		location_name=line[0]
		latitude=line[1]
		longitude=line[2]
		mut_values=line[3:]
		format_line(location_name, latitude, longitude, mut_values, output_table)