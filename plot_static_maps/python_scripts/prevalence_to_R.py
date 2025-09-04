'''
converts a prevalence table (as output by miptools or other programs) into a
table of mutations suitable for R analysis
'''
prevalence_table=snakemake.input.prevalence_table
output_file=open(snakemake.output.R_muts, 'w')
longitude_min=snakemake.params.longitude_min
longitude_max=snakemake.params.longitude_max
latitude_min=snakemake.params.latitude_min
latitude_max=snakemake.params.latitude_max

database={}
for row_number, row in enumerate(open(prevalence_table)):
	row=row.strip().split('\t')
	if row_number==0:
		mut_lookup={}
		mutations=[]
		location_type=row[0]
		for column_number, column in enumerate(row):
			if column_number>2:
				column=column.replace('-', '_')
				mut_lookup[column_number]=column
				mutations.append(column)
	elif row[0]!='overall':
		location=row[0]
		latitude=row[1]
		longitude=row[2]
		database.setdefault(location, {})
		database[location]['Latitude']=float(latitude)
		database[location]['Longitude']=float(longitude)
		for column_number, column in enumerate(row):
			if column_number>2:
				mutation=mut_lookup[column_number]
				print('column is', column)
				percent, frac=column.split(' ')
				percent=str(round(float(percent)*100, 1))
				size=str(int(frac[1:-1].split('/')[1]))
				database[location][mutation]=[percent, size]

header=[location_type, 'Latitude', 'Longitude']
for mutation in mutations:
	header.append(mutation+'_percent')
	header.append(mutation+'_size')
output_file.write(','.join(header)+'\n')

for location in database:
	latitude=database[location]['Latitude']
	longitude=database[location]['Longitude']
	if latitude>latitude_min and latitude<latitude_max and longitude>longitude_min and longitude<longitude_max:
		output_line=[location, str(latitude), str(longitude)]
		for mutation in mutations:
			output_line.append(database[location][mutation][0])
			output_line.append(database[location][mutation][1])
		output_file.write(','.join(output_line)+'\n')