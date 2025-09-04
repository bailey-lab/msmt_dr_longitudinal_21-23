'''
a little function to facilitate manual lookup and correction of misformatted
sample names - a user can lookup all samples that contain any substring within
any dataset to try and find candidate replacement values.
'''

import pickle
year='2021'
metadata=pickle.load(open(f'../kagera_stats_v5/sample_dicts/metadata_dict_{year}.pkl', 'rb'))
aa=pickle.load(open(f'../kagera_stats_v5/sample_dicts/AA_dict_{year}.pkl', 'rb'))
fastq=pickle.load(open(f'../kagera_stats_v5/sample_dicts/fastq_dict_{year}.pkl', 'rb'))
databases={'aa':aa, 'metadata':metadata, 'fastq':fastq}

def lookup(database, search_term):
	for key in sorted(database.keys()):
		if search_term in key:
			print(key)

menu_choice=''
print('year is', year)
while menu_choice!=('exit'):
	menu_choice=input('What would you like to do? Options are lookup, replace, and exit: ')
	if menu_choice=='lookup':
		database_name=input('what dictionary would you like to use?')
		database=databases[database_name]
		search_term=input('what search term would you like to use?')
		print('results are')
		lookup(database, search_term)