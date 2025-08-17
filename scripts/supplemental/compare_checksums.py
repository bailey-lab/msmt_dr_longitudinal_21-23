'''
code block from google gemini - checks if two files have the same checksums or
not to see if they're identical
'''
import hashlib

def calculate_file_sha256(filepath):
	"""Calculates the SHA-256 checksum of a file."""
	sha256_hasher = hashlib.sha256()
	with open(filepath, 'rb') as f:
		# Read file in chunks to handle large files efficiently
		for chunk in iter(lambda: f.read(4096), b''):
			sha256_hasher.update(chunk)
	return sha256_hasher.hexdigest()

def compare_checksums(filepath1, filepath2):
	same=False
	if calculate_file_sha256(filepath1)==calculate_file_sha256(filepath2):
		same=True
	return same

filepath1 = "/home/alfred/msmt_re_analysis_with_cleaned_metadata/v3_08-11-25_official_ms_github/coverage_AA_table.csv"
filepath2 = "/home/alfred/msmt_re_analysis_with_cleaned_metadata/v3_08-11-25_official_ms_github/msmt_dr_longitudinal_21-23/AA_tables/2023_AA_tables/coverage_AA_table.csv"
same=compare_checksums(filepath1, filepath2)
if same:
	print(f'filepath1 and filepath2 have identical contents')
else:
	print(f'filepath1 and filepath2 are different')