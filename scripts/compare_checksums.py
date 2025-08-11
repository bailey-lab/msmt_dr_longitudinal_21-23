'''
code block from google gemini
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

file1_path = "/nfs/jbailey5/baileyweb/asimkin/msmt_projects/msmt_longitudinal_DR/msmt_dr_longitudinal_21-23/2021/kagera_districts/DR_analyses.yaml"
file2_path = "/nfs/jbailey5/baileyweb/asimkin/msmt_projects/msmt_longitudinal_DR/msmt_dr_longitudinal_21-23/2021/kagera_districts/DR_analyses.yaml"
file_checksum = calculate_file_sha256(file_path)