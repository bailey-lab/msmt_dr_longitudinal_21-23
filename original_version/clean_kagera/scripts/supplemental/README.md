# Supplemental scripts
This folder is for scripts that weren't directly used in the pipeline but that
were useful for debugging, troubleshooting, and cleaning the data.

compare_checksums.py was used to make sure my various duplicated input AA tables
that I had on my hard drive were consistent with each other.

lookup_candidates.py was useful for browsing the metadata, AA, and fastq sample
names to make corrections to the fastq files that other parts of the pipeline
had flagged as problematic, and the outputs of this browsing were used to
produce the manual_corrections.

get_missense_stats.py was useful for testing various ways of cutting down the
output amino acid tables to a size that would easily fit on github.