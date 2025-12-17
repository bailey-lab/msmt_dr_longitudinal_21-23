# clean_kagera
The analysis for this manuscript originally included metadata, AA tables, and
fastq files that were not part of the final manuscript. This is because
sequencing is not yet complete for other regions of Tanzania outside of Kagera.
To create a clean dataset for manuscript deposition, we added some cleaning
scripts, which are included here:

clean_kagera_part1.smk was used to compare metadata, AA tables (output by
miptools), and original fastq files against each other. Most of the outputs are
stored in kagera_stats_v6. Stats on the number of "matching" samples from these
3 sources are located in "sample_comparison_stats". Any "odd" fastq files or AA
tables were added to odd_AA_tables and odd_fastq_files. sample_dicts holds
sample names of the 3 sources. A user is meant to manually examine any "odd"
entries and fix them up (sending the output to manual_corrections, as included
in this repo).

clean_kagera_part2.smk was used to perform all of the corrections from
manual_corrections. clean_kagera_part2.smk also includes a step that prints any
samples (after corrections were complete) that don't fit canonical formatting
expectations, which revealed a few more samples that were added to
manual_corrections for re-running of part2. Finally, because the original amino
acid tables included over 3,000 columns, the script trims down the amino acid
tables to include only the known drug resistance genes that are in this
manuscript. The final outputs are:
 - kagera_stats_v6/cleaned_filtered_AA_tables
 - kagera_stats_v6/fastq_filenames
 - kagera_stats_v6/cleaned_metadata

For purposes of the manuscript, "attempted" samples are those found in the
tables of cleaned_metadata, and represent samples that have both fastq files and
metadata available after applying all rounds of manual corrections. "Successful"
samples are those that have fastq files, metadata, and output AA tables. Only
successful samples have been deposited in sra.
