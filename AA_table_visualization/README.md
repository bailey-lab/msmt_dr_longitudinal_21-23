# AA table visualization
This pipeline is for filtering, summarizing, and visualizing amino acid changes
that have occurred across a sequencing dataset, and takes amino acid tables as
input. The documentation for this pipeline can be found here:
 https://github.com/simkin-bioinformatics/AA_table_visualization

## Inputs and outputs
I've run the pipeline for 2021, 2022, and 2023 datasets. Inputs include:
 - cleaned_metadata: copied from clean_kagera/kagera_stats_v6/cleaned_metadata
 - longitudinal_AA_tables: copied from clean_kagera/kagera_stats_v6/cleaned_
   filtered_AA_tables
 - variant_graphing.ipynb - copied from AA_table_visualization. Variables are
   filled out to match the 2021 3_1 dataset but can easily be substituted with
   other years and filtering thresholds

Outputs are generated for each year, in the format year_coverage_alternate_
output, where coverage and alternate are thresholds used to consider a mutation
'present' within a sample (if coverage and alternate UMI counts meet or exceed the
respective threshold numbers)', 'covered' (if coverage UMI count meets the
threshold but alternate UMI count does not) or 'not covered' (if neither
coverage nor alternate UMI counts meet the thresholds). Outputs in each 'year'
folder include:
 - prevalences_input_table.csv: This is an intermediate output in which the
 counts from the original AA table are converted into 1's if a sample/mutation
 passes filtering, and 0's if the sample/mutation does not pass filtering.
 - prevalence_summary.tsv: this is a summary of the prevalence of each mutation
 within each district.
 - within_sample_allele_frequencies.csv: these are the number of UMIs that
 contain the mutation within a sample divided by the number of total UMIs within
 the sample that 'cover' the genomic region.
 - k13_Arg561His_District.html: this is an example interactive plot of the data,
   for the Arg561His mutation of kelch13

