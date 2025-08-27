# Calculate Confidence Intervals

This step calculates confidence intervals on the prevalences reported in earlier
steps of the pipeline, specifically the
year_cov_alt_output/prevalence_summary.tsv files found in the
AA_table_visualization folder. The outputs are prevalence tables that also
contain confidence intervals, and are stored in confidence_interval_outputs.

The script is executed with
`bash
snakemake -s calculate_confidence_intervals.smk --cores 4 --use-conda
`
and assumes a hardcoded input folder of ../AA_table_visualization