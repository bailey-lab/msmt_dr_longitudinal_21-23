# Calculate Confidence Intervals

This step calculates confidence intervals on the prevalences reported in earlier
steps of the pipeline, specifically the
year_cov_alt_output/prevalence_summary.tsv files found in the
AA_table_visualization folder. The outputs are prevalence tables that also
contain confidence intervals, and are stored in confidence_interval_outputs.
Also outputs an aggregated table across years for candidate and validated k13
mutations, and a similar aggregated table for known drug resistance mutations.
Finally, there is a very hacky (and manual) script (scripts/mdr1_confints.py)
for calculating confidence intervals for the inverse prevalences of mdr1-86Y -
after inverse prevalences are calculated by hand, it outputs confidence
intervals associated with these prevalences.

The main snakemake script is executed with
`bash
snakemake -s calculate_confidence_intervals.smk --cores 4 --use-conda
`
and assumes a hardcoded input folder of ../AA_table_visualization

The more manual mdr1-86Y script is executed with
python3 mdr1_confints.py