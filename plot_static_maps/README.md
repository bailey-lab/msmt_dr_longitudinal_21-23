# Plot Static Maps
This step uses prevalences, shape files, and manually labeled features of Kagera
to produce high quality svg files for the Kagera region - this pipeline is the
backbone of most of the figures associated with this study.

# Inputs and Outputs
Inputs include:
 - the prevalences input tables generated during AA_table_visualization, e.g. 
   AA_table_visualization/2021_3_1_output/prevalence_summary.tsv
 - R_scripts/kagera_only_all_muts.R - the central R script that produces graphs
 - python_scripts/prevalence_to_R.py - converts prevalence tables to a format suitable for R
 - envs/R_environment.yaml - A list of packages needed by R
 - input_shape_files - Shape files that describe Tanzania and Kagera maps
Outputs include:
 - R_tables - a temporary file with mutations to graph reformatted by the python
   script into a format compatible with the R script.
 - kagera_graphs_year_cov_alt - a folder with png and svg plots for all
   mutations that were in the prevalence_summary for a given year, coverage
   threshold, and alternate threshold.

# Running the pipeline
This pipeline is somewhat manual. This is because the get_muts function in
kagera_plotting.smk wants a list of all mutations, but I'm not sure how to
modify the rule all so that 'all mutations' is defined differently for each
year and threshold. Because of this, I also didn't bother to send the converted
R tables to different folders so they will overwrite each other if multiple
years/thresholds are run in parallel.

 - Open the kagera_plotting.smk file and make sure the prevalence_file, year,
   and threshold steps are all set correctly for the dataset you want to graph.
 - activate a snakemake environment.
 - Run the pipeline with snakemake -s kagera_plotting.smk --use-conda
 - re-run the pipeline with different years, thresholds, and prevalence files.