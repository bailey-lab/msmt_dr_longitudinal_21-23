# msmt_dr_longitudinal_21-23
repo for a manuscript describing msmt longitudinal data from 2021-2023

## Step 1: clean metadata, AA tables, and fastq files
The analysis for this manuscript originally included metadata, AA tables, and
fastq files that were not part of the final manuscript. This is because
sequencing is not yet complete for other regions of Tanzania outside of Kagera,
so there is some metadata for samples that haven't been sequenced yet, and some
sequencing that is not ready for publication yet. To create a clean dataset for
manuscript deposition, we added some cleaning scripts, which are described in
the clean_kagera folder. Because these uncleaned input files are not part of
this study, some of the input files for these steps are not included in this
repo.

## Step 2: Plot Dynamic Maps
This step takes AA tables and metadata from step 1 and adds dynamic maps in
plotly. These maps are useful for interactive scrolling but less useful for
publication quality images. Calculates amino acid prevalences in the process.
Uses this github repo as a template:
https://github.com/simkin-bioinformatics/AA_table_visualization
A more detailed description (along with copies of input scripts) can be
found in the plot_dynamic_maps folder.

## Step 3: Plot Static Maps
This step takes uses the amino acid prevalences from step 2 to generate high
quality static maps in R. A more detailed description (along with copies of
input scripts) can be found in the plot_dynamic_maps folder.
