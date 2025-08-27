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

## Step 2: Calculate Prevalences and Plot Dynamic Maps
This step takes AA tables and metadata from step 1 and calculates the fraction
of samples that have each mutation at a given geographic location (in our case
a district of Tanzania, which is a level 2 identifier, akin to a county within
the USA). This step calculates amino acid prevalences in the process, and adds
dynamic maps in plotly. These maps are useful for interactive scrolling but less
useful for publication quality images. Uses this github repo as a template:
https://github.com/simkin-bioinformatics/AA_table_visualization
A more detailed description (along with copies of input scripts) can be
found in the AA_table_visualization folder.

## Step 3: Plot Static Maps
This step takes uses the amino acid prevalences from step 2 to generate high
quality static maps in R. A more detailed description (along with copies of
input scripts) can be found in the plot_static_maps folder.

## Step 4: Plot noncanonical K13
This step is similar to step 3 in that it uses the amino acid prevalences from
step 2 to generate high quality static maps in R. However, instead of treating
all of the mutations in the propeller domain of k13 as separate mutations, it
searches for novel k13 propeller domain mutations. This is accomplished by
creating a slightly modified AA table prior to generating static maps in R. A
more detailed description and script can be found in the plot_noncanonical_k13
folder.

## Step 5: generate confidence intervals
This is a very small step and adds confidence intervals to the prevalence tables
generated in step 2. A more detailed description (along with input and output
files) can be found in the AA_table_visualization folder.
