# Complex haplotypes
--------------------

This folder contains some custom code for analyzing complex haplotypes. Briefly,
the pipeline is as follows:

1. Run NFD_plus_GL.py to create specialized AA tables for haplotypes
2. Run variant_graphing.ipynb to convert these AA tables into mutation
   prevalences.
3. Repeat for 2021, 2022, and 2023.
4. Run calculate_confidence_intervals_v2.smk to convert these prevalences into a
   final table with confidence intervals.