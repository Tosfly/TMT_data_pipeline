# TMT_data_pipeline
Easy handle IP2 generated IP2 quant files.
1. process_16plex_tmt_multiple_groups.py is to handle data cleanup and calculate ANOVA and post hoc p-values for TMT-MS with more than two groups.
2. input files: TMT result in proList.csv (download from IP2 quant). Sample information in Sample_info.csv (mannually create, columns: TMT_channels, Sample_ids, Group, Sex). Put them in the same folder.
3. run this line: python process_16plex_tmt_multiple_groups.py --prolist proList.csv --sampleinfo Sample_info.csv
4. Please find the example files for the format.
5. The output file:
   a. Remove all columns whose header contains "avg".
   b. Remove all rows where the column "accession" contains "contaminant_".
   c. Remove all rows where the column "accession" contains "Reverse_".
   d. Add a new column called "GeneName" at the leftmost position.
   e. Use Sample_info.csv to rename the TMT intensity columns. in this format: Sample_Sample_ids_Group_Sex.
   f. Calculate p-values across all the groups using only columns with the renamed format. Groups are defined by the "Group" column in 
      Sample_info.csv (in this case four groups: A, B, C, D. The group number may be different but will always names in a similar way).
   g. Add new columns to show post hoc p-values between each pair of groups (e.g., A vs B, A vs C, etc.).
   h. Save the final processed data to a new CSV file.




