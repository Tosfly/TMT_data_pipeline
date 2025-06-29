# TMT_data_pipeline
Easy handle IP2 generated IP2 quant files.
1. process_16plex_tmt_multiple_groups.py is to handle data cleanup and calculate ANOVA and post hoc p-values for TMT-MS with more than two groups.

2. input files: TMT result in proList.csv (download from IP2 quant). Sample information in Sample_info.csv (mannually create). Put them in the same folder.

3. run this line: python process_16plex_tmt_multiple_groups.py --prolist proList.csv --sampleinfo Sample_info.csv
