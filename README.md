P2 TMT-MS: Multi-group Processing Overview
The process_16plex_tmt_multiple_groups.py script cleans IP2-generated TMT quant files and performs statistical analysis for experiments with more than two groups.
Features

Drops columns with headers containing "avg" (case-insensitive)
Removes contaminants, reverse hits, and low spectral count rows
Adds a leftmost GeneName column
Renames TMT intensity columns using sample metadata
Performs one-way ANOVA across groups
Computes post hoc p-values for each group pair
Adds ratio columns for each post hoc comparison
Excludes TMT channels not listed in metadata
Writes a cleaned, annotated results CSV

Inputs
Place both files in the same folder:

IP2 quant export
Example name: proList.csv (may be LMC.csv, LGS.csv, etc.)


Sample metadata
Filename: Sample_info.csv
Required columns:
TMT_channels
Sample_ids
Group
Sex





Example Sample_info.csv
TMT_channels,Sample_ids,Group,Sex
126,Mouse01,A,M
127N,Mouse02,A,F
127C,Mouse03,B,M
128N,Mouse04,B,F
128C,Mouse05,C,M
129N,Mouse06,C,FPens
129C,Mouse07,D,M
130N,Mouse08,D,F


Groups are defined in the Group column (e.g., A, B, C, D). Number and names may vary.

Command Line Usage
Basic Example
python process_16plex_tmt_multiple_groups.py --prolist LMC.csv --sampleinfo Sample_info.csv -o LMC_results.csv --spec-count-min 3

Arguments

--prolist: Path to the IP2 TMT results CSV (e.g., proList.csv)
--sampleinfo: Path to Sample_info.csv
-o: Output CSV filename
--spec-count-min: Minimum value for the spectral count filter

Processing Steps
Column Cleanup

Remove every column with "avg" in the header (case-insensitive).

Row Filters

Drop rows where:
Accession contains contaminant_
Accession contains Reverse_
Spectral count is less than --spec-count-min



Column Additions and Renaming

Insert a new leftmost column named GeneName
Rename TMT intensity columns to Sample_<Sample_ids> (e.g., Sample_Mouse01_A_M)

Channel Selection

Keep only TMT channels listed in Sample_info.csv

Statistics

Run one-way ANOVA across all groups using renamed intensity columns
Compute post hoc p-values for each group pair (e.g., A vs B, A vs C)
Add ratio columns for each pairwise comparison using mean(Group1) / mean(Group2)

Output

Save the processed table to the specified output CSV

Example Output Columns

GeneName
Original annotation columns from IP2 (excluding avg* columns)
Renamed sample intensity columns (e.g., Sample_Mouse01_A_M)
anova_pvalue
Pairwise comparison columns (e.g., A_vs_B_pvalue, A_vs_B_ratio)
Additional pairwise p-values and ratios for all group combinations

Notes

Ensure TMT_channels values in Sample_info.csv match TMT channel headers in the IP2 file.
Only samples listed in Sample_info.csv are analyzed.
Use a higher --spec-count-min for stricter filtering.

Quick Start
In the folder containing LMC.csv and Sample_info.csv:
python process_16plex_tmt_multiple_groups.py --prolist LMC.csv --sampleinfo Sample_info.csv -o LMC_results.csv --spec-count-min 3

Example files are included to illustrate expected formats.
