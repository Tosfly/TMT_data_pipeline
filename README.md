IP2 TMT-MS: Multi-group Processing
Overview

process_16plex_tmt_multiple_groups.py cleans IP2-generated TMT quant files and runs statistics for experiments with more than two groups. It:

Drops columns whose headers contain avg

Removes contaminants, reverse hits, and low spectral count rows

Adds a leftmost GeneName column

Renames TMT intensity columns using your sample metadata

Performs one-way ANOVA across groups

Computes post hoc p-values for each group pair

Adds ratio columns next to each post hoc comparison

Excludes TMT channels not listed in your metadata

Writes a cleaned, annotated results CSV

Inputs

Place both files in the same folder:

IP2 quant export
Example name: proList.csv (your file may be LMC.csv, LGS.csv, etc.)

Sample metadata
Filename: Sample_info.csv
Required columns:

TMT_channels

Sample_ids

Group

Sex

Example Sample_info.csv:

TMT_channels,Sample_ids,Group,Sex
126,Mouse01,A,M
127N,Mouse02,A,F
127C,Mouse03,B,M
128N,Mouse04,B,F
128C,Mouse05,C,M
129N,Mouse06,C,F
129C,Mouse07,D,M
130N,Mouse08,D,F


Groups are defined by the Group column and can be A, B, C, D, etc. The number and names may vary.

Command line usage

Basic example:

python process_16plex_tmt_multiple_groups.py \
  --prolist LMC.csv \
  --sampleinfo Sample_info.csv \
  -o LMC_results.csv \
  --spec-count-min 3


Arguments

--prolist Path to the IP2 TMT results CSV (for example proList.csv)

--sampleinfo Path to Sample_info.csv

-o Output CSV filename

--spec-count-min Minimum value for the spec count filter

Processing steps

Column cleanup
Remove every column whose header contains avg (case insensitive).

Row filters

Drop rows where accession contains contaminant_

Drop rows where accession contains Reverse_

Drop rows where spec count is less than --spec-count-min

Column additions and renaming

Insert a new leftmost column named GeneName

Rename TMT intensity columns using this format:
Sample_<Sample_ids>_<Group>_<Sex>
Example: Sample_Mouse01_A_M

Channel selection
Keep only TMT channels present in Sample_info.csv.

Statistics

Run one-way ANOVA across all groups using the renamed intensity columns

For each pair of groups (for example A vs B, A vs C, ...), compute post hoc p-values

Add ratio columns next to each pairwise comparison using mean(Group1) / mean(Group2)

Output
Save the fully processed table to the output CSV you specify.

Example output columns

GeneName

Original annotation columns retained from IP2 (minus any avg* columns)

Renamed sample intensity columns like Sample_Mouse01_A_M

anova_pvalue

For each pair, for example A_vs_B_pvalue and A_vs_B_ratio

Additional pairwise p-values and ratios for all group combinations

Notes

Ensure TMT_channels values in Sample_info.csv match the TMT channel headers in your IP2 file.

Only samples listed in Sample_info.csv are analyzed.

Consider a higher --spec-count-min for a stricter filter.

Quick start
# In the folder that contains LMC.csv and Sample_info.csv
python process_16plex_tmt_multiple_groups.py \
  --prolist LMC.csv \
  --sampleinfo Sample_info.csv \
  -o LMC_results.csv \
  --spec-count-min 3


Example files are included to illustrate the expected formats.
