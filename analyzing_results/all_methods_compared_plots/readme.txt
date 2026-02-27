This directory contains the data and code for recreating the "all methods compared" plots in the "Ankh-score produces better sequence alignments than AlphaFold3" manuscript.
The 4 CSV files: 

1. all_exampes_dataframe_BB.csv
2. all_exampes_dataframe_CDD.csv
3. OVER_05_all_exampes_dataframe_BB.csv
4. OVER_05_all_exampes_dataframe_CDD.csv

contain all results from the project for Ankh, US-align, Blosum45, Blosum50, Blosum62, Blosum80, Blosum90. 

Each row corresponds to a specific example, and the columns contain all 6 distances (SSP, SEQ, POS, CC, D, IA), for each method, as well as some identifying information for each example.
OVER_05_all_exampes_dataframe_BB.csv is a subset of  all_exampes_dataframe_BB.csv. It only contains rows/examples that have a min TM-score of 0.5 (score obtained from US-align).
OVER_05_all_exampes_dataframe_CDD.csv is a subset of  all_exampes_dataframe_CDD.csv. It only contains rows/examples that have a min TM-score of 0.5 (score obtained from US-align).

Once all required data is stored in the format of the above CSV files, the "all methods compared" plots can be easily generated with the code in sliding_window_plot.ipynb.