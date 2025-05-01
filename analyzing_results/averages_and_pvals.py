import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
import os
import math
import statistics

#this shows how the combined distance file for a model and be compared against the combined distance files from BLOSUM
#the same process is followed to compare the combined distance files for two models
#example
def main():

    model = pd.read_csv("model_file.csv")

    #this can be replaced to point to the combined distance file for another model
    #store all blosum files in a blosums folder
    blosum_path = "folder/with/all/blosum/distance/files"
    blosum_files_old = os.listdir(blosum_path)
    blosum_files = []
    for element in blosum_files_old:
        if element.endswith(".csv"):
            blosum_files.append(element)
    
    #example for 8 CDD domains, change this according to which example distances are being analyzed
    col_names = [ 'cd00012 SSP', 'cd00012 SEQ', 'cd00012 POS', 'cd00012 DD', 'cd00012 CC',
    'cd00083 SSP', 'cd00083 SEQ', 'cd00083 POS', 'cd00083 DD', 'cd00083 CC', 
    'cd00173 SSP', 'cd00173 SEQ', 'cd00173 POS', 'cd00173 DD', 'cd00173 CC',
    'cd00637 SSP', 'cd00637 SEQ', 'cd00637 POS', 'cd00637 DD', 'cd00637 CC',
    'cd01040 SSP', 'cd01040 SEQ', 'cd01040 POS', 'cd01040 DD', 'cd01040 CC', 
    'cd01068 SSP', 'cd01068 SEQ', 'cd01068 POS', 'cd01068 DD', 'cd01068 CC', 
    'cd10140 SSP', 'cd10140 SEQ', 'cd10140 POS', 'cd10140 DD', 'cd10140 CC',
    'cd21116 SSP', 'cd21116 SEQ', 'cd21116 POS', 'cd21116 DD', 'cd21116 CC' ]

    blosum_45 = 0
    blosum_50 = 0
    blosum_62 = 0
    blosum_80 = 0
    blosum_90 = 0
    tie_1 = 0
    m = 0
    b = 0

    current_element = ""
    row = []
    all_rows = []
    all_counts = []
    avg_list = []

    for element in col_names:
        domain = element.split(" ")[0]

        if domain != current_element:
            
            current_element = domain
            row = []
            better_counts = []
            row.append(domain) #element 0
            better_counts.append(domain)
            all_rows.append(row)
            all_counts.append(better_counts)


        model_column_values = model[element].tolist()
        best_averages = math.inf
        best_blosum = ""
        best_col_vals = []

        #this is to handle comparing multiple blosum files at once, but still works if only 1 (model) file is used 
        for blosum in blosum_files:
            #no outlier
            blousm_max = pd.read_csv(str(blosum_path)+"/"+str(blosum))
            blousm_column_values =  blousm_max[element].tolist() #migth still have outliers
            new_average  = mean_with_no_outliers(blousm_max, element)
            if new_average < best_averages:
                best_averages = new_average
                best_blosum = blosum
                best_col_vals = blousm_column_values
        
        #no outliers average for model
        model_average = mean_with_no_outliers(model, element)
        row.append(round(model_average, 6)) #element 1
        row.append(round(best_averages,6)) #element 2
        row.append(best_blosum.split(".csv")[0].split("_")[1])

        #again, if only 1 file is provided, still works
        if best_blosum == "blosum_45.csv":
            blosum_45 = blosum_45+1
            
        if best_blosum == "blosum_50.csv":
            blosum_50 = blosum_50+1
            
        if best_blosum == "blosum_62.csv":
            blosum_62 = blosum_62+1

        if best_blosum == "blosum_80.csv":
            blosum_80 = blosum_80+1
            
        if best_blosum == "blosum_90.csv":
            blosum_90 = blosum_90+1
        

        #counting how many times model is better
        if np.nanmean(model_column_values) == 0:
            model_column_values[0] = 0.0000001 #so the p-val function doesnt break
            
        better_model, better_blosum, tie = compare_two_lists(model_column_values,best_col_vals)
        m = m + better_model
        b = b + better_blosum
        tie_1 = tie_1 + tie

        #get stats
        better_counts.append(better_model)
        better_counts.append(better_blosum)
        better_counts.append(tie)
        avg = better_model/100
        avg_list.append(avg)

            
        if model_column_values == best_col_vals:
            statistic, p_value = 1, 1
            rounded_value = f"{p_value:.4e}"
            row.append(rounded_value)
        else:
            statistic, p_value = wilcoxon_ignore_nan(model_column_values, best_col_vals)
            if (statistic == None):
                continue
            rounded_value = f"{p_value:.4e}"
            row.append(rounded_value)



    # use to above info to make a CSV file
    #can change if using something other than ProtT5 or BLOSUM
    col_names_all = ["Domain", "CC ProtT5 Average", "CC BLOSUM Average", "CC Best BLOSUM", "CC P-Val", "DD ProtT5 Average", "DD BLOSUM Average", "DD Best BLOSUM", "DD P-Val", 
                    "POS ProtT5 Average", "POS BLOSUM Average", "POS Best BLOSUM", "POS P-Val", "SEQ ProtT5 Average", "SEQ BLOSUM Average", "SEQ Best BLOSUM", "SEQ P-Val",
                    "SSP ProtT5 Average", "SSP BLOSUM Average", "SSP Best BLOSUM", "SSP P-Val"]
    df = pd.DataFrame(columns=col_names_all)

    for element in all_rows:
        #again, can swap out ProtT5 and Blosum for other methods
        if len(element) == 21: #how many items each element should have
            new_row = {
            "Domain": element[0],
            "CC ProtT5 Average": element[17],
            "CC BLOSUM Average": element[18],
            "CC Best BLOSUM": element[19],
            "CC P-Val": element[20],
            "DD ProtT5 Average": element[13],
            "DD BLOSUM Average": element[14],
            "DD Best BLOSUM": element[15],
            "DD P-Val": element[16],
            "POS ProtT5 Average": element[9],
            "POS BLOSUM Average": element[10],
            "POS Best BLOSUM": element[11],
            "POS P-Val": element[12],
            "SEQ ProtT5 Average": element[5],
            "SEQ BLOSUM Average": element[6],
            "SEQ Best BLOSUM": element[7],
            "SEQ P-Val": element[8],
            "SSP ProtT5 Average": element[1],
            "SSP BLOSUM Average": element[2],
            "SSP Best BLOSUM": element[3],
            "SSP P-Val": element[4]
        }
        new_row_df = pd.DataFrame([new_row])
        df = pd.concat([df, new_row_df], ignore_index=True)

    #get the averages and p-values to compare two methods (in this example it compares 5 blosums and ProtT5)
    df.to_csv('averages_and_p_vals.csv')

    #making the better counts file 
    col_names_all = ["Domain", "CC ProtT5 Better", "CC BLOSUM Better", "CC Equal", "DD ProtT5 Better", "DD BLOSUM Better", "DD Equal", 
                "POS ProtT5 Better", "POS BLOSUM Better", "POS Equal", "SEQ ProtT5 Better", "SEQ BLOSUM Better", "SEQ Equal",
                "SSP ProtT5 Better", "SSP BLOSUM Better", "SSP Equal"]
                
    df2 = pd.DataFrame(columns=col_names_all)
    for element in all_counts:
        if len(element) == 16: #how many items each element should have
            new_row = {
                "Domain": element[0],
                "CC ProtT5 Better" : element[13], 
                "CC BLOSUM Better" : element[14], 
                "CC Equal" : element[15], 
                "DD ProtT5 Better" : element[10], 
                "DD BLOSUM Better" : element[11], 
                "DD Equal" : element[12], 
                "POS ProtT5 Better" : element[7], 
                "POS BLOSUM Better" : element[8], 
                "POS Equal" : element[9], 
                "SEQ ProtT5 Better" : element[4],
                "SEQ BLOSUM Better" : element[5], 
                "SEQ Equal" : element[6],
                "SSP ProtT5 Better" : element[1],
                "SSP BLOSUM Better" : element[2], 
                "SSP Equal" : element[3]
            }
            new_row_df = pd.DataFrame([new_row])
            df2 = pd.concat([df2, new_row_df], ignore_index=True)
    #this compares all tests for 2 methods head to head to see which one is "better"/ has a lower distance
    df2.to_csv('better_counts.csv')
        







#helper functions

def mean_with_no_outliers(data_frame, column):
    #drop NaN values from the column
    data_without_nan = data_frame[column].dropna()
    #calculate IQR for data2
    Q1_data2 = np.percentile(data_without_nan, 25)
    Q3_data2 = np.percentile(data_without_nan, 75)
    IQR_data2 = Q3_data2 - Q1_data2
    lower_bound_data2 = Q1_data2 - 1.5 * IQR_data2
    upper_bound_data2 = Q3_data2 + 1.5 * IQR_data2
    #filter out outliers for data2
    filtered_data2 = data_without_nan[(data_without_nan >= lower_bound_data2) & (data_without_nan <= upper_bound_data2)]
    #calculate mean for filtered data2
    mean_data2 = filtered_data2.mean()
    return mean_data2

def wilcoxon_ignore_nan(model_column_values, best_col_vals):
    #convert inputs to numpy arrays for element-wise operations
    model_column_values = np.array(model_column_values)
    best_col_vals = np.array(best_col_vals)
    #filter out pairs where either value is NaN
    valid_indices = ~np.isnan(model_column_values) & ~np.isnan(best_col_vals)
    filtered_model_values = model_column_values[valid_indices]
    filtered_best_values = best_col_vals[valid_indices]
    if np.array_equal(filtered_model_values, filtered_best_values):
        return None, None  
    #perform the Wilcoxon test
    statistic, p_value = wilcoxon(filtered_model_values, filtered_best_values)
    return statistic, p_value

def compare_two_lists(A,B):
    if len(A) != len(B):
        print("error: lists need to be the same length")
        return
    length_list = len(A)
    i = 0
    better_A = 0
    better_B = 0
    tie = 0
    for i in range(0,length_list):

        if A[i] < B[i]:
            better_A = better_A + 1
        if B[i] < A[i]:
            better_B = better_B + 1
        if A[i] == B[i]:
            tie = tie + 1
    return better_A, better_B, tie


if __name__ == "__main__":
    main()