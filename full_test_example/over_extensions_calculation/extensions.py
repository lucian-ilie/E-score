#move the completed domain_test folder(s) into the over_extension_calculations directory
#make sure the job file has all the domains as arguments
import sys
import re
import shutil
import numpy as np
import pandas as pd

#used to extract info from output files:
#the below is a modifed extract alignment function that now computes how much the alignments over extend
def extract_alignment_mod(protein1, protein2, ref1, ref2, algn1, algn2):
    
    ref1_no_gaps = ref1.replace('-','')
    ref2_no_gaps = ref2.replace('-','')
    algn1_no_gaps = algn1.replace('-','')
    algn2_no_gaps = algn2.replace('-','')

    left_top_extension = 0
    right_top_extension = 0
    left_bottom_extension = 0
    right_bottom_extension = 0

    #getting the needed indices
    #relative to protein 1
    #where is the reference
    ref1_in_protein1_start = protein1.find(ref1_no_gaps)
    ref1_in_protein1_end = ref1_in_protein1_start + len(ref1_no_gaps)

    #where is the algnment 
    algn1_in_protein1_start = protein1.find(algn1_no_gaps)
    algn1_in_protein1_end = algn1_in_protein1_start + len(algn1_no_gaps)

    #relative to protein 2
    #where is the reference
    ref2_in_protein2_start = protein2.find(ref2_no_gaps)
    ref2_in_protein2_end = ref2_in_protein2_start + len(ref2_no_gaps)

    #where is the algnment 
    algn2_in_protein2_start = protein2.find(algn2_no_gaps)
    algn2_in_protein2_end = algn2_in_protein2_start + len(algn2_no_gaps)

    left_top_extension = max(0,ref1_in_protein1_start-algn1_in_protein1_start)
    left_bottom_extension = max(0,ref2_in_protein2_start-algn2_in_protein2_start)
    
    right_top_extension = max(0,algn1_in_protein1_end-ref1_in_protein1_end)
    right_bottom_extension = max(0,algn2_in_protein2_end-ref2_in_protein2_end)
    
    return left_top_extension, right_top_extension, left_bottom_extension, right_bottom_extension

#helper to get a sequence from a list
def get_sequence(seq_list):
    sequence = ""
    for element in seq_list:
        parts = element.split()
        sequence = sequence + str(parts[3:][1])
        sequence = re.sub(r'\d+$', '', sequence) #incase any numbers are added
    return sequence

#we are getting these results based on the output files, to save us from having to re-compute the alignments, we extract the 
#information by parsing the output files 
def get_info_from_file(output_file_name):
    
    #needed as input to the extract alignment function
    protein1 = ""
    protein2 = ""
    ref1 = ""
    ref2 = ""
    algn1 = ""
    algn2 = ""
    #other helper variables 
    BLOSUM_options = ["BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90"]
    ref_1_lines = []
    ref_2_lines = []
    rawAlign_1_lines = []
    rawAlign_2_lines = []
    scoring_method_line = None
    following_lines = []
    ref_done = False #will let us know which section of the file we are in

    #read the file 
    with open(output_file_name, 'r') as file:
        lines = file.readlines()
    
    #parsing the file
    for i, line in enumerate(lines):

        if "Extracted Alignment:" in line:
            ref1 = get_sequence(ref_1_lines)
            ref2 = get_sequence(ref_2_lines)

            algn1 = get_sequence(rawAlign_1_lines)
            algn2 = get_sequence(rawAlign_2_lines)

            return protein1, protein2, ref1, ref2, algn1, algn2

        #geting the proteins:
        if lines[i].startswith("Protein 1:"):
            protein1 = lines[i + 2].strip()
        if lines[i].startswith("Protein 2:"):
            protein2 = lines[i + 2].strip()

        #getting the reference
        if line.startswith("Seq 1 :") and ref_done == False:
            ref_1_lines.append(line.strip())
        if line.startswith("Seq 2 : ") and ref_done == False:
            ref_2_lines.append(line.strip())
        
        #getting the raw alignment
        if line.startswith("Seq 1 :") and ref_done == True:
            rawAlign_1_lines.append(line.strip())
        if line.startswith("Seq 2 :") and ref_done == True:
            rawAlign_2_lines.append(line.strip())

        #entering section section
        if "Scoring Method" in line:
            ref_done = True

            scoring_method_line = line.strip() #store the line
            scoring_method = scoring_method_line.split(":")[1].lstrip() #keep only the part after :, and remove leading space
            #checking if the storing method is a BLOSUM matrix:
            if scoring_method in BLOSUM_options:
                following_lines = [lines[j].strip() for j in range(i+1, i+5)]
                #needed parameters:
                shift_down_or_blosum= str("_") + str(scoring_method)
                gap_open = following_lines[1].split(":")[1].lstrip()
                gap_extension = following_lines[2].split(":")[1].lstrip()
            #its a model
            else:
                following_lines = [lines[j].strip() for j in range(i+1, i+6)]
                #needed parameters:
                shift_down_or_blosum = following_lines[3].split(":")[1].lstrip() 
                gap_open = following_lines[1].split(":")[1].lstrip()
                gap_extension = following_lines[2].split(":")[1].lstrip()





### BLOSUM CSV FILE FOR EACH MSA IN MSA LIST

def get_blosum_rows(MSA,results_folder_path):
    rows = []
    #given the MSA, go into the MSA_test_1 folder
    #save the second part of all the names in the folder after "_test_1"
    items = os.listdir(str(results_folder_path) + "/" +str(MSA)+"_test_1")
    for element in items:
        rows.append(element.split("_test1_")[1])
    return rows


def blosum_csv_distances(msa,folder_old_data, folder_new_results):

    rows = get_blosum_rows(msa,folder_old_data)  #blosum45-90
    tests = [] #tests 1 - 100

    items = os.listdir(str(folder_old_data))
    for element in items:
        if "test" in element:
            tests.append(element)
    
    #FOR EACH BLOSUM
    for row in rows: #get all the tests for each blosum
        column_0 = row.split("_")[0]
        gap_open = row.split("_")[2]
        gap_extension = row.split("_")[3].split(".")[0]

        columns = ["test", "left top extension", "right top extension", "left bottom extension", "right bottom extension"]
        df = pd.DataFrame(columns=columns)
        file_name = str(folder_new_results) + "/"+ str(column_0) + "_" + str(gap_open) + "_" + str(gap_extension) + "_all_examples.csv"   

        #GET THE RESULTS FROM THE 100 TESTS 
        for test in tests:
            path = str(folder_old_data) + "/" + str(test)
            files = os.listdir(path)
            for element in files: #is is the result of the current blosum
                if column_0 in element:
                    test = element.split("_")[1]
                    path = str(path) + "/" + str(element)
                    protein1, protein2, ref1, ref2, algn1, algn2= get_info_from_file(path)     
                    left_top_extension, right_top_extension, left_bottom_extension, right_bottom_extension = extract_alignment_mod(protein1, protein2, ref1, ref2, algn1, algn2)
                    
                    new_row = pd.DataFrame([{
                        "test": test,
                        "left top extension": left_top_extension,
                        "right top extension": right_top_extension,
                        "left bottom extension": left_bottom_extension,
                        "right bottom extension": right_bottom_extension
                    }])
                    
                    df = pd.concat([df, new_row], ignore_index=True)
                    break

        df.to_csv(file_name, index=False)


def summarize_BLOSUM_for_MSA(MSA,path_over_data_for_this_msa, home_dir):

    columns = ["BLOSUM type", "Gap Open", "Gap Extension", "left top extension avg", "right top extension avg", "left bottom extension avg", "right bottom extension avg"]
    df = pd.DataFrame(columns=columns)
    files = os.listdir(str(path_over_data_for_this_msa))

    for element in files:
        if element.endswith(".csv"):
            data = pd.read_csv(str(path_over_data_for_this_msa) + "/" +str(element))

            new_row = pd.DataFrame([{
                                    "BLOSUM type": element.split("_")[0],
                                    "Gap Open" : element.split("_")[1],
                                    "Gap Extension": element.split("_")[2],
                                    "left top extension avg": round(data["left top extension"].mean(),3), 
                                    "right top extension avg": round(data["right top extension"].mean(),3),
                                    "left bottom extension avg": round(data["left bottom extension"].mean(),3),
                                    "right bottom extension avg": round(data["right bottom extension"].mean(),3),
                                }])
            df = pd.concat([df, new_row], ignore_index=True)

    df = df.sort_values(by='BLOSUM type', ascending=True)
    file_path = str(home_dir) + "/" + str(MSA) + "_summary.csv"
    
    df.to_csv(file_path, index=False)



def get_model_rows(MSA,results_folder_path):
    rows = []
    #given the MSA, go into the MSA_test_1 folder
    #save the second part of all the names in the folder after "_test_1"
    items = os.listdir(str(results_folder_path) + "/" +str(MSA)+"_test_1")
    for element in items:
        rows.append(element)
    return rows


#model functions

def model_csv_distances(msa,folder_old_data, folder_new_results):
    
    rows = get_model_rows(msa,folder_old_data)  
    print(rows)
    tests = [] #tests 1 - 100

    items = os.listdir(str(folder_old_data))
    for element in items:
        if "test" in element:
            tests.append(element)
    
    print(tests)
    print("\n")
    
    #FOR EACH PARAMETER SPACE
    for row in rows: #get all the tests for parameter set
        column_0 = row.split("_")[1] 
        gap_open = row.split("_")[4]
        gap_extension = row.split("_")[5].split(".txt")[0]

        columns = ["test", "left top extension", "right top extension", "left bottom extension", "right bottom extension"]
        df = pd.DataFrame(columns=columns)
        file_name = str(folder_new_results) + "/"+ str(column_0) + "_" + str(gap_open) + "_" + str(gap_extension) + "_all_examples.csv"   
     
     # GET THE RESULTS FROM THE 100 TESTS 
        for test in tests:

      
            path = str(folder_old_data) + "/" + str(test)
            files = os.listdir(path)
            print(files)
            print(element)
            for element in files: #is is the result of the current blosum
                if row.split("ProtT5_0.00000_")[1] in element:
                    test = element.split("_")[1]
                    path = str(path) + "/" + str(element)
                    protein1, protein2, ref1, ref2, algn1, algn2= get_info_from_file(path)     
                    left_top_extension, right_top_extension, left_bottom_extension, right_bottom_extension = extract_alignment_mod(protein1, protein2, ref1, ref2, algn1, algn2)
                    
                    new_row = pd.DataFrame([{
                        "test": test,
                        "left top extension": left_top_extension,
                        "right top extension": right_top_extension,
                        "left bottom extension": left_bottom_extension,
                        "right bottom extension": right_bottom_extension
                    }])
                    
                    df = pd.concat([df, new_row], ignore_index=True)
                    break

        df.to_csv(file_name, index=False)
    
def summarize_model_for_MSA(MSA,path_over_data_for_this_msa, home_dir):

    columns = ["Gap Open", "Gap Extension", "left top extension avg", "right top extension avg", "left bottom extension avg", "right bottom extension avg"]
    df = pd.DataFrame(columns=columns)
    files = os.listdir(str(path_over_data_for_this_msa))

    for element in files:
        if element.endswith(".csv"):
            data = pd.read_csv(str(path_over_data_for_this_msa) + "/" +str(element))

            new_row = pd.DataFrame([{
                                    "Gap Open" : element.split("_")[1],
                                    "Gap Extension": element.split("_")[2],
                                    "left top extension avg": round(data["left top extension"].mean(),3), 
                                    "right top extension avg": round(data["right top extension"].mean(),3),
                                    "left bottom extension avg": round(data["left bottom extension"].mean(),3),
                                    "right bottom extension avg": round(data["right bottom extension"].mean(),3),
                                }])
        
            df = pd.concat([df, new_row], ignore_index=True)

    df = df.sort_values(by='Gap Open', ascending=False)
    file_path = str(home_dir) + "/" + str(MSA) + "_summary.csv"
    
    df.to_csv(file_path, index=False)


    


def final_summary_msa(MSA_list, home_dir):

    if len(MSA_list) >= 2:
        first_data = pd.read_csv(str(home_dir) + "/"+str(MSA_list[0])+"_summary.csv")
        first_data = first_data.rename(columns={'left top extension avg' : str(MSA_list[0]) + ' left top extension avg',  'right top extension avg' : str(MSA_list[0]) +  ' right top extension avg', 'left bottom extension avg' :  str(MSA_list[0]) + ' left bottom extension avg', 'right bottom extension avg' :  str(MSA_list[0]) + ' right bottom extension avg'})
        
        for i in range(1,len(MSA_list)):
            data = pd.read_csv(str(home_dir) + "/"+str(MSA_list[i])+"_summary.csv")
            #drop first 3 columns
            data = data.drop(data.columns[:3], axis=1)
            data = data.rename(columns={'left top extension avg' : str(MSA_list[i]) + ' left top extension avg',  'right top extension avg' : str(MSA_list[i]) + ' right top extension avg', 'left bottom extension avg' :  str(MSA_list[i]) + ' left bottom extension avg', 'right bottom extension avg' :  str(MSA_list[i]) + ' right bottom extension avg'})

            first_data = pd.concat([first_data, data], axis=1)
            
        file_path = str(home_dir) + "/summary.csv"
        first_data = first_data.sort_values(by='Gap Open', ascending=False)
        first_data.to_csv(file_path, index=False)

    else:
        print("only one msa passed")











def final_summary_model(MSA_list, home_dir):

    if len(MSA_list) >= 2:
        first_data = pd.read_csv(str(home_dir) + "/"+str(MSA_list[0])+"_summary.csv")
        first_data = first_data.rename(columns={'left top extension avg' : str(MSA_list[0]) + ' left top extension avg',  'right top extension avg' : str(MSA_list[0]) +  ' right top extension avg', 'left bottom extension avg' :  str(MSA_list[0]) + ' left bottom extension avg', 'right bottom extension avg' :  str(MSA_list[0]) + ' right bottom extension avg'})
        
        for i in range(1,len(MSA_list)):
            data = pd.read_csv(str(home_dir) + "/"+str(MSA_list[i])+"_summary.csv")
            #drop first 3 columns
            data = data.drop(data.columns[:2], axis=1)
            data = data.rename(columns={'left top extension avg' : str(MSA_list[i]) + ' left top extension avg',  'right top extension avg' : str(MSA_list[i]) + ' right top extension avg', 'left bottom extension avg' :  str(MSA_list[i]) + ' left bottom extension avg', 'right bottom extension avg' :  str(MSA_list[i]) + ' right bottom extension avg'})

            first_data = pd.concat([first_data, data], axis=1)
            
        file_path = str(home_dir) + "/summary.csv"
        first_data = first_data.sort_values(by='Gap Open', ascending=False)
        first_data.to_csv(file_path, index=False)

    else:
        print("only one msa passed")




def main():

    home_dir = os.getcwd()

    MSA_list = str(sys.argv[1])
    MSA_list = MSA_list.split(" ")
    method = str(sys.argv[2])

    for element in MSA_list:

        #for each msa creating a folder to store the extension values for each test
        #in MSA_extensions, there will be 5 files, one for each blosum file storing all the extensios for that blosum
        folder_new_results = str(home_dir) + "/" + str(element) + "_extensions"
        folder_old_data = str(home_dir) + "/" + str(element) + "_tests" 
        
        if os.path.exists(folder_new_results):
            shutil.rmtree(folder_new_results)  #deletes the folder and its contents
        os.makedirs(folder_new_results)  #creates the folder

        if method == "BLOSUM":
            blosum_csv_distances(element,folder_old_data, folder_new_results)
        
            #get the averages for this msa
            summarize_BLOSUM_for_MSA(element,folder_new_results, home_dir)
        
        if method == "MODEL":
            model_csv_distances(element,folder_old_data, folder_new_results)

            #get the averages for this msa
            summarize_model_for_MSA(element,folder_new_results, home_dir)
    
    if method == "BLOSUM":
        final_summary_msa(MSA_list, home_dir)
    
    if method == "MODEL":

        #merge the above files into one document
        final_summary_model(MSA_list, home_dir)


if __name__ == "__main__":
    main()

