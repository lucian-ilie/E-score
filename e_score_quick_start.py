#1. required imports
import os
import sys
import warnings
import re
from Bio import SeqIO
from scipy import stats
import numpy as np
import pandas as pd
import torch
import sentencepiece
import shutil
from tqdm import tqdm
from transformers import T5Model, T5Tokenizer, T5EncoderModel, AutoTokenizer, AutoModelForSeq2SeqLM, T5ForConditionalGeneration
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu') 


#2. main function to run a single alignment example
def main():
    
    protein_seq1 = "AACD" #any protein sequences of capital letters
    protein_seq2 = "CDDA" #any protein sequence of capital letters
    
    gap_open = -1 #any value can be used
    gap_extension = -0.1  #any value can be used

    alignment_type = 'local' #can be 'local', 'global' or 'semi-global'
    int_shift = 0 #only matters if local alignment is being used
    #this shifts down all the cosine similarity values down

    #code for initalizing other embedding model
    model, tokenizer = Ankh_initialize()


    #computing the alignment:
    if alignment_type == 'local':
        aligned_seq1, aligned_seq2, max_score = affine_local_alignment_model_and_traceback(protein_seq1, protein_seq2, gap_open, gap_extension, model, tokenizer, int_shift)

    elif alignment_type == 'global':
        aligned_seq1, aligned_seq2, max_score = affine_global_alignment_escore_and_traceback(protein_seq1, protein_seq2, gap_open, gap_extension, model, tokenizer)

    elif alignment_type == 'semi-global':
        aligned_seq1, aligned_seq2, max_score = affine_semi_global_alignment_escore_and_traceback(protein_seq1, protein_seq2, gap_open, gap_extension, model, tokenizer)
    
    #printing out the results
    print("the computed alignment: ")
    print(aligned_seq1)
    print(aligned_seq2)

    print("with the score of: " + str(max_score))



#3. model initalization (need to change to personal file path)
def Ankh_initialize():

    ankh_model = T5ForConditionalGeneration.from_pretrained("/change/path/to/model", local_files_only=True)
    ankh_tokenizer = AutoTokenizer.from_pretrained("/change/path/to/model", local_files_only=True)
    ankh_model = ankh_model.eval()
    ankh_model = ankh_model.to(device)

    return ankh_model, ankh_tokenizer


def get_embs_Ankh(seq, ankh_model, tokenizer):
    
    inputs = tokenizer([seq], padding=True, return_tensors="pt") #pytorch format 
    inputs = {key: value.to(device) for key, value in inputs.items()}

    with torch.no_grad(): #no gradients since we are not training
        outputs = ankh_model.encoder(**inputs)
        last_hidden_states = outputs.last_hidden_state 
    
    sequence_embedding = last_hidden_states[0, :-1, :]
    return sequence_embedding



#4. affine local alignment model and traceback with e-score scoring function
def affine_local_alignment_model_and_traceback(seq1, seq2, g_open = -1, g_ext = -0.1, model = "", tokenizer = "", shift = 0):
    
    #1) initialization
    m = len(seq1); n = len(seq2)
    #3 working matrices:
    #A --> best alignment of seq1[1-i] and seq2[1-j] that aligns seq1[i] and seq2[j]
    #B --> best alignment of seq1[1-i] and seq2[1-j] that aligns gap with seq2[j] (gap in seq1)
    #C --> best alignment of seq1[1-i] and seq2[1-j] that aligns seq1[i] with gap (gap in seq2)
    A = np.zeros([m+1,n+1])
    B = np.zeros([m+1,n+1])
    C = np.zeros([m+1,n+1])

    #3 traceback matrices
    A_tb = np.zeros([m+1,n+1])
    B_tb = np.zeros([m+1,n+1])
    C_tb = np.zeros([m+1,n+1])

    #initializing matrices
    A[0,0] = 0 #needed for recursion
    A[1:,0] = -np.inf
    A[0,1:] = -np.inf
    B[0:,0] = -np.inf #-inf along left 
    C[0,0:] = -np.inf #-inf along top 

    #needed for computing embeddings
    Model = model
    Model_tokenizer = tokenizer

    #get Ankh embedded sequences
    emb1 = get_embs_Ankh(seq1, model, tokenizer).cpu().numpy()
    emb2 = get_embs_Ankh(seq2, model, tokenizer).cpu().numpy()
    cos = torch.nn.CosineSimilarity(dim=0)


    #2 ) Dynamic Programing
    for i in range(1, m + 1): #length of first sequence
        for j in range(1, n + 1): #length of second sequence

            sim = cos(torch.tensor(emb1[i - 1] , dtype = torch.float32), torch.tensor(emb2[j - 1] , dtype = torch.float32)).item()
            match_score = sim + shift #shift should be a negative number

            #filling up A
            A[i,j] = max(A[i-1,j-1] + match_score, B[i-1,j-1] + match_score, C[i-1,j-1] + match_score, 0) #need 0 since this is local alignment
            #store traceback info
            max_index = np.argmax([A[i-1,j-1] + match_score, B[i-1,j-1] + match_score, C[i-1,j-1] + match_score, 0])
            A_tb[i,j] = max_index #0 = from A, 1 = from B, 2 = from C, 3 = end

            #filling up B (gap in seq1)
            B[i,j] = max((A[i,j-1]+g_open+g_ext),(B[i,j-1]+g_ext),(C[i,j-1]+g_open+g_ext), 0)
            #store traceback info
            max_index = np.argmax([A[i,j-1] + g_open+g_ext, B[i,j-1] + g_ext, C[i,j-1] + g_open+g_ext, 0])
            B_tb[i,j] = max_index #0 = from A, 1 = from B, 2 = from C, 3 = end

            #filling up C (gap in seq2)
            C[i,j] = max((A[i-1,j]+g_open+g_ext),(B[i-1,j]+g_open+g_ext),(C[i-1,j]+g_ext), 0)
            #store traceback info
            max_index = np.argmax([A[i-1,j] + g_open+g_ext, B[i-1,j] + g_open+g_ext, C[i-1,j] + g_ext, 0])
            C_tb[i,j] = max_index #0 = from A, 1 = from B, 2 = from C, 3 = end


    #3) traceback
    #find max cell across all 3 matrices
    max_score = float('-inf')
    max_pos = None
    for i in range(m + 1):
        for j in range(n + 1):
            for matrix, tb_matrix in [(A, A_tb), (B, B_tb), (C, C_tb)]:
                if matrix[i, j] > max_score:
                    max_score = matrix[i, j]
                    max_pos = (i, j, tb_matrix, matrix)

    i, j, tb_matrix, matrix = max_pos
    aligned_seq1 = []
    aligned_seq2 = []

    while matrix[i, j] > 0:

        #currently in matrix A
        if tb_matrix is A_tb:

            if tb_matrix[i, j] == 0: #from A
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
                matrix = A
                tb_matrix = A_tb

            elif tb_matrix[i, j] == 1: #from A
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
                matrix = B
                tb_matrix = B_tb

            elif tb_matrix[i, j] == 2: #from C
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
                matrix = C
                tb_matrix = C_tb

            elif tb_matrix[i, j] == 3:
                break

        #currently in matrix B (gap in seq1)
        elif tb_matrix is B_tb:

            if tb_matrix[i, j] == 0:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
                matrix = A
                tb_matrix = A_tb

            elif tb_matrix[i, j] == 1:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
                matrix = B
                tb_matrix = B_tb

            elif tb_matrix[i, j] == 2:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
                matrix = C
                tb_matrix = C_tb

            elif tb_matrix[i, j] == 3:
                break

        #currently in matrix C (gap in seq2)
        elif tb_matrix is C_tb:

            if tb_matrix[i, j] == 0:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
                matrix = A
                tb_matrix = A_tb

            elif tb_matrix[i, j] == 1:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
                matrix = B
                tb_matrix = B_tb

            elif tb_matrix[i, j] == 2:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
                matrix = C
                tb_matrix = C_tb

            elif tb_matrix[i, j] == 3:
                break

    #reverse the alignment lists (as it is backwards, and make elements into a string)
    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2, max_score


#5. affine global alignment model and traceback with e-score scoring function
def affine_global_alignment_escore_and_traceback(seq1, seq2, g_open = -0.25, g_ext = -0.01, model = "", tokenizer = ""):
    
    #1) initialization
    m = len(seq1); n = len(seq2)
    #3 working matrices:
    #A --> best alignment of seq1[1-i] and seq2[1-j] that aligns seq1[i] and seq2[j]
    #B --> best alignment of seq1[1-i] and seq2[1-j] that aligns gap with seq2[j] (gap in seq1)
    #C --> best alignment of seq1[1-i] and seq2[1-j] that aligns seq1[i] with gap (gap in seq2)
    A = np.zeros([m+1,n+1])
    B = np.zeros([m+1,n+1])
    C = np.zeros([m+1,n+1])

    #3 traceback matrices
    A_tb = np.zeros([m+1,n+1])
    B_tb = np.zeros([m+1,n+1])
    C_tb = np.zeros([m+1,n+1])

    #initializing matrices
    A[0,0] = 0 #needed for recursion
    A[1:,0] = -np.inf
    A[0,1:] = -np.inf

    #changes from local to global found here:
    B[0:,0] = -np.inf #-inf along left
    B[0,1:] = g_open + g_ext * np.arange(0, n, 1) #along top

    C[0,0:] = -np.inf #-inf along top
    C[1:,0] =  g_open + g_ext * np.arange(0, m, 1) #along left

    #needed for computing embeddings
    Model = model
    Model_tokenizer = tokenizer

    #get Ankh embedded sequences
    emb1 = get_embs_Ankh(seq1, model, tokenizer).cpu().numpy()
    emb2 = get_embs_Ankh(seq2, model, tokenizer).cpu().numpy()
    cos = torch.nn.CosineSimilarity(dim=0)

    #2) DP
    for i in range(1, m + 1): #length of first sequence
        for j in range(1, n + 1): #length of second sequence

            sim = cos(torch.tensor(emb1[i - 1] , dtype = torch.float32), torch.tensor(emb2[j - 1] , dtype = torch.float32)).item()
            match_score = sim #no shift down for global

            #filling up A
            A[i,j] = max(A[i-1,j-1] + match_score, B[i-1,j-1] + match_score, C[i-1,j-1] + match_score) #change from local
            #store traceback info
            max_index = np.argmax([A[i-1,j-1] + match_score, B[i-1,j-1] + match_score, C[i-1,j-1] + match_score]) #change from local
            A_tb[i,j] = max_index #0 = from A, 1 = from B, 2 = from C

            #filling up B (gap in seq1)
            B[i,j] = max((A[i,j-1]+g_open+g_ext),(B[i,j-1]+g_ext),(C[i,j-1]+g_open+g_ext))  #change from local
            #store traceback info
            max_index = np.argmax([A[i,j-1] + g_open+g_ext, B[i,j-1] + g_ext, C[i,j-1] + g_open+g_ext])  #change from local
            B_tb[i,j] = max_index #0 = from A, 1 = from B, 2 = from C

            #filling up C (gap in seq2)
            C[i,j] = max((A[i-1,j]+g_open+g_ext),(B[i-1,j]+g_open+g_ext),(C[i-1,j]+g_ext))
            #store traceback info
            max_index = np.argmax([A[i-1,j] + g_open+g_ext, B[i-1,j] + g_open+g_ext, C[i-1,j] + g_ext])
            C_tb[i,j] = max_index #0 = from A, 1 = from B, 2 = from C


    #3) traceback
    #find max cell across all 3 matrices that is the bottom right corner
    A_val = A[m,n]
    B_val = B[m,n]
    C_val = C[m,n]

    max_score = max(A_val, B_val, C_val)

    if A_val == max_score:
        tb_matrix = A_tb
    elif B_val == max_score:
        tb_matrix = B_tb
    elif C_val == max_score:
        tb_matrix = C_tb

    aligned_seq1 = []
    aligned_seq2 = []

    while (i>0) or (j>0):

        #currently in matrix A
        if tb_matrix is A_tb:

            if tb_matrix[i, j] == 0: #from A

                if (i-1)>=0:
                    aligned_seq1.append(seq1[i-1])
                    i -= 1
                else:
                    aligned_seq1.append('-')

                if (j-1)>=0:
                    aligned_seq2.append(seq2[j-1])
                    j -= 1
                else:
                    aligned_seq2.append('-')

                tb_matrix = A_tb

            elif tb_matrix[i, j] == 1: #from B

                if (i-1)>=0:
                    aligned_seq1.append(seq1[i-1])
                    i -= 1
                else:
                    aligned_seq1.append('-')

                if (j-1)>=0:
                    aligned_seq2.append(seq2[j-1])
                    j -= 1
                else:
                    aligned_seq2.append('-')

                tb_matrix = B_tb

            elif tb_matrix[i, j] == 2: #from C

                if (i-1)>=0:
                    aligned_seq1.append(seq1[i-1])
                    i -= 1
                else:
                    aligned_seq1.append('-')

                if (j-1)>=0:
                    aligned_seq2.append(seq2[j-1])
                    j -= 1
                else:
                    aligned_seq2.append('-')

                tb_matrix = C_tb

        #currently in matrix B (gap in seq1)
        elif tb_matrix is B_tb:

            if tb_matrix[i, j] == 0:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
                tb_matrix = A_tb

            elif tb_matrix[i, j] == 1:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
                tb_matrix = B_tb

            elif tb_matrix[i, j] == 2:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
                tb_matrix = C_tb

        #currently in matrix C (gap in seq2)
        elif tb_matrix is C_tb:

            if tb_matrix[i, j] == 0:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
                tb_matrix = A_tb

            elif tb_matrix[i, j] == 1:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
                tb_matrix = B_tb

            elif tb_matrix[i, j] == 2:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
                tb_matrix = C_tb

    #reverse the sequences
    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))

    #remove leading '-' if both sequences start with '-'
    if aligned_seq1.startswith('-') and aligned_seq2.startswith('-'):
        aligned_seq1 = aligned_seq1[1:]
        aligned_seq2 = aligned_seq2[1:]

    return aligned_seq1, aligned_seq2, max_score

#6. affine semi-global alignment model and traceback with e-score scoring function
def affine_semi_global_alignment_escore_and_traceback(seq1, seq2, g_open = -0.25, g_ext = -0.01, model = "", tokenizer = ""):

    #1) initialization
    m = len(seq1); n = len(seq2)
    #3 working matrices:
    #A --> best alignment of seq1[1-i] and seq2[1-j] that aligns seq1[i] and seq2[j]
    #B --> best alignment of seq1[1-i] and seq2[1-j] that aligns gap with seq2[j] (gap in seq1)
    #C --> best alignment of seq1[1-i] and seq2[1-j] that aligns seq1[i] with gap (gap in seq2)
    A = np.zeros([m+1,n+1])
    B = np.zeros([m+1,n+1])
    C = np.zeros([m+1,n+1])

    #3 traceback matrices
    A_tb = np.zeros([m+1,n+1])
    B_tb = np.zeros([m+1,n+1])
    C_tb = np.zeros([m+1,n+1])

    #initializing matrices
    A[0,0] = 0 #needed for recursion
    A[1:,0] = -np.inf
    A[0,1:] = -np.inf

    #changes from semi-global to global found here:
    B[0:,0] = -np.inf #-inf along left (undefined)
    C[0,0:] = -np.inf #-inf along top (undefined)

    #compute embeddings
    Model = model
    Model_tokenizer = tokenizer

    #initalizing the embedding method do we want to use
    #get Ankh embedded sequences
    emb1 = get_embs_Ankh(seq1, model, tokenizer).cpu().numpy()
    emb2 = get_embs_Ankh(seq2, model, tokenizer).cpu().numpy()
    cos = torch.nn.CosineSimilarity(dim=0)

    #2) DP
    for i in range(1, m + 1): #length of first sequence
            for j in range(1, n + 1): #length of second sequence

                sim = cos(torch.tensor(emb1[i - 1] , dtype = torch.float32), torch.tensor(emb2[j - 1] , dtype = torch.float32)).item()
                match_score = sim #no shift down for global

                #filling up A
                A[i,j] = max(A[i-1,j-1] + match_score, B[i-1,j-1] + match_score, C[i-1,j-1] + match_score) #change from local
                #store traceback info
                max_index = np.argmax([A[i-1,j-1] + match_score, B[i-1,j-1] + match_score, C[i-1,j-1] + match_score]) #change from local
                A_tb[i,j] = max_index #0 = from A, 1 = from B, 2 = from C

                #filling up B (gap in seq1)
                B[i,j] = max((A[i,j-1]+g_open+g_ext),(B[i,j-1]+g_ext),(C[i,j-1]+g_open+g_ext))  #change from local
                #store traceback info
                max_index = np.argmax([A[i,j-1] + g_open+g_ext, B[i,j-1] + g_ext, C[i,j-1] + g_open+g_ext])  #change from local
                B_tb[i,j] = max_index #0 = from A, 1 = from B, 2 = from C

                #filling up C (gap in seq2)
                C[i,j] = max((A[i-1,j]+g_open+g_ext),(B[i-1,j]+g_open+g_ext),(C[i-1,j]+g_ext))
                #store traceback info
                max_index = np.argmax([A[i-1,j] + g_open+g_ext, B[i-1,j] + g_open+g_ext, C[i-1,j] + g_ext])
                C_tb[i,j] = max_index #0 = from A, 1 = from B, 2 = from C


    #3) traceback
    #THIS IS DIFFERENT from global:
    #find the maximum value in the last column or last row
    #find the maximum value in the last column or last row for each matrix
    A_last_col = A[:, n]
    B_last_col = B[:, n]
    C_last_col = C[:, n]

    A_last_row = A[m, :]
    B_last_row = B[m, :]
    C_last_row = C[m, :]

    #get the maximum score and its position from any of the last columns or rows
    A_max_col_val = np.max(A_last_col)
    A_max_row_val = np.max(A_last_row)
    B_max_col_val = np.max(B_last_col)
    B_max_row_val = np.max(B_last_row)
    C_max_col_val = np.max(C_last_col)
    C_max_row_val = np.max(C_last_row)

    A_max_val = max(A_max_col_val, A_max_row_val)
    B_max_val = max(B_max_col_val, B_max_row_val)
    C_max_val = max(C_max_col_val, C_max_row_val)

    #determine the overall maximum score
    max_score = max(A_max_val, B_max_val, C_max_val)

    #store the matrix and position of the maximum value
    if A_max_val == max_score:
        if A_max_col_val >= A_max_row_val:
            tb_matrix = A_tb
            position = (np.argmax(A_last_col), n)  #from last column
        else:
            tb_matrix = A_tb
            position = (m, np.argmax(A_last_row))  #from last row

    elif B_max_val == max_score:
        if B_max_col_val >= B_max_row_val:
            tb_matrix = B_tb
            position = (np.argmax(B_last_col), n)  #from last column
        else:
            tb_matrix = B_tb
            position = (m, np.argmax(B_last_row))  #from last row

    elif C_max_val == max_score:
        if C_max_col_val >= C_max_row_val:
            tb_matrix = C_tb
            position = (np.argmax(C_last_col), n)  #from last column
        else:
            tb_matrix = C_tb
            position = (m, np.argmax(C_last_row))  #from last row

    #starting postion
    i, j = position
    #now that we know where to start, max value from the last column or row in any of the 3 matrices
    aligned_seq1 = []
    aligned_seq2 = []

    #we need to add in the free gaps
    #check if gaps need to be added on the right side before starting traceback
    if i < m:
        for k in range(m, i, -1):  #add gaps until the end of seq1 is aligned
            aligned_seq1.append(seq1[k-1])
            aligned_seq2.append('-')

    if j < n:
        for l in range(n, j, -1):  #add gaps until the end of seq2 is aligned
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[l-1])


    while (i>0) or (j>0):

        #currently in matrix A
        if tb_matrix is A_tb:

            if tb_matrix[i, j] == 0: #from A

                if (i-1)>=0:
                    aligned_seq1.append(seq1[i-1])
                    i -= 1
                else:
                    aligned_seq1.append('-')

                if (j-1)>=0:
                    aligned_seq2.append(seq2[j-1])
                    j -= 1
                else:
                    aligned_seq2.append('-')

                tb_matrix = A_tb

            elif tb_matrix[i, j] == 1: #from B

                if (i-1)>=0:
                    aligned_seq1.append(seq1[i-1])
                    i -= 1
                else:
                    aligned_seq1.append('-')

                if (j-1)>=0:
                    aligned_seq2.append(seq2[j-1])
                    j -= 1
                else:
                    aligned_seq2.append('-')

                tb_matrix = B_tb

            elif tb_matrix[i, j] == 2: #from C

                if (i-1)>=0:
                    aligned_seq1.append(seq1[i-1])
                    i -= 1
                else:
                    aligned_seq1.append('-')

                if (j-1)>=0:
                    aligned_seq2.append(seq2[j-1])
                    j -= 1
                else:
                    aligned_seq2.append('-')

                tb_matrix = C_tb

        #currently in matrix B (gap in seq1)
        elif tb_matrix is B_tb:

            if tb_matrix[i, j] == 0:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
                tb_matrix = A_tb

            elif tb_matrix[i, j] == 1:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
                tb_matrix = B_tb

            elif tb_matrix[i, j] == 2:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
                tb_matrix = C_tb

        #currently in matrix C (gap in seq2)
        elif tb_matrix is C_tb:

            if tb_matrix[i, j] == 0:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
                tb_matrix = A_tb

            elif tb_matrix[i, j] == 1:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
                tb_matrix = B_tb

            elif tb_matrix[i, j] == 2:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
                tb_matrix = C_tb

    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))

    #remove leading '-' if both sequences start with '-'
    if aligned_seq1.startswith('-') and aligned_seq2.startswith('-'):
        aligned_seq1 = aligned_seq1[1:]
        aligned_seq2 = aligned_seq2[1:]

    return aligned_seq1, aligned_seq2, max_score

if __name__ == "__main__":
    main()