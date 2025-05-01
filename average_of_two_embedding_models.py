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


#2. main function to run a single alignment example, using the average embeddings made by models ProtT5 and Ankh
#example
def main():
    
    protein_seq1 = "AACD" #any protein sequences of capital letters
    protein_seq2 = "CDDA" #any protein sequence of capital letters
    
    gap_open = -1 #any value can be used
    gap_extension = -0.1  #any value can be used

    #code for initalizing 2 different embedding models
    model_1, tokenizer_1 = Ankh_initialize()
    model_2, tokenizer_2 = ProtT5_initialize()


    alignment, max_score = affine_local_alignment_model_and_traceback(protein_seq1, protein_seq2, gap_open, gap_extension, model_1, tokenizer_1, model_2, tokenizer_2)

    aligned_seq1 = alignment[0]
    aligned_seq2 = alignment[1]

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


def ProtT5_initialize():

    ProtT5 = T5EncoderModel.from_pretrained("/change/path/to/model", local_files_only=True)
    ProtT5_tokenizer = T5Tokenizer.from_pretrained("/change/path/to/model", do_lower_case=False)
    ProtT5 = ProtT5.half()
    ProtT5 = ProtT5.eval()
    ProtT5 = ProtT5.to(device)

    return ProtT5, ProtT5_tokenizer

def get_embs_T5(seq, T5, ProtT5_tokenizer):

  seq = re.sub(r"[UZOB]", "X", seq)
  seq = " ".join(list(seq))
  inputs = ProtT5_tokenizer([seq], padding=True, return_tensors="pt")
  inputs = {key: value.to(device) for key, value in inputs.items()}

  with torch.no_grad(): 
    outputs = T5(**inputs) 
    last_hidden_states = outputs.last_hidden_state 
  
  sequence_embedding = last_hidden_states[0, :-1, :]
  return sequence_embedding



#4. compute the local alignment using the average of two embedding models

def l2_normalize(vectors, axis=-1, epsilon=1e-12):
    norm = np.linalg.norm(vectors, axis=axis, keepdims=True)
    return vectors / (norm + epsilon)


#MODEL DP AND TRACEBACK
def affine_local_alignment_model_and_traceback(seq1, seq2, g_open = -1, g_ext = -0.25, model = "", tokenizer = "", model2 = "", tokenizer2 = ""):

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

  cos = torch.nn.CosineSimilarity(dim=0)

  #compute embeddings Ankh - size 1536 embeddings
  model_1 = model
  model_tokenizer_1 = tokenizer
  emb1_ankh = get_embs_Ankh(seq1, model_1, model_tokenizer_1).cpu().numpy() 
  emb2_ankh = get_embs_Ankh(seq2, model_1, model_tokenizer_1).cpu().numpy() 

  #compute embeddings ProtT5 - size 1024 embeddings
  model_2 = model2
  model_tokenizer_2 = tokenizer2
  emb1_protT5 = get_embs_T5(seq1, model_2, model_tokenizer_2).cpu().numpy()
  emb2_protT5 = get_embs_T5(seq2, model_2, model_tokenizer_2).cpu().numpy()

  #make same dim
  emb1_ankh = emb1_ankh[:, :1024] 
  emb2_ankh = emb2_ankh[:, :1024] 
  #now all embeddings have shape (N, 1024)
  
  #stack them 
  vectors_seq1 = np.array([emb1_ankh, emb1_protT5])  #NumPy array
  vectors_seq2 = np.array([emb2_ankh, emb2_protT5])  #NumPy array
  #the results will be arrays of shape (2, N, 1024)

  #normalize
  vectors_seq1 = l2_normalize(vectors_seq1, axis=2) #axis = 2 actually means axis 1 here, and axis = 1 means axis 0 here
  vectors_seq2 = l2_normalize(vectors_seq2, axis=2)  

  #compute average embedding
  #averaging the embeddings from the two models for each position in the sequence
  average_embedding_seq1 = np.mean(vectors_seq1, axis=0) #across the column
  average_embedding_seq2 = np.mean(vectors_seq2, axis=0)


  #2 ) DP
  for i in range(1, m + 1): #length of first sequence
      for j in range(1, n + 1): #length of second sequence

        sim = cos(torch.tensor(average_embedding_seq1[i - 1] , dtype = torch.float32), torch.tensor(average_embedding_seq2[j - 1] , dtype = torch.float32)).item()
        match_score = sim

        #filling up A
        A[i,j] = max(A[i-1,j-1] + match_score, B[i-1,j-1] + match_score, C[i-1,j-1] + match_score, 0)
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

    #currently in matrix C
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

  aligned_seq1 = ''.join(reversed(aligned_seq1))
  aligned_seq2 = ''.join(reversed(aligned_seq2))

  return [aligned_seq1, aligned_seq2], max_score


if __name__ == "__main__":
    main()