#!/change/this/to/point/to/your/virtual_env/bin/python
import os
import warnings
import re
from Bio import SeqIO
from scipy import stats
import blosum as bl
import numpy as np
import pandas as pd
import torch
import sentencepiece
import shutil
from tqdm import tqdm
from transformers import T5Model
from transformers import T5Tokenizer, T5EncoderModel, AutoTokenizer, AutoModelForSeq2SeqLM
from transformers import T5ForConditionalGeneration
import sys

torch.cuda.empty_cache()
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


#EXTRACTING ALIGNMENTS BASED OFF REFERENCE
def intersect_ranges(range1, range2):
    start1, end1 = range1
    start2, end2 = range2

    #find the intersection
    start_intersect = max(start1, start2)
    end_intersect = min(end1, end2)

    #check if they actually overlap
    if start_intersect < end_intersect:
        return (start_intersect, end_intersect)
    else:
        return None, None


def find_nth_non_dash(s, n):
    count = 0  # To count non-dash characters
    for i, char in enumerate(s):
        if char != '-':
            count += 1
            if count == n:
                return i
    return -1

def count_non_dashes_before_index(s, index):
    count = 0
    for i in range(index):
        if s[i] != '-':
            count += 1
    return count

def extract_alignment(protein1, protein2, ref1, ref2, algn1, algn2):

    ref1_no_gaps = ref1.replace('-','')
    ref2_no_gaps = ref2.replace('-','')
    algn1_no_gaps = algn1.replace('-','')
    algn2_no_gaps = algn2.replace('-','')

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

    #we now need to figure out the range of overlaps for the first algnment sequence
    range1 = (ref1_in_protein1_start,ref1_in_protein1_end) #where reference is
    range2 = (algn1_in_protein1_start, algn1_in_protein1_end) #where alignment is
    algnment1_intersecting_reference_start, algnment1_intersecting_reference_end = intersect_ranges(range1, range2)

    #worst case:
    if algnment1_intersecting_reference_start == None:
        algn1_final = str(ref1_no_gaps) + str("-"*len(ref2_no_gaps))
        algn2_final = str("-"*len(ref1_no_gaps)) + str(ref2_no_gaps)
        return algn1_final, algn2_final

    len_overlap_1 = algnment1_intersecting_reference_end - algnment1_intersecting_reference_start

    #we need to figure out the indices in the first algnment the overlap corresponds to
    algnment1_index_start = algnment1_intersecting_reference_start - algn1_in_protein1_start
    algnment1_index_start = find_nth_non_dash(algn1, algnment1_index_start) + 1

    num = count_non_dashes_before_index(algn1, algnment1_index_start)

    #the ending index of the algnment covered, we need to find the len_overlap_1th AA in the algnment
    algnment1_index_end = find_nth_non_dash(algn1, num+len_overlap_1+1)

    if algnment1_index_end == -1:
        algnment1_index_end = len(algn1)

    #this is the "usuable range" of the first algnment
    usable_range_algn1 = (algnment1_index_start,algnment1_index_end)

    #same thing for the bottom sequence
    #we now need to figure out the range of overlaps for the first algnment sequence
    range1 = (ref2_in_protein2_start,ref2_in_protein2_end) #where reference is
    range2 = (algn2_in_protein2_start, algn2_in_protein2_end)

    algnment2_intersecting_reference_start, algnment2_intersecting_reference_end = intersect_ranges(range1, range2)

    if algnment2_intersecting_reference_start == None:
        algn1_final = str(ref1_no_gaps) + str("-"*len(ref2_no_gaps))
        algn2_final = str("-"*len(ref1_no_gaps)) + str(ref2_no_gaps)
        return algn1_final, algn2_final

    len_overlap_2 = algnment2_intersecting_reference_end - algnment2_intersecting_reference_start
    #we need to figure out the indices in the first algnment the overlap corresponds to
    algnment2_index_start = algnment2_intersecting_reference_start - algn2_in_protein2_start
    #the start index is algnment2_index_start none gap characters
    algnment2_index_start = find_nth_non_dash(algn2, algnment2_index_start) + 1

    num = count_non_dashes_before_index(algn2, algnment2_index_start)
    algnment2_index_end = find_nth_non_dash(algn2, num+len_overlap_2+1)

    if algnment2_index_end == -1:
        algnment2_index_end = len(algn2)

    usable_range_algn2 = (algnment2_index_start,algnment2_index_end)

    #intersect the ranges in the algnment
    usable_range_start, usable_range_end = intersect_ranges(usable_range_algn1, usable_range_algn2)

    if usable_range_start == None:
        algn1_final = str(ref1_no_gaps) + str("-"*len(ref2_no_gaps))
        algn2_final = str("-"*len(ref1_no_gaps)) + str(ref2_no_gaps)
        return algn1_final, algn2_final

    overlap_algn1 = algn1[usable_range_start:usable_range_end]
    overlap_algn2 = algn2[usable_range_start:usable_range_end]

    #we now have the overlapping part of the algnment that covers the reference
    #we now have to check if there are any parts of the algnment not covering the reference
    #easiest way is to remove gaps from the overlaping parts and find them in the references without gaps and get the top left and right, and bottom left and right

    overlap_algn1_no_gaps = overlap_algn1.replace("-","")
    overlap_algn2_no_gaps = overlap_algn2.replace("-","")

    #we fully covered the algnment
    if (len(overlap_algn1_no_gaps) == len(ref1_no_gaps)) and (len(overlap_algn2_no_gaps) == len(ref2_no_gaps)):
        return overlap_algn1, overlap_algn2

    top_left = ""
    top_right = ""
    bottom_left = ""
    bottom_right = ""

    #getting the rest of the first alignment sequence
    algn1_start_in_ref = ref1_no_gaps.find(overlap_algn1_no_gaps)
    top_left = ref1_no_gaps[:algn1_start_in_ref]
    top_right = ref1_no_gaps[algn1_start_in_ref + len(overlap_algn1_no_gaps):]

    #getting the rest of the second alignment sequence
    algn2_start_in_ref = ref2_no_gaps.find(overlap_algn2_no_gaps)
    bottom_left = ref2_no_gaps[:algn2_start_in_ref]
    bottom_right = ref2_no_gaps[algn2_start_in_ref + len(overlap_algn2_no_gaps):]

    #reconstructing the alignment
    algn1_final = str(top_left) + str("-"*len(bottom_left)) + str(overlap_algn1) + str(top_right) + str("-"*len(bottom_right))
    algn2_final = str("-"*len(top_left)) + str(bottom_left) + str(overlap_algn2) + str("-"*len(top_right)) + str(bottom_right)

    return algn1_final, algn2_final


#ankh 
def Ankh_initialize():
    ankh_model = T5ForConditionalGeneration.from_pretrained("/change/to/point/to/ANKH", local_files_only=True)
    ankh_tokenizer = AutoTokenizer.from_pretrained("/change/to/point/to/ANKH", local_files_only=True)

    # ankh_model = ankh_model.half()
    ankh_model = ankh_model.eval()
    ankh_model = ankh_model.to(device)

    return ankh_model, ankh_tokenizer


def get_embs_Ankh(seq, ankh_model, tokenizer):
    
    # seq = re.sub(r"[UZOB]", "X", seq)
    # seq = " ".join(list(seq))

    inputs = tokenizer([seq], padding=True, return_tensors="pt") #pytorch format 
    inputs = {key: value.to(device) for key, value in inputs.items()}

    #get the embeddings (last hidden state)
    with torch.no_grad(): #no gradients since we are not training
        # outputs = ankh_model(**inputs) #send inputs through the model
        outputs = ankh_model.encoder(**inputs)

        last_hidden_states = outputs.last_hidden_state #get an embedding for each token
    
    sequence_embedding = last_hidden_states[0, :-1, :]
    
    return sequence_embedding



#MODEL DP AND TRACEBACK
def affine_local_alignment_model_and_traceback(seq1, seq2, g_open = -1, g_ext = -0.25, model = "", tokenizer = ""):

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

  #compute embeddings
  Model = model
  Model_tokenizer = tokenizer

  emb1 = get_embs_Ankh(seq1, model, tokenizer).cpu().numpy()
  emb2 = get_embs_Ankh(seq2, model, tokenizer).cpu().numpy()
  cos = torch.nn.CosineSimilarity(dim=0)

  #2 ) DP
  for i in range(1, m + 1): #length of first sequence
      for j in range(1, n + 1): #length of second sequence

        sim = cos(torch.tensor(emb1[i - 1] , dtype = torch.float32), torch.tensor(emb2[j - 1] , dtype = torch.float32)).item()
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

#DISTANCES:
def lists_X(algn):
  #get the first and second seqeuences in the alignment
  prot1 = algn[0]
  prot2 = algn[1]

  #counter use to print out the indices in the encodings
  c1 = 0
  c2 = 0
  #encoded alignment needed for ssp distance
  d_ssp1 = []
  d_ssp2 = []
  #encoded alignment needed for seq distance
  d_seq1 = []
  d_seq2 = []
  #encooded alignment needed for pos distance
  d_pos1 = []
  d_pos2 = []

  #go through each character in the proteins
  for i in range(len(prot1)): #prot1 and prort2 have the same length since they are aligned to each other

    #if the current character is a gap (in the first protein)
    if prot1[i] == "-":
      d_ssp1.append("-") #ssp distance "ignores" gaps so in the encoded version we just leave it
      d_seq1.append("G-" + str(1)) #seq distance treats all gaps in a protein the same, and since this is protein 1, encode the gap as "G-1"
      d_pos1.append("G-" + str(1) + "-" + str(c1)) #pos distance keeps positsional info about the gaps, its has a number which is the index of the sequence to the left of the gap

    #the current character in the first protein is not a gap
    elif prot1[i] != "-":
      c1 += 1 #increment the index of which character in protein 1 we are at
      #all encodings treate the non-gap characters the same, just print the index
      d_ssp1.append("S-" + str(1) + "-" + str(c1))
      d_seq1.append("S-" + str(1) + "-" + str(c1))
      d_pos1.append("S-" + str(1) + "-"+ str(c1))

    #repeat above steps for the second protein sequence
    if prot2[i] == "-":
      d_ssp2.append("-")
      d_seq2.append("G-" + str(2))
      d_pos2.append("G-" + str(2) + "-" + str(c2))

    elif prot2[i] != "-":
      c2 += 1
      d_ssp2.append("S-" + str(2) + "-" + str(c2))
      d_seq2.append("S-" + str(2) + "-" + str(c2))
      d_pos2.append("S-" + str(2) + "-" + str(c2))

  #return the 3 encoded versions of the alignment
  return [d_ssp1 , d_ssp2] , [d_seq1 , d_seq2] , [d_pos1 , d_pos2]


def get_h(x1 , x2 , type = "ssp"):
  h1 = [] #will store the other characters in the second sequence relative to the first sequence
  h2 = [] #will store the other characters in the first sequence relative to the second sequence

  for i in range(len(x1)):
    #creating the homolgy set for the first sequence (will contain the characters in the other sequence)
    #checking the character is not a gap, since gaps do not have homology sets of their own
    if x1[i][0] == "S":
      if x2[i][0] != "-": #if the corresponding character in the second sequence is not a gap we add the character to the homology set of the current character in the first sequence
        h1.append(x2[i])
      elif type == "ssp":
        h1.append("-")
    #creating the homology set for the second sequence (fill contain the characters in the other sequence)
    if x2[i][0] == "S":
      if x1[i][0] != "-":
        h2.append(x1[i])
      elif type == "ssp":
        h2.append("-")

  return h1, h2

#this function computes the SSP distance (Jacard) given the homolgy sets for the sequences to be compared
def d_ssp(H_A , H_B):
  #UPDATED
  up = 0 #accumulate the number of common elements between both sequences, where neither is - (intersection of homology sets)
  down = 0 #number of unique elements between both sequences, where neither is -
  for i in range(len(H_A)):
    for j in range(len(H_A[i])):
      if H_A[i][j] == "-" or H_B[i][j] == "-":
        down += len(set([H_A[i][j]]).union([H_B[i][j]]))
      else:
        up += len(set([H_A[i][j]]).intersection([H_B[i][j]]))
        down += len(set([H_A[i][j]]).union([H_B[i][j]]))

  #if all postions have gaps, there is nothing to compare, this is the worst case alignment
  if down == 0:
    return 1
    
  return 1 - (up / down)

#this function computes the symmetric difference between the homolgy sets for a given alignment
#this is used to compute both the seq and pos distances
def d_X(H_A , H_B , n , m):
  d = 0
  for i in range(len(H_A)):
    for j in range(len(H_A[i])):
      up = len(set([H_A[i][j]]).symmetric_difference([H_B[i][j]]))
      down = len(set([H_A[i][j]])) + len(set([H_B[i][j]]))
      d += up / down

  return (1 / (n + m)) * d


##### GETTING THE 3 ORIGINAL DISTANCES #####
def get_3_dists(algn1, algn2):

  #first we need to encode the 4 sequences (2 alignments)
  #the lists_X function encodes the alignments in 3 ways
  #way 1: needed for ssp, gaps are ignored
  #way 2: needed for seq, all gaps in each sequence are treated the same
  #way 3: needed for pos, keep positional info about the gaps
  x_A = lists_X(algn1)
  x_B = lists_X(algn2)

  len_x1 = len(algn1[0].replace("-" , ""))
  len_x2 = len(algn1[1].replace("-" , ""))

  #accessing the ssp encoded first alignment
  x1_A = x_A[0][0]
  x2_A = x_A[0][1]  # ssp = 0 , seq = 1 , pos = 2
  #accessing the ssp encoded second alignment
  x1_B = x_B[0][0]
  x2_B = x_B[0][1]
  #getting the homology sets
  #the default type is ssp, so we dont have to specify it as a parameter here
  H_A = get_h(x1_A , x2_A) #the 2 sequences of the first alignment
  H_B = get_h(x1_B , x2_B) #the 2 sequences of the second alignment

  #H_A and H_B both store 2 lists, one which represents the homology sets for the first sequence
  #and the other one is for the homology sets of the second sequence

  #computing the ssp distance (Jaccard)
  ssp_dis = d_ssp(H_A , H_B)

  #repeat above sets to get the seq distance
  #get the seq specific encoded version of the alignment
  #sequences for the first alignment
  x1_A = x_A[1][0]
  x2_A = x_A[1][1]  # ssp = 0 , seq = 1 , pos = 2
  #sequences for the second alignment
  x1_B = x_B[1][0]
  x2_B = x_B[1][1]

  #get the homology sets for this encoding
  H_A = get_h(x1_A , x2_A , type = "seq")
  H_B = get_h(x1_B , x2_B , type = "seq")

  #compute the seq distance, this calls d_X which just computes the symmetric difference between the homology sets
  seq_dis = d_X(H_A , H_B , len_x1 , len_x2)

  #repeat the above steps to get the pos distance
  #getting the pos encoded version of the alignment
  x1_A = x_A[2][0]
  x2_A = x_A[2][1]  # ssp = 0 , seq = 1 , pos = 2

  x1_B = x_B[2][0]
  x2_B = x_B[2][1]

  #getting the homology sets
  H_A = get_h(x1_A , x2_A , type = "pos")
  H_B = get_h(x1_B , x2_B , type = "pos")

  #getting the distance, this calls d_X which just computes the symmetric difference between the homology sets
  pos_dis = d_X(H_A , H_B , len_x1 , len_x2)

  return float(ssp_dis) , float(seq_dis) , float(pos_dis)

#Distance 4: d_d relative displacement distance
#this gets the indices of the characters that are not gaps (helper function for d_d)
def get_no_gap_indexes(string):
  return [(i + 1) for i, s in enumerate(string) if '-' not in s]

#given two alignments, this function iterates over the indexes of non-gap positions in the first string (p1a1_no_gap_idxs) and calculates a tmp value
#which accumulates a measure of the difference between these positions in both alignments (algn1 and algn2).
def d_d(algn1, algn2):
  p1a1 = algn1[0]
  p2a1 = algn1[1]

  p1a2 = algn2[0]
  p2a2 = algn2[1]

  p1a1_no_gap_idxs = get_no_gap_indexes(p1a1)
  p2a1_no_gap_idxs = get_no_gap_indexes(p2a1)

  p1a2_no_gap_idxs = get_no_gap_indexes(p1a2)
  p2a2_no_gap_idxs = get_no_gap_indexes(p2a2)

  distance = 0
  tmp = 0
  for i in range(len(p1a1_no_gap_idxs)):
    tmp = 0
    for j in range(len(p2a1_no_gap_idxs)):
      tmp += abs((p1a1_no_gap_idxs[i] - p2a1_no_gap_idxs[j]) - (p1a2_no_gap_idxs[i] - p2a2_no_gap_idxs[j]))
    distance += tmp

  n = len(p1a1.replace("-" , ""))
  m = len(p2a1.replace("-" , ""))

  return float((1 / ((m * n) * (m + n))) * distance)


#Distance 5: d_cc closest context distance
#this helper function creates a list where non-dash characters are located within each sequence (helper function for get_d_p)
def aligned_to_indexed(seqs):
  no_dash = []
  positions = []
  for seq in seqs:
    no_dash.append(seq.replace("-" , ""))
    pos = []
    for i , char in enumerate(seq):
      if char != "-":
        pos.append(i)
    positions.append(pos)

  return no_dash, positions

def get_d_P(first_align, second_align , trace = False):
  first_prot = first_align[0].replace("-" , "")
  second_prot = first_align[1].replace("-" , "")

  P_1 = first_align[0]
  Q_1 = first_align[1]

  P_2 = second_align[0]
  Q_2 = second_align[1]

  P1_indexed = aligned_to_indexed([P_1])
  Q1_indexed = aligned_to_indexed([Q_1])

  P2_indexed = aligned_to_indexed([P_2])
  Q2_indexed = aligned_to_indexed([Q_2])

  final_d = 0
  for i in range(len(first_prot)):
    if trace:
      print("------------------------")
      print("Step : " + str(i + 1))

    current_char = first_prot[i]

    if trace:
      print("Current Char : " + current_char)

    char_d = 0

    P1_corresp_idx = P1_indexed[1][0][i]
    P2_corresp_idx = P2_indexed[1][0][i]

    opp_char_Q1 = Q_1[P1_corresp_idx]

    if trace:
      print("Opp Char : " + opp_char_Q1)

    start_letter_idx = P2_corresp_idx
    # print(start_letter_idx)
    target_letter = opp_char_Q1

    if target_letter == "-":
      char_d = abs(len(Q_1[:P1_corresp_idx].replace("-" , "")) - len(Q_2[:P2_corresp_idx].replace("-" , "")))
      final_d += char_d
      if trace:
        print("char d : "  + str(char_d))
      continue

    if target_letter == Q_2[start_letter_idx]:
      # d = 0
      if trace:
        print("char d : "  + str(0))
      continue

    found_char_right_idx = Q_2[start_letter_idx + 1 :].find(target_letter)
    # print(Q_2[start_letter_idx + 1 :])
    # print(found_char_right_idx)

    found_char_left_idx = Q_2[: start_letter_idx].find(target_letter)
    # print(Q_2[: start_letter_idx])
    # print(found_char_left_idx)


    if found_char_right_idx == -1:
      if trace:
        print("char left idx : " + str(found_char_left_idx))
      letters_count_left = len(Q_2[found_char_left_idx : start_letter_idx].replace("-" , ""))
      if trace:
        print("Lefties : " + Q_2[found_char_left_idx : start_letter_idx])
      char_d = letters_count_left
      final_d += char_d

    elif found_char_left_idx == -1:
      found_char_right_idx += (start_letter_idx + 1)
      if trace:
        print("char right idx : " + str(found_char_right_idx))
      letters_count_right = len(Q_2[start_letter_idx + 1 : found_char_right_idx + 1].replace("-" , ""))
      if trace:
        print("Righties : " + Q_2[start_letter_idx + 1 : found_char_right_idx + 1] )
      char_d = letters_count_right
      final_d += char_d

    else:
      found_char_right_idx += (start_letter_idx + 1)

      letters_count_right = len(Q_2[start_letter_idx + 1: found_char_right_idx  + 1].replace("-" , ""))
      if trace:
        print("char right idx : " + str(found_char_right_idx))
        print("Righties : " + Q_2[start_letter_idx + 1 : found_char_right_idx + 1] )
      letters_count_left = len(Q_2[found_char_left_idx : start_letter_idx].replace("-" , ""))
      if trace:
        print("char left idx : " + str(found_char_left_idx))
        print("Lefties : " + Q_2[found_char_left_idx : start_letter_idx])

      char_d = min(letters_count_right , letters_count_left)

      final_d += char_d

    if trace:
      print("char d : "  + str(char_d))

  return final_d


def d_cc(first_align , second_align):

  n = len(first_align[0].replace("-" , ""))
  m = len(first_align[1].replace("-" , ""))

  ### d_p_1
  x1 = first_align[0]
  y1 = first_align[1]

  x2 =  second_align[0]
  y2 = second_align[1]

  d_p_1 = get_d_P([x1, y1] , [x2 , y2])

  ### d_p_2
  x1 = second_align[0]
  y1 = second_align[1]

  x2 = first_align[0]
  y2 = first_align[1]

  d_p_2 = get_d_P([x1, y1] , [x2 , y2])

  ### d_q_1
  x1 = first_align[1]
  y1 = first_align[0]

  x2 = second_align[1]
  y2 = second_align[0]

  d_q_1 = get_d_P([x1, y1] , [x2 , y2])

  ### d_q_2
  x1 = second_align[1]
  y1 = second_align[0]

  x2 = first_align[1]
  y2 = first_align[0]

  d_q_2 = get_d_P([x1, y1] , [x2 , y2])

  ###total
  d_P = d_p_1 + d_p_2
  d_Q = d_q_1 + d_q_2

  total_d = d_P + d_Q

  #normalize
  return ((1 / (4 * n * m)) * total_d)

#OUTPUT FILE FUNCTIONS
def aligned_to_indexed(seqs):
  no_dash = []
  positions = []
  for seq in seqs:
    no_dash.append(seq.replace("-" , ""))
    pos = []
    for i , char in enumerate(seq):
      if char != "-":
        pos.append(i)
    positions.append(pos)

  return no_dash, positions



def get_visualization(prot1, prot2 , score , Type = "" , Model = "" , Model_Tokenizer = ""):

  MODELS_LIST = ["ProtT5" , "ProtBert" , "ProtAlbert" , "XLNet" , "ESM1b" , "ESM2", "ankh"]
  BLOSUMS_LIST = ["BLOSUM45" , "BLOSUM50" , "BLOSUM62" , "BLOSUM80" , "BLOSUM90"]
  cos = torch.nn.CosineSimilarity(dim=0)
  sim = 0

  seqs = [prot1 , prot2]
  no_dash , positions = aligned_to_indexed(seqs)

  model = Model
  tokenizer = Model_Tokenizer

  p1_revived = ""
  p2_revived = ""
  aligned_info = ""

  for i in range(len(prot1)):

    if i in positions[0]:
      p1_revived += prot1[i]
    else:
      p1_revived += "-"

    if i in positions[1]:
      p2_revived += prot2[i]
    else:
      p2_revived += "-"

    if p1_revived[-1] == p2_revived[-1]:
      aligned_info += p1_revived[-1]

    elif p1_revived[-1] == "-" or p2_revived[-1] == "-":
      aligned_info += " "

    elif p1_revived[-1] != p2_revived[-1]:
      if True:
        sim = 0.5

      if sim > 0:
         aligned_info += " "

  return p1_revived , aligned_info, p2_revived, score

def length_matcher(x , y , place = ""):
  length = 5

  if len(x) < length:
    spaces = abs(len(x) - length)

    if place == "Back":
      x = " " * spaces + x
    if place == "Front":
      x = x + " " * spaces

  if len(y) < length:
    spaces = abs(len(y) - length)

    if place == "Back":
      y = " " * spaces + y
    if place == "Front":
      y = y + " " * spaces

  return x, y

#Returns a single string that displays the alignment split nicely into segments for writing to a file.
def get_formatted_alignent(prot1, prot2, Type = "", score = 0 , Model = "" , Model_Tokenizer = "", segment_length = 60):
  p1_pos = 1
  p2_pos = 1
  aligned_gaps = ""
  output = ""
  p1_al , aligned_info , p2_al , al_score = get_visualization(prot1 , prot2, score, Type, Model, Model_Tokenizer)

  for j in range(int(len(prot1) / segment_length) + 1):
    p1_posix = p1_al[j * segment_length: (j + 1) * segment_length]
    p2_posix = p2_al[j * segment_length: (j + 1) * segment_length]

    if set(list(p1_posix)) == set("-"):
      back1 = p1_pos - 1
    else:
      back1 = p1_pos
    if set(list(p2_posix)) == set("-"):
      back2 = p2_pos - 1
    else:
      back2 = p2_pos

    p1_back_str, p2_back_str = length_matcher(str(back1) , str(back2) , place = "Front")

    for k in range(len(p1_posix)):
      if p1_posix[k] != "-":
        p1_pos += 1
      if p2_posix[k] != "-":
        p2_pos += 1

    p1_end_str, p2_end_str = length_matcher(str(p1_pos - 1) , str(p2_pos - 1) , place = "Back")
    aligned_gaps = " " * len(p1_back_str)

    output += ("Seq 1 : " + p1_back_str + " " + p1_al[j * segment_length: (j + 1) * segment_length] + " " + p1_end_str)
    output += ("\n")
    output += ("        "  +  aligned_gaps + " " + aligned_info[j * segment_length: (j + 1) * segment_length])
    output += ("\n")
    output += ("Seq 2 : "  + p2_back_str + " " + p2_al[j * segment_length: (j + 1) * segment_length] + " " + p2_end_str)
    output += ("\n")
    output += ("\n")

  return output




#getting the output file
def run_individual_method_tests(ID, example_pairs, gap_open, gap_extension, model="", tokenizer="", N=-1):

  examples = pd.read_csv(example_pairs) #each row is an example
  num_examples = examples.shape[0] #how many rows we have
  if N == -1:
      num_examples = examples.shape[0] #no limit specified
  else:
      if N > num_examples:
          print("there are only " + str(num_examples) + " examples.")
      else: num_examples = N #we will only use the first N rows


  #create a directory to store files for each of the 'num_examples' examples
  os.makedirs(os.path.join(os.getcwd(), f"{ID}_tests"), exist_ok=True)

  all_distances = pd.DataFrame(columns=["SSP", "SEQ", "POS", "DD", "CC"])

  for n in range(0,num_examples):

    #get the current row of the examples table
    current_example = examples.iloc[n] #columns: protein1_ID, protein2_ID, msa1_ID, msa2_ID, protein1, protein2, ref1, ref2, p1_mStartIndex, p2_mStartIndex, p1_mEndIndex, p2_mEndIndex
    #create a file to store the results of testing this example with an informative name
    protein1_ID = current_example["protein1_ID"]
    protein2_ID = current_example["protein2_ID"]
    msa1_ID = current_example["msa1_ID"]
    msa2_ID = current_example["msa2_ID"]
    protein1_seq =  current_example["protein1"]
    protein2_seq =  current_example["protein2"]
    ref1 = current_example["ref1"]
    ref2 = current_example["ref2"]
    p1_mStartIndex = current_example["p1_mStartIndex"]
    p2_mStartIndex = current_example["p2_mStartIndex"]
    p1_mEndIndex = current_example["p1_mEndIndex"]
    p2_mEndIndex = current_example["p2_mEndIndex"]

    directory_path = os.path.join(os.getcwd(), f"{ID}_tests/{ID}_test_{n + 1}")
    os.makedirs(directory_path, exist_ok=True)
    method = "Ankh"
    
    file_name = f"{ID}_test{n+1}_{method}_{gap_open:1.5f}_{gap_extension:1.5f}.txt"
    with open(file_name, 'w') as file: #create the file, if one with the same file name already exists, we override
            #write basic info about the current test to the top of file
            file.write("Test " + str(n + 1) + " of " + str(num_examples))

            #print the protein sequences
            file.write("\nProtein sequences in FASTA:\n")
            file.write("\nProtein 1:")
            file.write("\n>" + str(protein1_ID))
            file.write("\n" +str(protein1_seq))
            file.write("\nLength of protein sequence 1: " + str(len(str(protein1_seq))))
            file.write("\nContains MSA sequence1: " + str(msa1_ID) + "\n" + str(ref1) + "\n")
            file.write("At index: " + str(p1_mStartIndex) + ".." + str(p1_mEndIndex) + "\n")

            file.write("\n\nProtein 2:")
            file.write("\n>" + str(protein2_ID))
            file.write("\n" +str(protein2_seq))
            file.write("\nLength of protein sequence 2: " + str(len(str(protein2_seq))))
            file.write("\nContains MSA sequence2: " + str(msa2_ID) + "\n" + str(ref2) + "\n")
            file.write("At index: " + str(p2_mStartIndex) + ".." + str(p2_mEndIndex) + "\n")

            #msa info
            file.write("\n\nFormatted alignment as induced by reference MSA - golden standard alignment:")
            file.write("\n" + get_formatted_alignent(ref1, ref2))


            
            
            alignment, score = affine_local_alignment_model_and_traceback(protein1_seq, protein2_seq, gap_open, gap_extension, model = model, tokenizer = tokenizer)

            file.write("\nScoring Method: " + str(method))
            file.write("\nParameters:")
            file.write("\n  Gap Open Penalty: " + str(gap_open))
            file.write("\n  Gap Extension Penalty: " + str(gap_extension))
            file.write("\nScore: " + str(score))

            #writing the actual alignment
            file.write("\n\nRaw Alignment:")
            file.write("\n"+get_formatted_alignent(alignment[0], alignment[1], method, score, model, tokenizer))

            #extract the alignment
            file.write("\n\nExtracted Alignment:")

            extract_alignment1, extract_alignment2 = extract_alignment(protein1_seq, protein2_seq, ref1, ref2, alignment[0],alignment[1])
            file.write("\n" + get_formatted_alignent(str(extract_alignment1), str(extract_alignment2), method, score, model, tokenizer))

            #### at this point, we have the first methods alignment (raw, and extracted), Now we compute the 5 distances between the extracted alignment and reference #####
            alignment_generated = (extract_alignment1, extract_alignment2)
            #get all 5 distances:

            ssp , seq, pos = get_3_dists((ref1,ref2), alignment_generated)
            dd = d_d((ref1,ref2), alignment_generated)
            cc = d_cc((ref1,ref2), alignment_generated)

            file.write("\n\nDistances:")
            file.write(f"\n  SSP: {ssp:.6f}")
            file.write(f"\n  SEQ: {seq:.6f}")
            file.write(f"\n  POS: {pos:.6f}")
            file.write(f"\n  DD: {dd:.6f}")
            file.write(f"\n  CC: {cc:.6f}")

            new_row = {
                "SSP" : ssp,
                "SEQ" : seq,
                "POS" : pos,
                "DD"  : dd,
                "CC"  : cc
            }

            #save info as a dataframe so we can append it to the dataframe that will store all examples
            new_row_df = pd.DataFrame([new_row])
            if n == 0: #if its the first example we make it the distance dataframe since we shouldnt concat. to an empty dataframe
                all_distances =  new_row_df
            else:
                #adding the example in the dataframe
                all_distances = pd.concat([all_distances, new_row_df], ignore_index=True)

            shutil.move(file_name , directory_path)
    file.close()
  all_distances.to_csv(f"{ID}_{method}_{gap_open:.5f}_{gap_extension:.5f}_distances.csv")
  shutil.move(f"{ID}_{method}_{gap_open:.5f}_{gap_extension:.5f}_distances.csv", os.path.join(os.getcwd(), f"{ID}_tests/{ID}_distances"))
  return


def main():
  
  example_id = sys.argv[1]  #'cd01068'
  parms = pd.read_csv('parms.csv')
  model, tokenizer = Ankh_initialize()
  example_pairs = str(example_id)+"_examples.csv"

  gap_open = parms.iloc[0][0]
  gap_extension = parms.iloc[0][1]

  #run all the exampels in the example_pairs file
  run_individual_method_tests(example_id, example_pairs, gap_open, gap_extension, model, tokenizer)


if __name__ == "__main__":
    main()