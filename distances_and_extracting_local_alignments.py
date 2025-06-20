#1. required imports
import os
import sys
import warnings
import re
from Bio import SeqIO
from scipy import stats
import esm
import numpy as np
import pandas as pd
import torch
import shutil
from tqdm import tqdm


#2. extrating the computed alignment relative to the reference
#given two proteins (protein1, protein2) that were locally aligned, the computed alignment is [algn1, algn2]
#if the reference local alignment between protein1 and protein2 is [ref1, ref2], we can "extract" an alignment
#based on the computed alignment is [algn1, algn2] (details found in our paper), so the distances between the 2 
#alignments (reference, and computed) can be calculated

#helper functions:
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


#main function for extrating an alignment relative to a reference
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
    

#3. DISTANCES:
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


#example
def main():

    protein1 = "METDLNSQDRKDLDKFIKFFALKTVQVIVQARLGEKICTRSSSSPTGSDWFNLAIKDIPEVTHEAKKALAGQLPAVGRSMCVEISLKTSEGDSMELEIWCLEMNEKCDKEIKVSYTVYNRLSLLLKSLLAITRVTPMQEQATSSIAASSLPSSSERSSSSALHHELKEGMESDDEIRRVPEMGGEATGTTSASGRDGVSAAGQAQPSAGTQRKRGRSPADKENKRLKRLLRNRVSAQQARERKKAYLIDLEARVKELETKNAELEERLSTLQNENQMLRHILKNTTAGAQEGRKAYRLSRKQGHEYVILYRIYFGEVQLSGLGEGFQTVRVGTVGTPVGTITLSCAYRINLAFMST"
    protein2 = "MAYQLYRNTTLGNSLQESLDELIQSQQITPQLALQVLLQFDKAINAALAQRVRNRVNFMAAQEQEQEKQQVKTSTTSSLPSSSERSSSSAPNNLKEGGGVESDEEIRRVPEMGGGGGSASSGAGADERQGKEDGKQQGGGGGGAAAAGGGQEQAPPARKRGRSAGDKEQNRLKRLLRNRVSAQQARERKKAYMTELEAKAKDLELRNAELEQRVSTLQNENNTLRQILKNTTAHAGKRGGGGGGKGGDGGGGGKKHHFTKSRGSLNTYRFCDNVWTFVLNDVEFREVTELIKVDKVKIVACDGKNTGSNTTE"

    reference_alignment_seq1 = "MQEQATSSIAASSLPSSSERSSSSAL-HHELKEGMESDDEIRRVPEMGGEATG------TTSASGRDGVSAAGQ--AQPSAG--------TQRKRGRSPADKENKRLKRLLRNRVSAQQARERKKAYLIDLEARVKELETKNAELEERLSTLQNENQMLRHILKNTTAGA--QEGRK"
    reference_alignment_seq2 = "QEKQQVKTSTTSSLPSSSERSSSSAPNNLKEGGGVESDEEIRRVPEMGGGGGSASSGAGADERQGKEDGKQQGGGGGGAAAAGGGQEQAPPARKRGRSAGDKEQNRLKRLLRNRVSAQQARERKKAYMTELEAKAKDLELRNAELEQRVSTLQNENNTLRQILKNTTAHAGKRGGGG"

    computed_alignment_seq1 = "VEISLKTSEGDSMELEIWCLEMNEKCDKEIKVSYTVYNRLSLLLKSLLAIT----------RVTPMQEQATSSIAASSLPSSSERSSSSALHH-ELKEGMESDDEIRRVPEMGGEATGTTSASGRD----------------GVSAAGQAQPSAGTQRKRGRSPADKENKRLKRLLRNRVSAQQARERKKAYLIDLEARVKELETKNAELEERLSTLQNENQMLRHILKNTTAGAQEGR------------KAYRLSRKQGHEYVILYRIYFGEVQLSGLGEGFQTVRVGTVGTPVGTIT"
    computed_alignment_seq2 = "MAYQLYRNTTLGNSLQESLDELIQSQQITPQLALQVLLQFDKAINAALAQRVRNRVNFMAAQEQEQEKQQVKTSTTSSLPSSSERSSSSAPNNLKEGGGVESDEEIRRVPEMGGGGGSASSGAGADERQGKEDGKQQGGGGGGAAAAGGGQEQAPPARKRGRSAGDKEQNRLKRLLRNRVSAQQARERKKAYMTELEAKAKDLELRNAELEQRVSTLQNENNTLRQILKNTTAHAGKRGGGGGGKGGDGGGGGKKHHFTKSRGSLNTYRFCDNVWTFVLNDVEFREVTELIKVDKVK--I"

    #extract computed alignment so the distances can be used
    algn1_final, algn2_final = extract_alignment(protein1, protein2, reference_alignment_seq1, reference_alignment_seq2, computed_alignment_seq1, computed_alignment_seq2)
    print(algn1_final)
    print(algn2_final)
    print("\n")

    #get distances between reference alignment and extracted alignment
    alignment_extracted = [algn1_final, algn2_final]
    reference_alignment = [reference_alignment_seq1, reference_alignment_seq2]

    d_SSP, d_SEQ , d_POS = get_3_dists(alignment_extracted, reference_alignment)
    d_CC = d_cc(alignment_extracted, reference_alignment)
    d_D = d_d(alignment_extracted, reference_alignment)

    print("SSP: " + str(d_SSP))
    print("SEQ: " + str(d_SEQ))
    print("POS: " + str(d_POS))
    print("CC: " + str(d_CC))
    print("D: " + str(d_D))




if __name__ == "__main__":
    main()
