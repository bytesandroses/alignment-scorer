#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

# Importing relevant packages (Biopython, numpy, and os)
from Bio import AlignIO
import numpy as np
import os
from natsort import natsorted
import csv


# Provide the paths to the reference and realignment directories
refer_path = ("/home/isaac/Desktop/Test_Simulations/Results/A100_Samples/Reference_Alignments/")
reali_path = ("/home/isaac/Desktop/Test_Simulations/Results/A100_Samples/MAFFT_Alignments/")


# Create a list of the 100 simulations performed 
simulation_list = list(range(0, 100))


# Returns the path for the alignment file that will be analyzed
def Get_Alignment_File_Path(simulation_number, alignment_directory_path):
    alignment_files_list = os.listdir(alignment_directory_path)
    alignment_files_list = natsorted(alignment_files_list)
    alignment_filename = alignment_files_list[simulation_number]
    alignment_file_path = str(alignment_directory_path) + str(alignment_filename)
    return alignment_file_path


# Imports MSA file and sorts sequence IDs by the same order
def Open_MSA(file_path):
    alignment = AlignIO.read(file_path, "fasta")
    alignment.sort()
    return alignment


# Converts any lowercase nucleotides to uppercase
def Make_UpperCase(alignment):
    for i in range(100):
        sequence = alignment[i].seq
        sequence = sequence.upper()
        alignment[i].seq = sequence
    return alignment


# Calculates the number of columns in an alignment
def len_alignment_fn(alignment):
    len_alignment = len(alignment[0].seq)
    return len_alignment


# Represent MSA as an array of nucleotide positions
"""
For example, the alignment:
    
       [ , 0] [ , 1] [ , 2] [ , 3] [ , 4] [ , 5] ...... [ , 1574] [ , 1575]
[0, ]    A      T      G      A      -      T    ......     C        -
[1, ]    A      -      0      A      C      T    ......     -        T   
[2, ]    A      T      G      A      C      A    ......     -        -
[3, ]    -      -      -      G      C      A    ......     T        T
[4, ]    C      T      G      0      -      4    ......     T        T
[5, ]    A      -      -      G      T      A    ......     C        -
.....   ...    ...    ...    ...    ...    ...   ......    ...      ...
.....   ...    ...    ...    ...    ...    ...   ......    ...      ...
.....   ...    ...    ...    ...    ...    ...   ......    ...      ...
[99, ]   C      T      G      A      -      A    ......     C        T     
    

Would be represented as:
           
       [ , 0] [ , 1] [ , 2] [ , 3] [ , 4] [ , 5] ...... [ , 1574] [ , 1575]
[0, ]    1      2      3      4      0      5    ......    1500      0
[1, ]    1      0      0      2      3      4    ......     0       1500   
[2, ]    1      2      3      4      5      6    ......     0        0
[3, ]    0      0      0      1      2      3    ......    1500     1501
[4, ]    1      2      3      0      0      4    ......    1499     1500
[5, ]    1      0      0      2      3      5    ......    1500      0
.....   ...    ...    ...    ...    ...    ...   ......    ...      ...
.....   ...    ...    ...    ...    ...    ...   ......    ...      ...
.....   ...    ...    ...    ...    ...    ...   ......    ...      ...
[99, ]   1      2      3      4      0      5    ......    1499     1500 


In this array, rows [0, ] to [99, ] represent sequences 1-100, while columns 
[ , 0] to [ , 1574] represent the nucleotide positions. In the array, a value of 
0 represents a gap ("-") at that position. For example, the 0 at [0, 4] means 
that there is a gap in the 5th position/column of the first sequence. 

Values of 1 and greater means that the cell is occupied by a nucleotide, and 
not a gap. This number is also the gap free position. For example, [3, 3] has a 
value of 1, meaning that the 1st nucleotide in sequence 4 is in column 4 
(this is because the first 3 columns in this sequence are gaps).

The maximum value in each row is normally 1500, but it varies depending on how
many indels has experienced.

"""
def Get_Nucleotide_Positions(alignment, len_alignment):
    # Create a zero array with the same dimensions as the MSA
    reference_nucleotide_positions = np.zeros((100, len_alignment), dtype=int)
    
    # For loop through the MSA, assigning a 0 to cells with gaps, and 1-~1500 
    # for cells with nucleotides
    gap_count = int()
    for i in range(100):
        sequence = alignment[i].seq
        for j in range(len_alignment):
            residue = sequence[j]
            if residue == "-":
                reference_nucleotide_positions[i, j] = 0
            else:
                gap_count = (alignment[i].seq[0:j]).count("-")
                cell_column = j
                reference_nucleotide_positions[i, j] = cell_column + 1 - gap_count
    return reference_nucleotide_positions


# Calculates reference sequence nucleotide positions
def Get_Refer_Sequence_Nucleotide_Positions(alignment_positions, len_alignment):
    # Create a list, refer_seq_nuc_pos, to store the nucleotide positions of the 
    # reference sequence. Sequence ID 1 (alignment_positions[0, ]) will be used
    # as the reference
    refer_sequence = list(alignment_positions[0, ])
    refer_sequence_nuc_count = len_alignment - refer_sequence.count(0)
    refer_sequence_nuc_pos = [0] * refer_sequence_nuc_count
    
    j = 0
    for i in range (len_alignment):
        value = refer_sequence[i]
        if value > 0:
            refer_sequence_nuc_pos[j] = i
            j = j + 1
        else:
            value = value
    return refer_sequence_nuc_pos


# Calculate reference alignment score
def Get_Refer_Alignment_Column_Score(reference_alignment_nucleotide_positions, ref_seq_nuc_pos):
    refer_column_alignment_score = [0] * len(ref_seq_nuc_pos)
    for i in range(len(ref_seq_nuc_pos)):
        column = ref_seq_nuc_pos[i]
        column_gap_count = (list(reference_alignment_nucleotide_positions[:, column])).count(0)
        column_score = 100 - column_gap_count
        refer_column_alignment_score[i] = column_score
    return refer_column_alignment_score


# Score realignment sequence columns against reference sequence columns
def Get_Realignment_Column_Scores(reference_columns, realignment_columns, ref_seq_nuc_pos, reali_seq_nuc_pos):
    reali_col_score = [0] * len(ref_seq_nuc_pos) 
    for i in range(len(ref_seq_nuc_pos)):
        ref_col_seq = reference_columns[:, i] 
        reali_col_seq = realignment_columns[:, i]
        for j in range(100):
            ref_value = ref_col_seq[j]
            reali_value = reali_col_seq[j]
            if ref_value == 0:
                ref_value = ref_value
            elif ref_value > 0:
                if ref_value == reali_value:
                    score = reali_col_score[i] + 1
                    reali_col_score[i] = score
                elif ref_value != reali_value:
                    if reali_value == 0:
                        score = reali_col_score[i] - 5
                        reali_col_score[i] = score
                    elif reali_value != 0:
                        score = reali_col_score[i] - 1
                        reali_col_score[i] = score
    return reali_col_score
        

# Create a list to store the alignment score difference for all 100 alignments
alignment_score_difference = [0] * 100


# Create a csv to save the alignment scores into
with open('Alignment Score Results.csv', 'a', newline = '') as csvfile:
    fieldnames = ['Reference_Score', 'Realignment_Score', 'Alignment_Score_Difference']
    alignmentwriter = csv.DictWriter(csvfile, fieldnames = fieldnames)
    alignmentwriter.writeheader()


# Calculates the alignment score of reference alignment
def Alignment_Score_Differenence_Calculator(simulation_number):
    # Upload the reference and realignment files that will be analyzed
    refer_file_path = Get_Alignment_File_Path(simulation_number, refer_path)
    refer_alignment = Open_MSA(refer_file_path)
    refer_alignment_upper = Make_UpperCase(refer_alignment)
    
    reali_file_path = Get_Alignment_File_Path(simulation_number, reali_path)
    reali_alignment = Open_MSA(reali_file_path)
    reali_alignment_upper = Make_UpperCase(reali_alignment)
    
    
    # Get the lengths of the reference and realignment alignments
    len_refer = len_alignment_fn(refer_alignment_upper)
    len_reali = len_alignment_fn(reali_alignment_upper)
    
    
    # Get the column numbers of the reference sequence's nucleotide positions
    # (refer_sequence_positions), the nucleotide composition of that entire 
    # column (refer_columns), the score for each of those columns 
    # (refer_alignment_column_scores), and the reference alignment score (refer_alignment_score)
    refer_alignment_positions = Get_Nucleotide_Positions(refer_alignment_upper, len_alignment = len_refer)
    refer_sequence_positions = Get_Refer_Sequence_Nucleotide_Positions(refer_alignment_positions, len_refer)
    refer_columns = refer_alignment_positions[:, list(refer_sequence_positions)]
    refer_alignment_column_scores = Get_Refer_Alignment_Column_Score(refer_alignment_positions, refer_sequence_positions)
    
    
    # Get the column numbers of the reference sequence's nucleotide positions
    # within the realignment alignment (reali_sequence_positions), and the nucleotide 
    # composition of those column (reali_columns), the score for each of those columns 
    # (refer_alignment_column_scores)
    reali_alignment_positions = Get_Nucleotide_Positions(reali_alignment_upper, len_alignment = len_reali)
    reali_sequence_positions = Get_Refer_Sequence_Nucleotide_Positions(reali_alignment_positions, len_reali)
    reali_columns = reali_alignment_positions[:, list(reali_sequence_positions)]
    
    
    # Get the reference alignment score, the individual realignment column scores 
    # (reali_alignment_column_scores), the realignment score, and the alignment 
    # score difference (alignment_score_difference)
    reali_alignment_column_scores = Get_Realignment_Column_Scores(refer_columns, reali_columns, refer_sequence_positions, reali_sequence_positions)
    reference_alignment_score = sum(refer_alignment_column_scores)
    realignment_alignment_score = sum(reali_alignment_column_scores)
    alignment_score_difference[simulation_number] = reference_alignment_score - realignment_alignment_score
    print(alignment_score_difference[simulation_number])
    
    
    # Save the alignment score results to a csv file
    with open('Alignment Score Results.csv', 'a', newline = '') as csvfile:
        alignmentwriter = csv.writer(csvfile, delimiter = ',',
                                     quotechar = '|', quoting=csv.QUOTE_MINIMAL)
        alignmentwriter.writerow([reference_alignment_score, realignment_alignment_score, alignment_score_difference[simulation_number]])


for simulation_number in simulation_list: 
    Alignment_Score_Differenence_Calculator(simulation_number)



