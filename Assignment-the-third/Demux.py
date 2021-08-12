#!/usr/bin/env python


def convert_phred(letter):
#Takes as it's input each individual Phred score (character), and uses Phred33 encoding to convert it to an integer Quality Score.```
    QS = ord(letter) - 33
    return QS

def qual_score(phred_score):
#Takes a string of phred scores, converts them to QS number scores, and then finds the average QS of that string
    n=0
    for phred in phred_score:
        QS = convert_phred(phred)
        n+=QS
    return((n)/len(phred_score))

def rev_comp(sequence):
#Takes as it's input index2 and creates the reverse complement
    comp=""
    for i in sequence:
        if i == "A":
            comp += "T"
        elif i == "T":
            comp += "A"
        elif i == "G":
            comp += "C"
        elif i == "C":
            comp += "G"
        else:
            comp += "N"
    comp="".join(reversed(comp))
    return comp


def Index_Match(seq1, seq2):
#Takes as it's input index 1 and index 2, compares them, and if they are the same, returns True. If the indexes don't match, it returns False.
    if seq1 == seq2:
        return True
    else:
        return False
# # 48 files = one Read1 FASTQ file and one Read2 FASTQ file per matching index-pair (24 index pairs)
# #     - each read file needs to have the header adjusted to have the indexes on it (the read 1 and read 2 index)
# # 2 files = two FASTQ files with index-hopped reads-pairs
# #     - If the indexes in the header are not the same, bin it into the index-hopped reads-pairs file
# # 2 files = files undetermined (non-matching or low quality) index-pairs
# #     - The indexes don't match any of the 24 set indexes
# #     - Or if the indexes are low-quality (Ns in the index)
    
# # The algorithm should report:
# #     - the number of read-pairs with properly matched indexes (per index-pair) 
# #     - the number of read pairs with index-hopping observed 
# #     - the number of read pairs with unknown index(es). 


# <!-- 4. Pseudocode
# ```
# Open Read1, Index1, Index2, and Read2 - don't use with open
# FileNames:
# /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz  - 1452986940 lines
# /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz - 1452986940 linse
# /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
# /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz

import argparse
from os import write
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="A program to input Illumina files")
    parser.add_argument("-R1", help="Read1_file", required=True, type=str)
    parser.add_argument("-R2", help="Index1_file", required=True, type=str)
    parser.add_argument("-R3", help="Index2_file", required=True, type=str)
    parser.add_argument("-R4", help="Read2_file", required=True, type=str)
    return parser.parse_args()
     
args = get_args()

Read1 = args.R1 
Read2 = args.R2
Read3 = args.R3
Read4 = args.R4

#Read1 = "/projects/bgmp/nkearns/bioinformatics/Bi622/Demultiplex/Demultiplex/TEST-input_FASTQ/TEST_R1_input.fastq"

# Make a name for each file that makes sense:
Read1_file = gzip.open(Read1, "rt")
Read2_file = gzip.open(Read4, "rt")
Index1_file = gzip.open(Read2, "rt")
Index2_file = gzip.open(Read3, "rt")

#/projects/bgmp/nkearns/bioinformatics/Bi622/Demultiplex/Demultiplex/TEST-input_FASTQ

# Make a set of all the indexes
# sets are initialized with {}

Index_list = ["GTAGCGTA", "AACAGCGA", "CTCTGGAT", "CACTTCAC", "TATGGCAC", "TCGACAAG", "ATCGTGGT", "GATCTTGC", "CGATCGAT", "TAGCCATG","TACCGGAT", "GCTACTCT", "TGTTCCGT" , "TCTTCGAC", "TCGAGAGT", "AGAGTCCA", "GATCAAGG", "CGGTAATC", "CTAGCTCA", "ACGATCAG", "GTCCTAAG", "ATCATGCG", "TCGGATTC", "AGGATAGC"]

# Index_dict[key] = Index_R1 or Index_R2
# Index_dict[value] = file handle Index_R1.fastq or Index_R2.fastq
# to write to the file:
#   my_dict[Index_R1].write("Hi")

Index_dict = {} # initialize an empty dictionary
for i in Index_list: # loop through all the indexes
    Index_dict[str(i + "_R1")] = '' # set keys equal to Index_R1 or Index_R2
    Index_dict[str(i + "_R2")] = ''

for k in Index_dict: # loop through the Index dictionary
    #fh = '{}.fastq'.format(k) # create a file handle taking the key and adding ".fastq" to make the file handle. 
    fh = '{}.fastq'.format(k)
    #Index_dict[k] = open(fh, "w") # set the values of the dictioanry equal to all the file handles
    Index_dict[k] = open(fh, "wt")

R1_swapped = open("R1_swapped.fastq", "wt") #wt writes it out as text, not bytes
R2_swapped = open("R2_swapped.fastq", "wt")
R1_bad = open("R1_bad.fastq","wt")
R2_bad = open("R2_bad.fastq","wt")

num_lines = 0  # initialize a counter
for line in Read1_file: 
    num_lines +=1 # count the number of lines in the file
Read1_file.seek(0) # reset the file pointer to the start of the file

#     Create a sliding window of four lines so only the current record is held in memory as an array (record_array)
#     Use readline() to read in a line and store it as a variable, then store those variables in the array

R1_four = []    # [header, sequence, +, qscore_line]
R2_four = []    # [header, sequence, +, qscore_line]
I1_four = []    # [header, sequence, +, qscore_line]
I2_four = []    # [header, sequence, +, qscore_line]
readline_counter = 0
while readline_counter < num_lines: # for every line in the file
    R1_four = []    
    R2_four = []
    I1_four = []
    I2_four = [] 
    for j in range(4): # do the following for for lines
        Read1_line = Read1_file.readline().rstrip("\n") # temporarily stores the current line
        Read2_line = Read2_file.readline().rstrip("\n")
        Index1_line = Index1_file.readline().rstrip("\n")
        Index2_line = Index2_file.readline().rstrip("\n")
        Index2_line_RC = rev_comp(Index2_line)
        R1_four.append(Read1_line)
        R2_four.append(Read2_line)
        I1_four.append(Index1_line)
        I2_four.append(Index2_line_RC) # I2_four stores the REV COMP of the index, not the actual index
    readline_counter += 4    # increment by four to start at the next line of the file

    # print("at level of for loop \n") # if you do stuff here it processes every single record
# create new_header:
# new_header = current_header + total_index     
    Total_Index = I1_four[1] + "-" + I2_four[1]

# change the R1_four header to have Total_Index at the end
    R1_four[0] = R1_four[0] + " " + Total_Index 
    R2_four[0] = R2_four[0] + " " + Total_Index

# Check if indexes contain N by looking at the last 17 characters of the new header
# If Total_Index contains an N -> write the corresponding record from Read1 to R1_bad.out
    if "N" in R1_four[0][-18:]: 
        for i in R1_four:
            R1_bad.write(i + "\n")

    if "N" in R2_four[0][-18:]: 
        for i in R2_four:
            R2_bad.write(i + "\n")

# Check if Read Qscore is below the cutoff(30), or if Index quality line is below the cutoff(30)
    Read1_QS = qual_score(R1_four[3])
    Read2_QS = qual_score(R2_four[3])
    Index1_QS = qual_score(I1_four[3])
    Index2_QS = qual_score(I2_four[3])

    if Read1_QS < 20:
        for i in R1_four:
            R1_bad.write(i + "\n")
        
    if Read2_QS < 20:
        for i in R2_four:
            R2_bad.write(i + "\n")

    if Index1_QS < 30:
        for i in R1_four:
            R1_bad.write(i + "\n")
        
    if Index2_QS < 30:
        for i in R2_four:
            R2_bad.write(i + "\n")

    #print("Index 1 = " + str(I1_four[1]))
    #print("Index 2 = " + str(I2_four[1]))

    if Index_Match(I1_four[1], I2_four[1]) == False and (I1_four[1] in Index_list and I2_four[1] in Index_list):
        for i in R1_four:
            R1_swapped.write(i + "\n")
        for i in R2_four:
            R2_swapped.write(i + "\n")


    if Index_Match(I1_four[1], I2_four[1]) == True and (I1_four[1] in Index_list and I2_four[1] in Index_list): # The Rev comp will be in the list because it's the same as index1, which is in the list
        for i in R1_four:
            #R1_handle = Index_dict[str(I1_four[1] + "_R1")]
            Index_dict[str(I1_four[1] + "_R1")].write(i + "\n")
            #print(i)
            # <gzip _io.BufferedWriter name='GTAGCGTA_R1.fastq.gz' 0x2aaab2256d90>

        for i in R2_four:
            Index_dict[str(I1_four[1] + "_R2")].write(i + "\n")
            #print(i)
    elif (I1_four[1] not in Index_list) or (I2_four[1] not in Index_list):
        for i in R1_four:
            R1_bad.write(i + "\n")
        for i in R2_four:
            R2_bad.write(i + "\n")

#print("at the level of while loop \n") # if you do stuff here it only processes the last record

Read1_file.close()
Read2_file.close()
Index1_file.close()
Index2_file.close()

for k in Index_dict: # loop through the Index dictionary
    Index_dict[k].close()


#         Check if index Qscore is below the cutoff(20), or if Read quality line is below the cutoff(30):
#         Else if either index is below the cutoff (20) or read quality is below cutoff (30) > write the corresponding record from Read1 to R1_bad.out

#         Use Index_Match to check if indexes match:
#         If both indexes are in the known set of indexes: 
#             If Index1 is the same as (==) Index2 > return TRUE and write the corresponding record from Read1 to index#_R1.out
#             else if Index1 is different from (!=)  Index2 > return FALSE and write the corresponding record from Read1 to R1_swapped.out
#         Else (one or both of the indexes is not in the known set of indexes):
#             Write the corresponding record to R1_bad.out


#     For Read2:
#         create new_header:
#         new_header = current_header + total_index 

#         Check if indexes contain N:
#         If Index1 or Index2 contains an N > write the corresponding record from Read2 to R2_bad.out

#         Check if index Qscore is below the cutoff(20), or if Read quality line is below the cutoff(30):
#         Else if either index is below the cutoff (20) or read quality is below cutoff (30) > write the corresponding record from Read2 to R2_bad.out

#         Use Index_Match to check if indexes match:
#         If both indexes are in the known set of indexes: 
#             If Index1 is the same as (==) Index2 > return TRUE and write the corresponding record from Read2 to index#_R2.out
#             else if Index1 is different from (!=)  Index2 > return FLASE and write the corresponding record from Read1 to R2_swapped.out
#         Else (one or both of the indexes is not in the known set of indexes):
#             Write the corresponding record to R2_bad.out

# Close all the files!


