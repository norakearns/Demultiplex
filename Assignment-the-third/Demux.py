#!/usr/bin/env python

# Function Definitions
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

# # DEMUX Algorithm Description
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

# FileNames:
# /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz  - 1452986940 lines
# /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz - 1452986940 lines
# /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
# /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz

import argparse
from os import write
import gzip

def get_args(): # generalize code using argparse
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

# Make a name for each file that makes sense:
Read1_file = gzip.open(Read1, "rt") # write to the file with text, not binary
Read2_file = gzip.open(Read4, "rt")
Index1_file = gzip.open(Read2, "rt")
Index2_file = gzip.open(Read3, "rt")


# Make a list of all the indexes
Index_list = ["GTAGCGTA", "AACAGCGA", "CTCTGGAT", "CACTTCAC", "TATGGCAC", "TCGACAAG", "ATCGTGGT", "GATCTTGC", "CGATCGAT", "TAGCCATG","TACCGGAT", "GCTACTCT", "TGTTCCGT" , "TCTTCGAC", "TCGAGAGT", "AGAGTCCA", "GATCAAGG", "CGGTAATC", "CTAGCTCA", "ACGATCAG", "GTCCTAAG", "ATCATGCG", "TCGGATTC", "AGGATAGC"]

# Index_dict[key] = Index_R1 or Index_R2
# Index_dict[value] = file handle Index_R1.fastq or Index_R2.fastq
# to write to the file: my_dict[Index_R1].write("Hi")

Index_dict = {} # initialize an empty dictionary
for i in Index_list: # loop through all the indexes
    Index_dict[str(i + "_R1")] = '' # set keys equal to Index_R1 or Index_R2
    Index_dict[str(i + "_R2")] = ''

for k in Index_dict: # loop through the Index dictionary
    # create a file handle taking the key and adding ".fastq" to make the file handle. This line is for test file
    fh = '{}.fastq'.format(k) 
    # set the values of the dictioanry equal to all the file handles and open them to write text
    Index_dict[k] = open(fh, "wt")

R1_swapped = open("R1_swapped.fastq", "wt") # Swapped reads file
R2_swapped = open("R2_swapped.fastq", "wt")
R1_bad = open("R1_bad.fastq","wt") # Low quality read, or reads containing Ns, file
R2_bad = open("R2_bad.fastq","wt")
Demux_algorithm_report = open("Demultiplex_algorithm_report.out","wt") # Output report file
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
# initialize the counters that are going to hold statistics that will be written to the output report    
N_count = 0
low_qual = 0
swapped_count = 0
correct_index = 0
unknown_index = 0
readline_counter = 0
while readline_counter < num_lines: # for every line in the file
    R1_four = []    # clear the four line array before moving on to the next record
    R2_four = []
    I1_four = []
    I2_four = [] 
    for j in range(4): # do the following for for lines
        Read1_line = Read1_file.readline().rstrip("\n") # temporarily stores the current line
        Read2_line = Read2_file.readline().rstrip("\n")
        Index1_line = Index1_file.readline().rstrip("\n")
        Index2_line = Index2_file.readline().rstrip("\n")
        Index2_line_RC = rev_comp(Index2_line)
        R1_four.append(Read1_line) # adds the line to the R1_four array
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
# Find the quality scores of all the reads
    Read1_QS = qual_score(R1_four[3])
    Read2_QS = qual_score(R2_four[3])
    Index1_QS = qual_score(I1_four[3])
    Index2_QS = qual_score(I2_four[3])

# Check if indexes contain N by looking at the last 17 characters of the new header
# If Total_Index contains an N -> write the corresponding record from Read1 to R1_bad.out
    if "N" in R1_four[0][-18:] or "N" in R2_four[0][-18:]:
        N_count += 1 
        for i in R1_four:
            R1_bad.write(i + "\n")
        for i in R2_four:
            R2_bad.write(i + "\n")    

# Check if Read Qscore is below the cutoff(30), or if Index quality line is below the cutoff(30). 
# If so, write to the R1 or R2_bad file
    elif Read1_QS < 20 or Index1_QS < 30 or Read2_QS < 20 or Index2_QS < 30:
        low_qual += 1
        for i in R1_four:
            R1_bad.write(i + "\n")
        for i in R2_four:
            R2_bad.write(i + "\n")
# Check it the Indexes do not match. If Index_Match returns False, the indexes don't match.
# If Index_Match returns False, write the records to the swapped files.
    elif Index_Match(I1_four[1], I2_four[1]) == False and (I1_four[1] in Index_list and I2_four[1] in Index_list):
        swapped_count += 1
        for i in R1_four:
            R1_swapped.write(i + "\n")
        for i in R2_four:
            R2_swapped.write(i + "\n")

# Check it the Indexes do match. If Index_Match returns True, the indexes do match.
# If Index_Match returns False, write the records to the swapped files.
    elif Index_Match(I1_four[1], I2_four[1]) == True and (I1_four[1] in Index_list and I2_four[1] in Index_list): # The Rev comp will be in the list because it's the same as index1, which is in the list
        correct_index += 1
        for i in R1_four:
            #R1_handle = Index_dict[str(I1_four[1] + "_R1")]
            Index_dict[str(I1_four[1] + "_R1")].write(i + "\n")
            #print(i)
            # <gzip _io.BufferedWriter name='GTAGCGTA_R1.fastq.gz' 0x2aaab2256d90>
        for i in R2_four:
            Index_dict[str(I1_four[1] + "_R2")].write(i + "\n")
            #print(i)
 # Catch-all which grabs all the reads with unknown indexes at the end.   
    else:
        unknown_index += 1
        for i in R1_four:
            R1_bad.write(i + "\n")
        for i in R2_four:
            R2_bad.write(i + "\n") 

Demux_algorithm_report.write("Reads with Ns = " + str(N_count) + "\n")
Demux_algorithm_report.write("Reads with Low Quality = " + str(low_qual) + "\n")
Demux_algorithm_report.write("Reads with Swapped Index = " + str(swapped_count) + "\n")
Demux_algorithm_report.write("Reads with Correct Index = " + str(correct_index)+ "\n")
Demux_algorithm_report.write("Reads with Unknown Index = " + str(unknown_index)+ "\n")

#print("at the level of while loop \n") # if you do stuff here it only processes the last record

# Close all the files!
Read1_file.close()
Read2_file.close()
Index1_file.close()
Index2_file.close()
Demux_algorithm_report.close()

for k in Index_dict: # loop through the Index dictionary
    Index_dict[k].close()



