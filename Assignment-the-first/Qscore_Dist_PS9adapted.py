#!/usr/bin/env python

import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description="A program to input coverage limit")
    parser.add_argument("-p", help="positions", required=True, type=int) #the number of nucleotides in the read. 101 for reads, 8 for indexes.
    parser.add_argument("-f", help="file name", required=True, type=str) #the path to the file on talapas
    parser.add_argument("-o", help="output file name", required=True, type=str) #the name of the output histogram file: R1_hist.png, R2_hist.png, R3_hist.png, R4_hist.png
    parser.add_argument("-r", help="Read or Index", required=True, type=str) #Read1, Read2, Index1, or Index2
    return parser.parse_args()
     
args = get_args()
print(args.f, args.p, args.o, args.r)

filename = args.f 
positions = args.p
output = args.o
read = args.r

def convert_phred(letter):
    QS = ord(letter) - 33
    return QS

all_qscores = np.zeros(positions, dtype = int) #create an array of 0s
n = 0
record_num = 0
with gzip.open(filename, "rt") as fh: #use gzip to open the file because it is HUGE
    for line in fh:     #each line in the file
        n+=1            #increment the line
        if n%4==0:      #only grab the Qscore line
            stripped = line.strip('\n') #strip off the newline 
            bp = 0      #set the position of the base pair to 0
            for i in stripped:  #for each phred score in the quality line
                qscore = convert_phred(i) #convert it to a quality score
                all_qscores[bp] += qscore  #add on that value to the value in the correponsing position in the all_qscores array (take a running total of the qscores at a given base position)
                bp += 1    #increment the position of the base pair by 1
            record_num += 1 #incremet the number of records by 1

all_qscores = all_qscores/record_num #take the mean of each qscores by dividing the running total by the number of records.

print(all_qscores)


plt.xlabel('# Base Pair')
plt.ylabel('Mean Quality Score')
plt.title("Mean Quality Score of Each Base Position for " + read)
plt.bar(range(len(all_qscores)), all_qscores)
plt.savefig(output)
