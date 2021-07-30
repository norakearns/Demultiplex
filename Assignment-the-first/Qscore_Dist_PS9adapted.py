#!/usr/bin/env python

import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description="A program to input coverage limit")
    parser.add_argument("-p", help="positions", required=True, type=int)
    parser.add_argument("-f", help="file name", required=True, type=str)
    parser.add_argument("-o", help="file name", required=True, type=str)
    parser.add_argument("-r", help="file name", required=True, type=str)
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

all_qscores = np.zeros(positions, dtype = int)
n = 0
record_num = 0
with gzip.open(filename, "rt") as fh:
    for line in fh:     #each line in the file
        n+=1            
        if n%4==0:      #only grab the Qscore line
            stripped = line.strip('\n') #strip off the newline 
            bp = 0
            for i in stripped:
                qscore = convert_phred(i)
                all_qscores[bp] += qscore
                bp += 1
            record_num += 1

all_qscores = all_qscores/record_num

print(all_qscores)


plt.xlabel('# Base Pair')
plt.ylabel('Mean Quality Score')
plt.title("Mean Quality Score of Each Base Position for " + read)
plt.bar(range(len(all_qscores)), all_qscores)
plt.savefig(output)
