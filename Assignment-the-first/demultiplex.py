#!/usr/bin/env python

import argparse
def get_args():
    parser = argparse.ArgumentParser(description="demultiplex")
    parser.add_argument("-l", help="length", type=int)
    parser.add_argument("-f", help="specify the filename")
    parser.add_argument("-o")
    return parser.parse_args()

args=get_args()
# 1294_S1_L008_R1_001.fastq.gz
# 1294_S1_L008_R2_001.fastq.gz
# 1294_S1_L008_R3_001.fastq.gz
# 1294_S1_L008_R4_001.fastq.gz
import numpy as np
import bioinfo
import gzip
lane1 = np.zeros(args.l)
with gzip.open(args.f, "rt") as fh:
    #open fastq file
    i = 0 
    for line in fh:
        #loop record
        i+=1

        line = line.strip('\n')
        #print(my_list)
        if i%4 == 0:
            #for character position in each line, sum the converted phred scores:
            #counter_array= (i//4) #this line addition for using in numpy since we cant append
            for ind, letter in enumerate(line):
                #index each letter in the line, resulting in an index and the value
                score=bioinfo.convert_phred(letter)
                lane1[ind]+=score #this is instead of list.append used in lists
        # if i % 500000 == 0:
            # print("working on line", i)
#print(lane1)
mean = np.zeros(args.l,float)
for k, score in enumerate(lane1):
    #for each score in lane1, index and run:
    mean[k] = score/(i/4)


import matplotlib.pyplot as plt
x=range(args.l)
y=mean
fig, ax=plt.subplots()
ax.bar(x,y)
#only need to specify x axis3
plt.xlabel('Read Position')
plt.ylabel('Mean')
plt.title('Mean Quality Score')
plt.savefig(args.o+"png")