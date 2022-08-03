#!/usr/bin/env python
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="demultiplex")
    #parser.add_argument("-l", help="length", type=int)
    parser.add_argument("-read1", help="specify the filename")
    parser.add_argument("-index1", help="specify the filename")
    parser.add_argument("-index2", help="specify the filename")
    parser.add_argument("-read2", help="specify the filename")
    parser.add_argument("-known_indexes", help="specify the filename")
    #parser.add_argument("-o")
    return parser.parse_args()

args=get_args()
# 1294_S1_L008_R1_001.fastq.gz
# 1294_S1_L008_R2_001.fastq.gz
# 1294_S1_L008_R3_001.fastq.gz
# 1294_S1_L008_R4_001.fastq.gz
import numpy as np
import bioinfo
import gzip
known_indexes=[]
with open(args.known_indexes, "r") as fh:
    for line in fh:
        #print(line)
        known_indexes.append(line.strip())
R1_matched={}
R2_matched={}
for index in known_indexes:
    R1_matched[index]=open(f"./outputpractice/{index}_matched.fastq", "w")
    R2_matched[index]=open(f"./outputpractice/{index}_matched.fastq", "w")
#print(known_indexes)
read1fh=open(args.read1, "r")
index1fh=open(args.index1, "r")
index2fh=open(args.index2, "r")
read2fh=open(args.read2, "r")
unknown_R1=open("./outputpractice/unknown_R1.fastq", "w")
unknown_R2=open("./outputpractice/unknown_R2.fastq", "w")
unmatched_R1=open("./outputpractice/hopped_R1.fastq", "w")
unmatched_R2=open("./outputpractice/hopped_R2.fastq", "w")
unknown_counter=0
unmatched_counter=0
matched_counter=0
matched_dict=[]
while True:
    r1_record=[read1fh.readline().strip(), read1fh.readline().strip(), read1fh.readline().strip(), read1fh.readline().strip()]
    r2_record=[read2fh.readline().strip(), read2fh.readline().strip(), read2fh.readline().strip(), read2fh.readline().strip()]
    i1_record=[index1fh.readline().strip(), index1fh.readline().strip(), index1fh.readline().strip(), index1fh.readline().strip()]
    i2_record=[index2fh.readline().strip(), index2fh.readline().strip(), index2fh.readline().strip(), index2fh.readline().strip()]
    if(r1_record[0] == ""):
        break
    index1=i1_record[1]
    R1_header=r1_record[0]
    R2_header=r2_record[0]
    index2=i2_record[1]
    if "N" in {index1} or "N" in {index2}:
        index2_rev=bioinfo.reverse_complement(index2)
        R1_header=R1_header+'_'+index1+'_'+index2_rev
        R2_header=R2_header+'_'+index1+'_'+index2_rev
    #print(R1_header)
        unknown_R1.write(R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3]+'\n')
        unknown_R2.write(R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3]+'\n')
        unknown_counter+=1
#print(unknown_counter)
    else:
        if index1 not in known_indexes or index2 not in known_indexes:
            index2_rev=bioinfo.reverse_complement(index2)
            R1_header=R1_header+'_'+index1+'_'+index2_rev
            R2_header=R2_header+'_'+index1+'_'+index2_rev
            unknown_R1.write(R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3]+'\n')
            unknown_R2.write(R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3]+'\n')
            unknown_counter+=1
        else:
            written = False
            for value in i1_record[1]:
                convert_value=bioinfo.convert_phred(value)
                print(convert_value)
                if convert_value < 30:
                    index2_rev=bioinfo.reverse_complement(index2)
                    R1_header=R1_header+'_'+index1+'_'+index2_rev
                    R2_header=R2_header+'_'+index1+'_'+index2_rev
                    unknown_R1.write(R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3]+'\n')
                    unknown_R2.write(R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3]+'\n')
                    unknown_counter+=1
                    written=True
                    break
            if not written:
                for value in i2_record[1]:
                    convert_value=bioinfo.convert_phred(value)
                    if convert_value < 30:
                        index2_rev=bioinfo.reverse_complement(index2)
                        R1_header=R1_header+'_'+index1+'_'+index2_rev
                        R2_header=R2_header+'_'+index1+'_'+index2_rev
                        unknown_R1.write(R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3]+'\n')
                        unknown_R2.write(R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3]+'\n')
                        unknown_counter+=1
                    else:
                        index2_rev=bioinfo.reverse_complement(index2)
                        if index2_rev != index1:
                            R1_header=R1_header+'_'+index1+'_'+index2_rev
                            R2_header=R2_header+'_'+index1+'_'+index2_rev
                            unmatched_R1.write(R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3]+'\n')
                            unmatched_R2.write(R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3]+'\n')
                            unmatched_counter+=1
                        else:
                            R1_header=R1_header+'_'+index1+'_'+index2_rev
                            R2_header=R2_header+'_'+index1+'_'+index2_rev
                            R1_matched[index1]= (R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3])
                            R2_matched[index2_rev] = (R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3])
                            matched_counter+=1        
print(unknown_counter)
print(unmatched_counter)
print(matched_counter)

                        

# create a set of known indexes from the R1 file
# with all 4 files open
#     for the nth record, create temporary arrays ([header, seq, +, qualityscore]) for each of the 4 files (result=four arrays with header in position 1, seq in position 2, + in position 3, qualscore position 4)
#         is there an "n" present in index1 or index2 of the nth record?
#             if yes: send record from original read file to unknown file designated for that read (1 or 2), create reverse compliment variable and add to header
#                 unknown_counter+=1
#             if no:for qscore line of nth record, are any of the phred values in qscore line < 30?
#                 if yes: send record from original read file to unknown file designated for that read (1 or 2), create reverse compliment variable and add to header
#                     unknown_counter+=1
#                 if no: is index 1 in set of known indexes?
#                     if not: write record from original read file to unknown file designated for that read, create reverse_compliment variable and add to header
#                         unknown_counter+=1
#                     if yes: create reverse complement of index2, is that reverse complement in the set?
#                         if not: write to unknown file for that read, add reverse compliment variable to header
#                             unknown_counter+=1
#                         if yes: is reverse compliment of index2 equal to that read's index1?
#                             if yes:write to matched file for that read, add reverse compliment to header
#                                 matched_counter+=1
#                             if not:write to unmatched file for that read, add reverse compliment to header
#                                 unmatched_counter+=1