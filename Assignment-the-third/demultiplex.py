#!/usr/bin/env python
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="demultiplex")
    parser.add_argument("-output", help="specify the filename")
    parser.add_argument("-read1", help="specify the filename")
    parser.add_argument("-index1", help="specify the filename")
    parser.add_argument("-index2", help="specify the filename")
    parser.add_argument("-read2", help="specify the filename")
    parser.add_argument("-known_indexes", help="specify the filename")
    return parser.parse_args()

args=get_args()
import numpy as np
import bioinfo
import gzip
import itertools
index_hopped_dict={}
index_matched_dict={}
known_indexes=[]
with open(args.known_indexes, "r") as fh:
    for line in fh:
        #print(line)
        known_indexes.append(line.strip())
# index_perm = itertools.permutations(known_indexes,2)
# for index1, index2 in index_perm:
#     key = index1+'_'+index2
#     index_hopped_dict[key]=0
index_prod= itertools.product(known_indexes, repeat=2)
for index1, index2 in index_prod:
    key = index1+'_'+index2
    index_matched_dict[key]=0
#print(index_hopped_dict)
#print(index_matched_dict)
R1_matched={}
R2_matched={}
for index in known_indexes:
    R1_matched[index]=gzip.open(f"./{args.output}/{index}_matched_R1.fastq.gz", "wt")
    R2_matched[index]=gzip.open(f"./{args.output}/{index}_matched_R2.fastq.gz", "wt")
#print(known_indexes)
read1fh=gzip.open(args.read1, "rt")
index1fh=gzip.open(args.index1, "rt")
index2fh=gzip.open(args.index2, "rt")
read2fh=gzip.open(args.read2, "rt")
unknown_R1=gzip.open(f"./{args.output}/unknown_R1.fastq.gz", "wt")
unknown_R2=gzip.open(f"./{args.output}/unknown_R2.fastq.gz", "wt")
unmatched_R1=gzip.open(f"./{args.output}/hopped_R1.fastq.gz", "wt")
unmatched_R2=gzip.open(f"./{args.output}/hopped_R2.fastq.gz", "wt")
outputstats=open(f"./{args.output}/stats.txt", "wt")
unknown_counter=0
hopped_counter=0
matched_counter=0
lowquality=0
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
    if "N" in index1 or "N" in index2:
        index2_rev=bioinfo.reverse_complement(index2)
        #print("trigger n")
        R1_header=R1_header+'_'+index1+'_'+index2_rev
        R2_header=R2_header+'_'+index1+'_'+index2_rev
    #print(R1_header)
        unknown_R1.write(R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3]+'\n')
        unknown_R2.write(R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3]+'\n')
        unknown_counter+=1
#print(unknown_counter)
    else:
        index2_rev=bioinfo.reverse_complement(index2)
        if index1 not in known_indexes or index2_rev not in known_indexes:
            #print("index1 or index2 unknown")
            R1_header=R1_header+'_'+index1+'_'+index2_rev
            R2_header=R2_header+'_'+index1+'_'+index2_rev
            unknown_R1.write(R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3]+'\n')
            unknown_R2.write(R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3]+'\n')
            unknown_counter+=1
        else:
            written = False
            for value in i1_record[3]:
                convert_value=bioinfo.convert_phred(value)
                #print(convert_value)
                if convert_value < 30:
                    #print("index1 lowquality")
                    R1_header=R1_header+'_'+index1+'_'+index2_rev
                    R2_header=R2_header+'_'+index1+'_'+index2_rev
                    unknown_R1.write(R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3]+'\n')
                    unknown_R2.write(R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3]+'\n')
                    unknown_counter+=1
                    lowquality+=1
                    written=True
                    break
            if not written:
                for value in i2_record[3]:
                    convert_value=bioinfo.convert_phred(value)
                    if convert_value < 30:
                        #print("index2 low quality")
                        R1_header=R1_header+'_'+index1+'_'+index2_rev
                        R2_header=R2_header+'_'+index1+'_'+index2_rev
                        unknown_R1.write(R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3]+'\n')
                        unknown_R2.write(R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3]+'\n')
                        unknown_counter+=1
                        lowquality+=1
                        written=True
                        break
                if not written:
                    if index2_rev != index1:
                        #print("hopped")
                        R1_header=R1_header+'_'+index1+'_'+index2_rev
                        R2_header=R2_header+'_'+index1+'_'+index2_rev
                        unmatched_R1.write(R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3]+'\n')
                        unmatched_R2.write(R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3]+'\n')
                        hopped_counter+=1
                        index_matched_dict[index1+'_'+index2_rev]+=1
                    else:
                        #print("matched:", index1)
                        R1_header=R1_header+'_'+index1+'_'+index2_rev
                        R2_header=R2_header+'_'+index1+'_'+index2_rev
                        R1_matched[index1].write(R1_header+'\n'+r1_record[1]+'\n'+r1_record[2]+'\n'+r1_record[3])
                        R2_matched[index2_rev].write(R2_header+'\n'+r2_record[1]+'\n'+r2_record[2]+'\n'+r2_record[3])
                        matched_counter+=1 
                        index_matched_dict[index1+'_'+index2_rev]+=1       
outputstats.write("unknown"+"\t"+ "hopped"+"\t"+ "matched"+"\t"+"low_qual"+"\n")
outputstats.write(str(unknown_counter)+'\t'+str(hopped_counter)+'\t'+str(matched_counter)+'\t'+str(lowquality)+'\n')
# print(unknown_counter)
# print(hopped_counter)
# print(matched_counter)
# print(lowquality)
# outputstats.write('Hopped Indexes')
# for key, value in index_hopped_dict.items():
#     print(key, value)
#     outputstats.write(str(key)+'\t'+str(value)+'\n')
outputstats.write('Indexes'+"\n")
for key, value in index_matched_dict.items():
    print(key, value, str(round(value/(matched_counter+hopped_counter+unknown_counter)*100,2))+"%",sep="\t")
    outputstats.write(str(key)+'\t'+str(value)+"\t"+str(round(value/(matched_counter+hopped_counter+unknown_counter)*100,2))+"%"+'\n')
for index in known_indexes:
        unknown_R1.close()
        unknown_R2.close()
read1fh.close()
index1fh.close()
index2fh.close()
read2fh.close()
unknown_R1.close()
unknown_R2.close()
unmatched_R1.close()
unmatched_R2.close()
outputstats.close()
#for test files first record has N in index1 (unknown), second record has matching indexes, third record hopped (but matching known so it makes it through loop) 