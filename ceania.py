# -*- coding: utf-8 -*-
"""
@author: earlm.lanl.gov
"""
import os
import concurrent.futures
import threading
import sys
import time

"""
USAGE:
python3 ceania.py ~/PATH/to/Multifasta_alignment.fa #CPUs_to_use > ~/PATH/to/Multifasta_alignment.ani
"""

FILE = sys.argv[1]
threads = int(sys.argv[2])
FILE_base = os.path.splitext(FILE)[0]
out_FILE=FILE_base+".ANI"
lock = threading.Lock()


#imports multifasta (aligned or not) into a dict with header/name as key
def import_fasta(FILE):
    with open(FILE, newline='') as multifasta:
        seq_dict = {}
        seq_names=[]
        for line in multifasta:
            if line[0] == ">" :
                 seq_name = line.replace(">", "").replace("\n","")
                 seq_dict[seq_name] = ""
                 seq_names.append(seq_name)
            else:
                seq_dict[seq_name]=seq_dict[seq_name] + line.replace("\n","")
    return seq_dict

# function to find the ANI of two sequences (seq1/2 = names and seq1/2_seq = sequence)
def ani(seq1,seq2,seq1_seq,seq2_seq):
    I = 0 # total length of alignment
    J = 0 # match counter
    k = 0 # no-gap alignment length counter
    while I < len(seq1_seq):
        if seq1_seq[I] == "-" or seq2_seq[I] == "-":
            None
        elif seq1_seq[I] == seq2_seq[I]:
            J+=1
            k+=1
        else:
            k+=1
        I+=1
    no_gaps=J/float(k)
    gaps=J/float(I)
    return [seq1, seq2, no_gaps, gaps]

if __name__ == '__main__':
    # import aligned multifasta
    seq_dict = import_fasta(FILE)

    # Make a list of all the comparisons we would like to make (remove redundancy)
    ANI_list=[]
    for I in seq_dict.keys():
        for J in seq_dict.keys():
            # remove self comparisons
            # and redundant half of comparisons (by only keeping sorted order...I know, smart)
            if [I,J] == sorted([I,J]) and I != J:
                ANI_list.append([I,J])
    # calculate ANI and write to file as seq1\tseq2\tANI
    print("seq1\tseq2\tANI_no_gap\tANI_gaps")
    ANI_gaps=[]
    ANI_no_gaps=[]
    start = time.perf_counter()
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        results = [executor.submit(ani, I[0],I[1],seq_dict[I[0]],seq_dict[I[1]]) for I in ANI_list]
        for f in concurrent.futures.as_completed(results):
            print(str(f.result()[0])+ "\t" + str(f.result()[1])+"\t" + str(f.result()[2])+"\t" + str(f.result()[3]))
            #apples=f.result()[2]/5
            #print(str(apples))
            ANI_no_gaps.append(f.result()[2])
            ANI_gaps.append(f.result()[3])
    finish = time.perf_counter()
    print("### average ani with gaps is: " + str(sum(ANI_gaps) / len(ANI_gaps)))
    print("### minimum ani with gaps is: " + str(min(ANI_gaps)))
    print("### average ani without gaps is: " + str(sum(ANI_no_gaps) / len(ANI_no_gaps)))
    print("### minimum ani without gaps is: " + str(min(ANI_no_gaps)))

    #print ("True ANI calc took " + str(finish-start) + " seconds")


