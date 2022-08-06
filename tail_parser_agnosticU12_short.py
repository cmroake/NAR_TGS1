'''
Caitlin Roake
This script takes text output from hte RACE_parser_12 ***this only works for human U12 RNAs
and categorizes the tails into six bins-Mature (annotated end), short (<mature), Ext (>mature), and also counts number of extended and oligo A or T residues
It outputs a tab-delim file with the sequences and their bins
(sequences.txt)and outputs a summary csv file (RACE_summary.csv)
with some summary statistics for each file read. IT also outputs a
RACE_species.txt file  with the unique species found in that sample.
IT places all these files into a new directory called output, with
a timestamp.
'''

import os
import sys
import csv
import re
import datetime
import numpy as np
	

now = datetime.datetime.now()  # makes unique output folder name
genomic_sequence = '''ACCCTATTCACGCCTAAAAAGTAGACTGACTGTGGGGTGGTCGTGTTTTTTGTTTCTTGTTGGTAGGTGGTGAATG
CGTTTTTTTCGTTGTTTTCTCCGTTACTCAGGCTGCCAGTTGCTTGGCAGTCTTGTCGCTGGCTGTGGACGCTCTGCACT'''

gene_3prime_end = "CGGGATGCCTGGGAGTTGCGATCTGCCCG" # far end of U4atac 3'.

root_dir = input("Input folder location")
output_dir = root_dir + "\\output" + str(now.minute) + str(now.second)
os.mkdir(output_dir)
output_dir = output_dir + "\\"
output_header = "RACE\tpolyA\tpolyU\tExt\tAbn\tA_len\tExt_len\tShort_len\tU_len\n"
output_summary = open(output_dir+"U2waterfall_RACE_summary.csv", "wt")
summary_writer = csv.writer(output_summary, delimiter=',',
                            lineterminator='\n')
summary_header = ["File", "TotalReads", "Mature", "M%", "Short", "S%",
                  "Ext", "Ext%", "PolyA", "PolyA%","PolyU","PolyU%", "AtailMean","Atailmed", "EtailMean",
                  "Etailmed","Utailmed","Utail_mature","Utail_extended"]
summary_writer.writerow(summary_header)

def check_short(short_seq):
    '''This method takes a sequences shorter than gene_3prime_end and makes
    sure it is u4atac and not just garbage. Returns 'short' if it's a true
    short sequence and 'error' if it's garbage.'''
    i = 0
    for letter in short_seq:
        
        if letter == gene_3prime_end[i]:
            i += 1
        else:
            return ('error')
    return ('short')

for file_name in sorted(os.listdir(root_dir)):
    handle = file_name
    if handle.endswith(".txt"):
        graphing_name = re.sub(r"\..*", "", handle,)
        file_name = os.path.join(root_dir, file_name)
        output_sequences = open(output_dir + str(graphing_name) +
                                "U4atacwaterfall_sequences.txt", "wt")
        output_sequences.write(output_header)
        RACE_species = open(output_dir + str(graphing_name) +
                           "U4atacwaterfall_RACE_species.txt", "w")
        short, mature_num, extended_num, polya_num, polyu_num = 0, 0, 0, 0, 0
        totalAs, totalReads, totalExt, Utail_mature, Utail_ext = 0, 0, 0, 0, 0
        A_median = []
        E_median = []
        U_median = []
        RACE_set = set()
        with open(file_name, "r") as RACE_file:
            for line in RACE_file:
                RACE = line.rstrip('\n')
                #print(str(RACE))
                if len(RACE) == len(gene_3prime_end) and RACE == gene_3prime_end:  # a nonextended species, 'mature'
                    #print (str(len(gene_3prime_end))
                    #print (str(len(RACE)))
                    mature_num += 1
                    totalReads += 1
                    sequences_row = "\t".join([RACE, '0', '0', '0', '0', '0','0','0','0']) 
                    sequences_row = sequences_row + "\n"
                    output_sequences.write(sequences_row)
                    continue
                if len(RACE) < len(gene_3prime_end):  # less than annoteated 3' end, 'short'
                    if check_short(RACE) == 'short':
                        short_len = len(gene_3prime_end)-len(RACE)
                        short += 1
                        totalReads +=1
                        sequences_row = "\t".join([RACE, '0', '0', '0', '0', '0','0',str(short_len), '0', '0']) 
                        sequences_row = sequences_row + "\n"
                        output_sequences.write(sequences_row)
                        continue
                    else:
                        continue
                i, poly_a, poly_u, extended, abnormal, A_len, ext_len, U_len = 0, 0, 0, 0, 0, 0, 0, 0
                for letter in RACE[len(gene_3prime_end):]:  
                    if letter == genomic_sequence[i] and poly_a == 0 and poly_u == 0:
                        extended = 1
                        ext_len += 1
                        i += 1
                        continue
                    elif letter == "A":  # not genomic it's oligoA
                        poly_a = 1
                        A_len += 1
                        i += 1
                        continue
                    elif letter == "T": # not genomic it's oligoU
                        poly_u = 1
                        U_len += 1
                        continue
                    else:  # not genomic and not an A, 'abnormal'
                        abnormal = 1
                        break
                if abnormal == 0:
                    
                    
                    RACE_set.update([RACE])
                    totalReads += 1
            
                    if extended == 1:
                        
                        totalExt = totalExt + ext_len
                        E_median.append(ext_len)
                        extended_num += 1
                    if poly_a == 1 and A_len > 1:
                        polya_num += 1
                        totalAs = totalAs + A_len
                        A_median.append(A_len)
                   
                    if poly_u == 1:
                        polyu_num += 1
                        U_median.append(U_len)

                    if poly_u == 1 and extended == 0:
                        Utail_mature += 1
                        
                    if poly_u == 1 and extended == 1:
                        Utail_ext += 1
                        
                        
                sequences_row = "\t".join([RACE,
                                           str(poly_a),
                                           str(poly_u),
                                           str(extended),
                                           str(abnormal),
                                           str(A_len),
                                           str(ext_len), "0", str(U_len)])
                sequences_row = sequences_row + "\n"
                
                output_sequences.write(sequences_row)
    else:
        continue

    if totalReads > 0:  # as long as the files had good sequence
        RACE_list = "\n".join(sorted(RACE_set))
        RACE_species.write(RACE_list)
        RACE_species.close()
        Atail_avg = 0
        if polya_num > 0 :
            Atail_avg = totalAs / polya_num
        Etail_avg = totalExt / extended_num
        A_median = np.median(A_median)
        E_median = np.median(E_median)
        U_median = np.median(U_median)
        summary_row = [str(graphing_name), str(totalReads),
                       str(mature_num), str(mature_num/totalReads),
                       str(short), str(short/totalReads),
                       str(extended_num), str(extended_num/totalReads), str(polya_num),
                       str (polya_num/totalReads), str(polyu_num), str(polyu_num/totalReads),
                       str(Atail_avg), str(A_median),
                       str(Etail_avg),
                       str(E_median), str(U_median), str(Utail_mature/totalReads), str(Utail_ext/totalReads)]
        summary_writer.writerow(summary_row)
output_sequences.close()
output_summary.close()




