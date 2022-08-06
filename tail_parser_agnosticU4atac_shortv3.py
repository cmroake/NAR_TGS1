'''
Caitlin Roake
This script takes text output from the RACE_parser_u4 ***this only works for human U4atac RNAs
and categorizes the tails into six bins-Mature (annotated end), short (<mature), Ext (>mature),
and also counts number of extended and oligo A or T residues
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

now = input("Input unique run ID") # makes unique output folder name
genomic_sequence = '''AAAACCTGTTTTCATAGACTTATCAGTTCAAACAGCAGTAATTCGTAAAT
AAACTAGTACTTTGTGGTTAAACCAGTAGAGGGTGCACAAGACGCGTGGTTTTAGTGTCG
CAAGTAAAGTTCTTTCAGTTTTTGCGGTGATTAAACGGGAAGGATTTCAACAGAACTTCA
CCCCTTTAACTTTACGCCGATCATCAACTGTTCTGTAACTTCCTTATGCTGAGTCAGCCC'''

gene_3prime_end = "TGGTGCAATTTTTGGAAAAATG" # far end of U4atac 3'. 

root_dir = input("Input folder location")
output_dir = root_dir + "\\output" + str(now) 
os.mkdir(output_dir)
output_dir = output_dir + "\\"
output_header = "RACE\tpolyA\tpolyU\tExt\tA_len\tExt_len\tShort_len\tU_len\n"
output_summary = open(output_dir+"U4atacwaterfall_RACE_summary.csv", "wt")
summary_writer = csv.writer(output_summary, delimiter=',',
                            lineterminator='\n')
summary_header = ["File", "CountedReads", "Mature", "M%", "Short", "S%",
                  "Ext", "Ext%", "PolyA", "PolyA%","PolyU","PolyU%", "AtailMean","Atailmed", "EtailMean",
                  "Etailmed","Utailmed","Utail_mature","Utail_extended","Uncounted"]
summary_writer.writerow(summary_header)



def check_short(short_seq):
    '''This method takes a sequences shorter than gene_3prime_end and makes
    sure it is U4 and not just garbage. Returns 'short' if it's a true
    short sequence and 'error' if it's garbage.'''
    i = 0
    for letter in short_seq:
        
        if letter == gene_3prime_end[i]:
            i += 1
        else:
            return ('error')
    return ('short')
    

def check_mature(mature_seq):
    '''this method takes a sequence that is the same length as gene_3prime_end and makes
    sure it is U4 and not just garbage. Returns 'mature' if it's a true mature sequence
    and 'error' if it's just garbage.'''
    if len(mature_seq) == len(gene_3prime_end) and mature_seq == gene_3prime_end:
        return ('mature')
    if len(mature_seq) == len(gene_3prime_end):
        return ('nonmatch')
    else:
        return('immature')
    

def check_extension(extended_seq):
    i = poly_a = poly_u = extended = A_len = ext_len = U_len = 0
    if extended_seq[:len(gene_3prime_end)] != gene_3prime_end:
        return ('abnormal')
    for letter in extended_seq[len(gene_3prime_end):]:  
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
        else:
            return ('abnormal')

    return ([poly_a, poly_u, extended, A_len, ext_len, U_len])

'''this section goes through the txt files from the RACE parser
and assigns each transcript to a bin'''

for file_name in sorted(os.listdir(root_dir)):
    handle = file_name
    count = 0
    if handle.endswith(".txt"):
        graphing_name = re.sub(r"\..*", "", handle,)
        file_name = os.path.join(root_dir, file_name)
        output_sequences = open(output_dir + str(graphing_name) +
                                "U4atacwaterfall_sequences.txt", "wt")
        output_sequences.write(output_header)
        RACE_species = open(output_dir + str(graphing_name) +
                           "U4atacwaterfall_RACE_species.txt", "w")
        short = mature_num = extended_num = extended_only = polya_num = polyu_num = 0
        totalAs = totalReads = totalExt = Utail_mature = Utail_ext = 0
        A_median = []
        E_median = []
        U_median = []
        RACE_set = set()
        uncounted = 0
        with open(file_name, "r") as RACE_file:
            for line in RACE_file:
                count += 1
                totalReads += 1
                RACE = line.rstrip('\n')
                if check_mature(RACE) == 'mature':  # a nonextended species, 'mature'
                    mature_num += 1
                    sequences_row = "\t".join([RACE, '0', '0', '0', '0', '0','0','0','0']) 
                    sequences_row = sequences_row + "\n"
                    output_sequences.write(sequences_row)
                    continue
                if check_mature(RACE) == 'nonmatch':
                    uncounted += 1
                    continue
                if len(RACE) < len(gene_3prime_end):  # less than annotated 3' end, 'short'
                    if check_short(RACE) == 'short':
                        short_len = len(gene_3prime_end)-len(RACE)
                        short += 1
                        sequences_row = "\t".join([RACE, '0', '0', '0', '0','0',str(short_len), '0', '0']) 
                        sequences_row = sequences_row + "\n"
                        output_sequences.write(sequences_row)
                        continue
                    else:
                        uncounted += 1
                        continue
                if len(RACE) > len(gene_3prime_end):
                    transcript_desc = check_extension(RACE) 
                    if transcript_desc =='abnormal':
                        uncounted += 1
                        continue
                    if transcript_desc[2] == 1 or transcript_desc[0] == 1 or transcript_desc[1] == 1:
                        extended_num += 1
                        RACE_set.update([RACE])
                    if transcript_desc[2] == 1:
                        totalExt = totalExt + transcript_desc[4]
                        E_median.append(transcript_desc[4])
                    if transcript_desc[0] == 1 and transcript_desc[3] > 1:
                        polya_num += 1
                        totalAs = totalAs + transcript_desc[3]
                        A_median.append(transcript_desc[3])
                    if transcript_desc[1] == 1:
                        polyu_num += 1
                        U_median.append(transcript_desc[5])
                        if transcript_desc[2] == 1:
                            Utail_ext += 1
                        if transcript_desc[2] == 0:
                            Utail_mature += 1
                
                sequences_row = "\t".join([RACE,
                                           str(transcript_desc[0]),
                                           str(transcript_desc[1]),
                                           str(transcript_desc[2]),
                                           str(transcript_desc[3]),
                                           str(transcript_desc[4]),
                                           "0",
                                           str(transcript_desc[5])])
                sequences_row = sequences_row + "\n"
                output_sequences.write(sequences_row)
                
    else:
        continue


    if totalReads > 0:  # as long as the files had good sequence
        counted_reads = totalReads - uncounted
        print (str(count)+str(file_name) + "uncounted " + str(uncounted) + "totalReads" + str(totalReads))
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
        summary_row = [str(graphing_name),
                       str(counted_reads),
                       str(mature_num),
                       str(mature_num/counted_reads*100),
                       str(short),
                       str(short/counted_reads*100),
                       str(extended_num),
                       str(extended_num/counted_reads*100),
                       str(polya_num),
                       str (polya_num/counted_reads*100),
                       str(polyu_num),
                       str(polyu_num/counted_reads*100),
                       str(Atail_avg),
                       str(A_median),
                       str(Etail_avg),
                       str(E_median),
                       str(U_median),
                       str(uncounted)]
        summary_writer.writerow(summary_row)
output_sequences.close()
output_summary.close()





