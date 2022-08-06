'''
Caitlin Roake
Given a directory full of fastq RACE files, goes through each file
and tests each sequence to see if it contains a RACE product (defined
as having the ligated adapter, and the 3' end of dmenU2). If so,
takes the 'tail' between the 3' end and the adapter and outputs
it into a text file with each sequence on a new line. This program in
particular is designed to find sequences that are shorter than the full
length dmenU2
'''

import os
import re
from Bio.Seq import Seq

root_dir = input("Enter a file location")
counter = 0
adapter_count = 0
forward_hit = 0
reverse_hit = 0
adapter = Seq("CTGTAGGCACCATCAATCGTTACGTAG")  # sequence of the adapter
gene_3prime = Seq("CTGTCACGGGTTGGCCCGGT")  # 3'U2 distinguishing
gene_3prime_end = Seq("CCGCCGGGATTTCGGCCCAAC") # far end of U2 3'


def RACE(sequence, adapter, gene_3prime):
    '''
    takes a sequence from a fastq file, an adapter sequence, and the 3' gene
    sequence, and returns the 'tail' sequence between the adapter and
    the  3' end. If there's missing sequence between those two,
    returns the string 'short'
    '''
    first_cut = sequence.split(adapter)
    first_cut_left = first_cut[0]
    second_cut = first_cut_left.split(gene_3prime)
    if len(second_cut) > 1:
        return (second_cut[1])
    else:
        return ("short")


for file_name in sorted(os.listdir(root_dir)):
    handle = re.sub(r"fastq", "", file_name)
    output_file = open(str(handle) + "RACEoutput.txt", "w")
    file_name = os.path.join(root_dir, file_name)
    with open(file_name, "r") as fastq:
        counter = 0
        for line in fastq:
            counter += 1
            if counter == 2:
                sequence = Seq(line.rstrip("\r\n"))
                revcomp_sequence = sequence.reverse_complement()
                if adapter in sequence:
                    adapter_count += 1
                    if gene_3prime in sequence:
                        forward_hit += 1
                        output = RACE(sequence, adapter, gene_3prime)
                        output_file.write(str(output) + '\n')
                elif adapter in revcomp_sequence:
                    adapter_count += 1
                    if gene_3prime in revcomp_sequence:
                        reverse_hit += 1
                        output = RACE(revcomp_sequence, adapter, gene_3prime)
                        output_file.write(str(output) + '\n')
                else:
                    continue
            elif counter == 4:
                counter = 0
            else:
                continue
    print (str(adapter_count) + '\t' + str(forward_hit) + '\t' + str(reverse_hit))
    output_file.close()
