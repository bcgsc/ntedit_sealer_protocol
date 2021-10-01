#!/usr/bin/env python3

import argparse
import re


'''
Convert ACTG sequences less than k to lower case
'''

def extend_gaps(reads_filename, args):
    "Extend the gap lengths of the reads if flanks are less than k"
    with open(reads_filename, 'r') as fin:
        for line in fin:
            line = line.strip()
            if line[0] == ">":
                print(line)
            else:
                seq = ""
                groups_seq = re.findall(r"([ACTG]+|[Nn]+|[actgUNMRWSYKVHDBunmrwsykvhdb]+)", line)
                for i, new_seq in enumerate(groups_seq):
                    if new_seq[0] == "N":
                        seq += new_seq
                    else:
                        if len(new_seq) < args.k:
                            seq += new_seq.lower()
                        else:
                            seq += new_seq
                seq = seq.strip("Nn")
                if seq == "":
                    seq = "N"
                print(seq)

def main():
    parser = argparse.ArgumentParser(description="Soft-mask ACTG sequences less than supplied k value")
    parser.add_argument("-k", help="K-mer size", required=True, type=int)
    parser.add_argument("FASTA", help="Input fasta file for masking (or - to read from stdin)")
    args = parser.parse_args()

    reads = "/dev/stdin" if args.FASTA == "-" else args.FASTA

    extend_gaps(reads, args)

if __name__ == "__main__":
    main()