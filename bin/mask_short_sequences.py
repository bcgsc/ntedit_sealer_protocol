#!/usr/bin/env python3

import argparse
import re


'''
Convert ACTG sequences less than k to lower case or "N"
'''

def extend_gaps(reads_filename, args):
    "Extend the gap lengths of the reads if flanks are less than k"
    with open(reads_filename, 'r') as fin:
        for line in fin:
            line = line.strip()
            if line[0] == ">":
                print(line)
            else:
                if len(line) < 2*args.k:
                    line = line.upper()
                else:
                    line = line[:args.k].upper() + line[args.k:-1*args.k] + line[args.k*-1:].upper()
                seq = ""
                groups_seq = re.findall(r"([ACTG]+|[Nn]+|[actgUNMRWSYKVHDBunmrwsykvhdb]+)", line)
                for i, new_seq in enumerate(groups_seq):
                    if new_seq[0] == "N":
                        seq += new_seq
                    else:
                        if len(new_seq) < args.k:
                            if args.n:
                                seq += "N"*len(new_seq)
                            else:
                                seq += new_seq.lower()
                        else:
                            seq += new_seq
                seq = seq.strip("Nn")
                if seq == "":
                    seq = "N"
                print(seq)

def main():
    parser = argparse.ArgumentParser(description="Hard or soft-mask ACTG sequences less than supplied k value")
    parser.add_argument("-n", help="Hard-mask regions less than k", action="store_true")
    parser.add_argument("-s", help="Soft-mask regions less than k", action="store_true")
    parser.add_argument("-k", help="K-mer size", required=True, type=int)
    parser.add_argument("FASTA", help="Input fasta file for masking (or - to read from stdin)")
    args = parser.parse_args()

    if not args.n and not args.s:
        parser.error("Either -h or -s must be set")
    if args.n and args.s:
        parser.error("Both -h and -s cannot be set - choose one to hard mask OR soft mask")

    reads = "/dev/stdin" if args.FASTA == "-" else args.FASTA

    extend_gaps(reads, args)

if __name__ == "__main__":
    main()