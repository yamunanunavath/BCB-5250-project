#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def interleave_fastq(forward_file, reverse_file, output_file):
    with open(forward_file, 'r') as f1, open(reverse_file, 'r') as f2, open(output_file, 'w') as out:
        for record1, record2 in zip(SeqIO.parse(f1, "fastq"), SeqIO.parse(f2, "fastq")):
            SeqIO.write(record1, out, "fastq")
            SeqIO.write(record2, out, "fastq")

# Define the file names (paths where the files are uploaded)
forward_file = 'bacterium_R1.fastq'
reverse_file = 'bacterium_R2.fastq'
output_file = 'interleaved_output.fastq'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Interleave Paired-End FASTQ Files')
    parser.add_argument('forward_file', help='Forward reads FASTQ file')
    parser.add_argument('reverse_file', help='Reverse reads FASTQ file')
    parser.add_argument('output_file', help='Output interleaved FASTQ file')
    args = parser.parse_args()
    interleave_fastq(args.forward_file, args.reverse_file, args.output_file)
