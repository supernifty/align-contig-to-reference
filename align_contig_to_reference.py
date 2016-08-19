#!/usr/bin/env python
'''
  aligns a contig to a reference
'''

import argparse
import os
import random
import sys

import Bio.SeqIO
import Bio.SeqRecord

def log(msg):
    sys.stderr.write('{0}\n'.format(msg))

def align_contig_to_reference(kmer, contig, reference):
    '''
        generate fastq from contig and align to reference
    '''
    log('generating fastq...')
    fastq_filename = 'tmp{}'.format(random.randint(1,1000000))
    with open(fastq_filename, 'w') as fastq_fh:
        counter = 0
        for record in Bio.SeqIO.parse(open(contig, "rU"), "fasta"):
            for i in range(0, len(record.seq) - kmer + 1):
                substring = record.seq[i:i+kmer]
                fastq_rec = Bio.SeqRecord.SeqRecord(substring, id="{}".format(counter), description="{}/{}".format(record.id, i))
                fastq_rec.letter_annotations["phred_quality"] = [40] * len(substring)
                Bio.SeqIO.write(fastq_rec, fastq_fh, "fastq")
                counter += 1

    log('indexing reference...')
    os.system('bwa index {}'.format(reference))
    
    log('aligning to reference...')
    os.system('bwa mem -M -t 8 -k 19 {} {}'.format(reference, fastq_filename))

def main():
    '''
        execute via command line
    '''
    parser = argparse.ArgumentParser(description='Align a contig to a reference')
    parser.add_argument('--kmer', type=int, default=100, help='length of reads to extract')
    parser.add_argument('--contig', required=True, help='contig to align from')
    parser.add_argument('--reference', required=True, help='reference to align to')
    args = parser.parse_args()
    align_contig_to_reference(args.kmer, args.contig, args.reference)

if __name__ == '__main__':
    main()
