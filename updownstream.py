#!/usr/bin/env python
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

__author__ = "Gene Blanchard"
__email__ = "me@geneblanchard.com"

def process_fasta(fasta):
    with open(fasta, 'r') as fasta_h:
        contigs = []
        for record in SeqIO.parse(fasta_h, "fasta"):
            contigs.append(record)
    return contigs


def main():

    START_CODONS = ["atg", "gtg", "ttg", "att"]

    # Argument Parser
    parser = argparse.ArgumentParser(description='Grab a specified number of bases up and downstream from genome')
    # Input file
    parser.add_argument('-i', '--input', dest='input', required=True, help='The genomic features files')
    # Input file
    parser.add_argument('-f', '--fasta', dest='fasta', required=True, help='The fasta file to search')
    # Upstream
    parser.add_argument('-u', '--up', dest='upstream', type=int, required=True, help='The number of bases to grab upstream')
    # Downstream
    parser.add_argument('-d', '--down', dest='downstream', type=int, required=True, help='The number of bases to grab downstream')
    # Output
    parser.add_argument('-o', '--outfile', dest='outfile', required=True, help='The basename for the output files')

    # Parse arguments
    args = parser.parse_args()
    infile = args.input
    fasta = args.fasta
    upstream = args.upstream
    downstream = args.downstream
    outfile = args.outfile

    # Process the fasta
    contigs = process_fasta(fasta)
    # Output files
    posuphit = "{}_pos_up_hit.fa".format(outfile)
    posuphit_h = open(posuphit, 'w')
    print "open"
    posupmiss = "{}_pos_up_miss.fa".format(outfile)
    posupmiss_h = open(posupmiss, 'w')
    print "open"
    posdownhit = "{}_pos_down_hit.fa".format(outfile)
    posdownhit_h = open(posdownhit, 'w')
    print "open"
    posdownmiss = "{}_pos_down_miss.fa".format(outfile)
    posdownmiss_h = open(posdownmiss, 'w')
    print "open"
    neguphit = "{}_neg_up_hit.fa".format(outfile)
    neguphit_h = open(neguphit, 'w')
    print "open"
    negupmiss = "{}_neg_up_miss.fa".format(outfile)
    negupmiss_h = open(negupmiss, 'w')
    print "open"
    negdownhit = "{}_neg_down_hit.fa".format(outfile)
    negdownhit_h = open(negdownhit, 'w')
    print "open"
    negdownmiss = "{}_neg_down_miss.fa".format(outfile)
    negdownmiss_h = open(negdownmiss, 'w')
    print "open"
    # Process genomic record
    with open(infile, 'r') as infile_h:
        for line in infile_h:
            if not line.startswith("genome_id"):
                data = line.rstrip('\n').split('\t')
                feature_type = data[5]
                patric_id = data[6]
                start = int(data[9])
                end = int(data[10])
                strand = data[11]
                product = data[14]
                # Forward strand
                if strand == '+':
                    for sequence in [str(contig.seq) for contig in contigs]:
                        start = start - 1
                        # Find start codon
                        if start - upstream < 0:
                            upstart = 0
                        else:
                            upstart = start - upstream
                        if sequence[start:start + 3] in START_CODONS:
                            # posuphit_h
                            upstream_seq = sequence[upstart:start + 3]
                            output = ">{} | {} | {}\n{}\n".format(patric_id, feature_type, product, upstream_seq)
                            posuphit_h.write(output)
                            # posdownhit_h
                            downstream_seq = sequence[start:start + downstream + 1]
                            output = ">{} | {} | {}\n{}\n".format(patric_id, feature_type, product, downstream_seq)
                            posdownhit_h.write(output)
                        else:
                            # posupmiss_h
                            upstream_seq = sequence[upstart:start + 3]
                            output = ">{} | {} | {}\n{}\n".format(patric_id, feature_type, product, upstream_seq)
                            posupmiss_h.write(output)
                            # posdownmiss_h
                            downstream_seq = sequence[start:start + downstream + 1]
                            output = ">{} | {} | {}\n{}\n".format(patric_id, feature_type, product, downstream_seq)
                            posdownmiss_h.write(output)
                elif strand == '-':
                    for sequence in [str(contig.seq.reverse_complement()) for contig in contigs]:
                        # Find start codon
                        rstart = len(sequence) - end
                        if rstart - upstream < 0:
                            upstart = 0
                        else:
                            upstart = rstart - upstream
                        if sequence[rstart:rstart + 3] in START_CODONS:
                            # neguphit_h
                            upstream_seq = sequence[upstart:rstart + 3]
                            output = ">{} | {} | {}\n{}\n".format(patric_id, feature_type, product, upstream_seq)
                            neguphit_h.write(output)
                            # negdownhit_h
                            downstream_seq = sequence[rstart:rstart + downstream + 1]
                            output = ">{} | {} | {}\n{}\n".format(patric_id, feature_type, product, downstream_seq)
                            negdownhit_h.write(output)
                        else:
                            # negupmiss_h
                            upstream_seq = sequence[upstart:rstart + 3]
                            output = ">{} | {} | {}\n{}\n".format(patric_id, feature_type, product, upstream_seq)
                            negupmiss_h.write(output)
                            # negdownmiss_h
                            downstream_seq = sequence[rstart:rstart + downstream + 1]
                            output = ">{} | {} | {}\n{}\n".format(patric_id, feature_type, product, downstream_seq)
                            negdownmiss_h.write(output)
    posuphit_h.close()
    posupmiss_h.close()
    posdownhit_h.close()
    posdownmiss_h.close()
    neguphit_h.close()
    negupmiss_h.close()
    negdownhit_h.close()
    negdownmiss_h.close()



if __name__ == '__main__':
    main()
