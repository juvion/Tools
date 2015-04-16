#!/usr/bin/python
__author__ = 'ju'
'''
task: compute conserved codon pair counts and inputed codon pair conserved counts across the alignment.
input: ICP14.in
output: gene_conservation.out
usage: gene_conservation.py ICP14.in gene_conservation.out
'''




import sys
import glob
from os.path import basename
from Bio import AlignIO
import os

WORKDIR = os.getcwd().split('analysis')[0]
# from Bio.Seq import Seq
from Bio.Alphabet import Gapped, IUPAC

#read in codon pair list (inhibitory codon pairs)
with open(sys.argv[1]) as fi:
    codon_pairs_list = [line.strip().upper().replace('U', 'T').replace('-','')[0:6] for line in fi]

#open write out file
fo = open(sys.argv[2], 'w')

#write the headers
header = 'gene_id\tgene_codons_num\tconserved_cp_count\tconservation_rate\tconserved_icp_count\tconserved_icp_report\n'
fo.write(header)

#read in alignment files
alignment_files = glob.glob(WORKDIR + 'ref/JohnstonDataset/coding/*codon.clustalw')
for alignment_file in alignment_files:
    # get the parse of alignment
    alignment = AlignIO.read(alignment_file, "clustal", alphabet=Gapped(IUPAC.unambiguous_dna))
    # remove the first start codon and last stop codon for shift0.
    # start_codon [AAG codon2 ... codonN] stop_codon
    alignment_shift0 = alignment[:, 3:-3]
    #start_codon [AG codon2 ... codonN+A] stop_codon
    alignment_shift1 = alignment[:, 4:-3] + alignment[:, 3:4]
    #start_codon [G codon2 ... codonN+AA] stop_codon
    alignment_shift2 = alignment[:, 5:-3] + alignment[:, 3:5]

    #get the name of the align ID (gene name) from the alignment file
    alignment_file_name = basename(alignment_file)
    align_id = alignment_file_name.split('_')[1]
    align_len = alignment.get_alignment_length()

    conserved_cp_count = 0
    conserved_icp_count = 0
    conserved_icp_report = []
    #iterate along the sequence for all codon pairs.
    for i in range(0, align_len - 3, 3):
        #create a list of algnment species number of codon pairs at each codon pair position. A column of codon pairs.
        align_codon_pairs = [alignment[j].seq[i:i + 6].tostring() for j in range(len(alignment))]
        if len(set(align_codon_pairs)) == 1:
            conserved_cp_count += 1
            if align_codon_pairs[0] in codon_pairs_list:
                conserved_icp_count += 1
                conserved_icp_report.append('.'.join([str(i/3), align_codon_pairs[0]]))
    conserved_icp_report = ', '.join(conserved_icp_report)

    #get gene_id from align_id, gene_codons_num from align_len (need to subtract initial Met and
    # stop codon) first.
    (gene_id, gene_codons_num) = (align_id, (align_len - 6) / 3)
    conservation_rate = conserved_cp_count / float(gene_codons_num)
    # write gene name, gene codons number, conserved codon pair count, conservation rate
    fo.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(gene_id, gene_codons_num, conserved_cp_count,
                                                round(conservation_rate, 4), conserved_icp_count,
                                                conserved_icp_report))
fo.close()
