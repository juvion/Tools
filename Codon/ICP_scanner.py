#!/usr/bin/python

'''scan each gene's alignment,
take a box size as input,
if any species has ICP been detected in the box
write down the box location and code of the conservation.
10010, first and fourth species have ICP found.
'''


# from collections import defaultdict
import sys
import glob
import collections
from os.path import basename
from Bio import AlignIO

# from Bio.Seq import Seq
from Bio.Alphabet import Gapped, IUPAC

#read in the list of codon pairs.
with open(sys.argv[1]) as fi:
    codon_pairs_list = [line.strip().upper().replace('U', 'T').replace('-','')[0:6] for line in fi]

#define a conservation relaxing window size, in + windows size frame, the same/is ICP will be considered as conserved.
window_size = 10 * 3 #set 10 codons as the window size
#read in the alignement.
alignment_files = glob.glob('/Users/ju/Works/Codon/ref/JohnstonDataset/coding/*codon.clustalw')
# alignment_files = glob.glob('*codon.clustalw')

#open codon pair list file, and put the codon pairs into a list, convert to upper case and replace 'U' as 'T'.
alignment_counter = 0
for alignment_file in alignment_files:
    #get the parse of alignment
    alignment = AlignIO.read(alignment_file, "clustal", alphabet=Gapped(IUPAC.unambiguous_dna))
    #get the name of the align ID (gene name) from the alignment file
    alignment_file_name = basename(alignment_file)
    align_id = alignment_file_name.split('_')[1]
    align_len = alignment.get_alignment_length()


    #create two dictionaries to store icp ID and counts for each species
    map_icp_ids = {}
    map_icp_counts = {}
    #iterate along species
    for j in range(len(alignment)):
        #get the species name
        species = alignment[j].id.split('_')[0]
        #convert gene sequence into a list of codon pairs along the sequence.
        gene_codon_pairs = [alignment[j].seq[i:i+6].tostring() for i in range(0, align_len - 3, 3)]
        
        #process the identity and count of the ICPs.
        id_report = []
        codon_index = 0
        icp_count = 0
        for codon_pair in gene_codon_pairs:
            if codon_pair in codon_pairs_list:
                id_report.append(str(codon_index) + '.' + codon_pair)
                icp_count += 1
            codon_index += 1

        #put the count of found ICP into the dictionary.
        map_icp_counts.update({species:str(icp_count)})
        #put the identity and index of found ICP into the dictionary.
        map_icp_ids.update({species:','.join(id_report)})
    
    #generate output.
    count_report = []
    info_report = []
    for species in ['Scer', 'Sbay', 'Smik', 'Skud', 'Spar']:
        if species in map_icp_counts.keys():
            count_report.append(map_icp_counts[species])
            if len(map_icp_ids[species]) == 0:
                info_report.append('NA')
            else:
                info_report.append(map_icp_ids[species])
        else:
            count_report.append('0')
            info_report.append('NA')
    print align_id + '\t' + '\t'.join(count_report) + '\t' + '\t'.join(info_report)
