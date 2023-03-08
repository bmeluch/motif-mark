#!/usr/bin/env python

import cairo
import argparse

# GitHub: https://github.com/bmeluch/motif-mark

# Goals: Python script using object-oriented code to visualize motifs on sequences

class Motif:
    def __init__(self):
        '''asdf'''
        pass

    def create_regex(self):
        '''asdf'''
        pass

    def find_motif(self, gene: Gene):
        '''asdf'''
        pass

class Gene:
    def __init__(self, gene_name: str, gene_seq: str):
        '''asdf'''
        self.name = gene_name
        self.seq = gene_seq
        self.motifs = dict() # key: motif name, value: list of start positions
        pass

    def detect_regions(self, motif: Motif):
        '''asdf'''
        # something about regexes
        pass

class Figure():
    pass

class Gene_Vis():
    pass


def get_filepaths() -> tuple:
    parser = argparse.ArgumentParser(description="Get files containing genes, motifs")
    parser.add_argument("-f", help="FASTA file containing gene sequences", required=True)
    parser.add_argument("-m", help="Text file containing motifs, one per line", required=True)
    args = parser.parse_args()
    gene_filepath: str = args.f
    motif_filepath: str = args.m
    return gene_filepath, motif_filepath


def make_oneline_fasta(input_file: str, output_file: str):
    '''Given a FASTA file, creates a new FASTA file with newlines
    removed from sequences'''

    header: str = ""
    seq: str = ""

    with open(input_file, 'r') as in_fa, open(output_file, 'w') as out_fa:
        for line in in_fa:

            # Every time you hit a header line, remove newlines, write to file, get new record
            if line[0]==">":
                seq = seq.replace("\n", "")
                if header and seq:
                    out_fa.write(header+'\n')
                    out_fa.write(seq+'\n')
                header = line.strip()
                seq = ""

            else:
                seq += line
        # Finish writing out last record
        seq = seq.replace("\n", "")
        out_fa.write(header+'\n')
        out_fa.write(seq+'\n')
    return(output_file+" Sequences converted to single line")


gene_filepath, motif_filepath = get_filepaths()

# extract filenames from file paths
gene_fname: str = "test_filename"

# make gene objects from fasta file, set of gene objects

make_oneline_fasta(gene_filepath, gene_fname)

# make motif objects from motif file, set of motif objects

# search each motif object in each gene
# save hit start position


#     Capable of handling:
#         Motifs with ambiguous nucleotides (see https://en.wikipedia.org/wiki/Nucleic_acid_notation)
#         Handle overlapping motifs
#         Denote introns/exons





# Output
#   Same prefix as input file (e.g. Figure_1.fa -> Figure_1.png)
#   All features (motifs, introns, exons) should be to scale
#   One png per FASTA file
#   Well-labeled, with key