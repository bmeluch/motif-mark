#!/usr/bin/env python

import cairo
import argparse
import re

# GitHub: https://github.com/bmeluch/motif-mark
# Goals: Python script using object-oriented code to visualize motifs on sequences

# IUPAC degenerate bases code
degenerate_bases: dict = {'A': 'A', 'T': 'TU', 'G': 'G', 'C': 'C', 'U': 'TU',
                          'R': 'AG', 'Y': 'CT', 'S': 'GC','W': 'AT',
                          'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
                          'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}

color_codes: list = [(0.6, 0, 0, 0.8), (0, 0.6, 0, 0.8), (0, 0, 0.6, 0.8), 
                     (0.6, 0.6, 0, 0.8), (0.6, 0, 0.6, 0.8)]

class Motif:
    def __init__(self, motif_seq: str, motif_color: tuple):
        '''Defines a Motif object containing a nucleotide sequence 
        and a designated color for figure rendering.'''
        self.seq = motif_seq.upper()
        self.color: tuple = motif_color
        pass

    def create_regex(self):
        '''Adds a "regex" attribute  to the Motif object. 
        The "regex" attribute is a compiled regular expression pattern 
        incorporating degenerate nucleotides according to IUPAC notation.'''
        pattern: str = ""
        for base in self.seq:
            pattern += ("["+degenerate_bases[base]+"]")
        self.regex = re.compile(pattern)
        pass


class Gene:
    def __init__(self, gene_name: str, gene_seq: str):
        '''Defines a Gene object with a name (eg. from FASTA header) and
        a sequence attribute. Initializes empty structures to hold identified
        motif and exon positions.'''
        self.name = gene_name
        self.seq = gene_seq
        self.motifs = dict() # key: motif name, value: list of start positions
        self.exons: list = []
        pass

    def detect_exons(self):
        '''Identifies start and end positions of exons in the gene sequence
        based on capital letters in the string. Saves exons as a list of tuples
        containing start and end positions in the gene.'''
        index: int = 0
        in_exon: bool = False
        exon_start: int = 0
        exon_end: int = 0
        for ch in self.seq:
            if ch.isupper() and not in_exon:
                exon_start = index
                in_exon = True
            elif ch.isupper() and in_exon:
                pass
            elif not ch.isupper() and in_exon:
                exon_end = index
                self.exons.append((exon_start, exon_end))
                in_exon = False
            elif not ch.isupper() and not in_exon:
                pass
            index += 1
        pass

    def find_motif(self, motif: Motif):
        '''Using the "regex" attribute of Motif objects, finds all start 
        positions in the gene sequence for a given motif.'''
        found_motifs = motif.regex.finditer(self.seq)
        start_pos: list = []

        for f_motif in found_motifs:
            start_pos.append(f_motif.start())
        if start_pos:
            self.motifs[motif.seq] = start_pos
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
    removed from sequences.'''

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
    pass


def main():
    gene_filepath, motif_filepath = get_filepaths()

    # Extract filename from file path
    gene_fname = re.match(r'(.+)\..+$', gene_filepath).group(1)

    # Create a set of Gene objects from a FASTA file
    gene_set = set()

    oneline_fasta_filepath: str = gene_fname + "_oneline.fasta"
    make_oneline_fasta(gene_filepath, oneline_fasta_filepath)

    with open(oneline_fasta_filepath, 'r') as gene_file:
        for line in gene_file:
            if line[0]==">":
                header: list = line[1:].strip().split()
            
            else:
                seq: str = line.strip()
                gene_set.add(Gene(header[0], seq))

    # Detect exons, then make gene sequence capitalized for easy searching
    # Also detect longest gene length for image
    longest_gene: int = 0
    for gene in gene_set:
        gene.detect_exons()
        gene.seq = gene.seq.upper()
        if len(gene.seq) > longest_gene:
            longest_gene = len(gene.seq)

    # Create a set of Motif objects from a text file
    motif_set = set()
    motif_count: int = 0
    with open(motif_filepath, 'r') as motif_file:
        for line in motif_file:
            motif_set.add(Motif(line.strip(), color_codes[motif_count]))
            motif_count += 1

    # Search each Motif object in each Gene object
    for motif in motif_set:
        motif.create_regex()
        for gene in gene_set:
            gene.find_motif(motif)

    # Create figure displaying motif positions on genes
    HEIGHT_PER_GENE: int = 125
    LEFT_MARGIN: int = 50
    HEIGHT: int = (len(gene_set)+1)*HEIGHT_PER_GENE
    WIDTH: int = longest_gene + 100

    with cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT) as surface:
        
        ctx = cairo.Context(surface)
        save_point = tuple()

        ctx.set_source_rgb(1,1,1)
        ctx.rectangle(0, 0, WIDTH, HEIGHT)
        ctx.fill() # Fill background with opaque white

        gene_count: int = 0
        for gene in gene_set:
            gene_count += 1
            # Draw genes as lines
            ctx.move_to(LEFT_MARGIN, (gene_count*HEIGHT_PER_GENE)-(HEIGHT_PER_GENE/2))
            ctx.rel_line_to(len(gene.seq), 0)
            ctx.set_source_rgb(0, 0, 0)
            ctx.set_line_width(3)
            ctx.stroke()
            
            # Add gene names under gene lines
            ctx.move_to(LEFT_MARGIN, (gene_count*HEIGHT_PER_GENE)-(HEIGHT_PER_GENE/4))
            ctx.set_source_rgb(0, 0, 0)
            ctx.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            ctx.set_font_size(14)
            ctx.show_text(gene.name)

            # Draw exons as transparent gray lines over the gene
            if gene.exons:
                for exon in range(0,len(gene.exons)):
                    ctx.set_source_rgba(0.1, 0.1, 0.1, 0.7) # transparent dark gray
                    exon_start: int = gene.exons[exon][0]
                    exon_length: int = (gene.exons[exon][1] - exon_start)
                    ctx.move_to(LEFT_MARGIN + exon_start, (gene_count * HEIGHT_PER_GENE) - (HEIGHT_PER_GENE / 2))
                    ctx.rel_line_to(exon_length, 0)
                    ctx.set_line_width(20)
                    ctx.stroke()

            # Draw motifs as semitransparent blocks using Motif object color
            for motif in motif_set:
                ctx.set_source_rgba(*motif.color)
                if motif.seq in gene.motifs:
                    for pos in gene.motifs[motif.seq]:
                        ctx.move_to(LEFT_MARGIN + pos, (gene_count*HEIGHT_PER_GENE)- (HEIGHT_PER_GENE / 2))
                        ctx.rel_line_to(len(motif.seq),0)
                        ctx.set_line_width(30)
                        ctx.stroke()

        # Display motif colors in a key at the bottom
        ctx.move_to(LEFT_MARGIN, (gene_count + 1) * HEIGHT_PER_GENE - (HEIGHT_PER_GENE / 2))
        ctx.set_source_rgb(0, 0, 0)
        ctx.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(14)
        ctx.show_text("KEY")

        ctx.move_to(LEFT_MARGIN, (gene_count + 1) * HEIGHT_PER_GENE - (HEIGHT_PER_GENE * 0.375))

        for motif in motif_set:
            ctx.set_source_rgb(0, 0, 0)
            ctx.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            ctx.set_font_size(14)
            ctx.show_text(motif.seq)

            ctx.rel_move_to(20, 0)
            save_point = ctx.get_current_point()
            ctx.set_source_rgba(*motif.color)
            ctx.rectangle(*ctx.get_current_point(), 10, -10)
            ctx.fill()
            ctx.move_to(*save_point)
            ctx.rel_move_to(100,0)

        # Add exon color to key
        ctx.move_to(LEFT_MARGIN, (gene_count + 1) * HEIGHT_PER_GENE - (HEIGHT_PER_GENE / 4))
        ctx.set_source_rgb(0, 0, 0)
        ctx.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(14)
        ctx.show_text("EXON")

        ctx.rel_move_to(20,0)
        ctx.set_source_rgba(0.1, 0.1, 0.1, 0.7)
        ctx.rectangle(*ctx.get_current_point(), 10, -10)
        ctx.fill()

        surface.write_to_png (gene_fname+".png")

if __name__ == "__main__":
    main()