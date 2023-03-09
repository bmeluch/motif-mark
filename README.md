# motif-mark
Motif Mark assignment for BI625

```motif-mark-oop.py``` is a Python script which uses object-oriented code to produce a visualization of nucleotide motifs on gene sequences.

## Run Requirements
The script must be run in a conda environment with ```pycairo``` installed.

## Input Requirements
* ```-f```: A FASTA file containing gene sequences
   * Up to 10 gene sequences
   * Each sequence up to 1000 nucleotides in length
   * Exons denoted as CAPITAL nucleotides, otherwise lowercase nucleotides
* ```-m```: A text file containing motif sequences
   * Up to 5 motifs
   * One motif sequence per line
   * Motifs should use IUPAC nucleotide notation

## Output
* The script will output a png figure displaying genes, exons, and motifs to scale.
  * Genes are represented by a black line in the final image.
  * Exons are represented by a semitransparent gray rectangle at the corresponding position on the gene line.
  * Motifs are indicated by colored semitransparent markings along the black gene line.

## Considerations
* The script will identify motifs with degenerate nucleotides that are coded according to IUPAC notation.
* The script will identify overlapping motifs; the overlap will be visible with semitransparent colors on the final figure.
