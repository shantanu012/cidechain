#!/usr/bin/python

from Bio.Align.Applications import ClustalOmegaCommandline
in_file = "unaligned.fasta"
out_file = "aligned.fasta"
out_fmt = "clustal"

clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, outfmt=out_fmt, verbose=True, auto=True)
clustalomega_cline()

# clustal omega shell command example :
# clustalo -i unaligned.fasta -o aligned.fasta --auto -v --outfmt clustal

