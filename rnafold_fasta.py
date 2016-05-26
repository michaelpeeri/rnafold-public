# Calculate the RNA folding energy for each protein contained in a fasta file, using the Vienna RNA package.
from Bio import SeqIO
import RNA
for seq in SeqIO.parse(open("./data/GCF_000002595.1_v3.0_rna.fna"), 'fasta'):
  struct, energy = RNA.fold(str(seq.seq))
  print("{0}\t{1}".format(seq.id, energy))
