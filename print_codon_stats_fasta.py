# 
from __future__ import print_function
import sys
from Bio import SeqIO
from collections import Counter

fastafile = sys.argv[1]


def allCodons():
    alphabet = ('a','c','g','t')
    for c1 in alphabet:
        for c2 in alphabet:
            for c3 in alphabet:
                yield( ''.join((c1, c2, c3)) )

"""
Generate all codons in a CDS sequence
"""
def codons(seq):
    assert(len(seq) % 3 == 0)
    for i in [3*a for a in range(len(seq)/3)]:
        codon = str.lower(seq[i:i+3])
        assert(len(codon)==3)
        yield(codon)

class LastStreamElement(object):
    def process(self, item):
        pass

    def setNext(self, nextobj):
        assert(False)

class StreamElement(object):
    def __init__(self):
        self._next = LastStreamElement()
    
    def setNext(self, nextobj):
        self._next = nextobj
    
    def handle(self, item):
        self.process(item)
        self._next.process(item)

class CodonsFastaSource(StreamElement):
    def run(filename):
        for seq in SeqIO.parse(open(filename, 'r'), 'fasta'):
            for codon in codons(str(seq.seq)):
                self.handle(codon)


def fastaSource(filename):
    for seq in SeqIO.parse(open(filename, 'r'), 'fasta'):
        for codon in codons(str(seq.seq)):
            yield codon
        


counts = Counter(dict([(a, 0) for a in allCodons()]))
codonCount = 0
for codon in fastaSource(fastafile):
    counts.update((codon,))
    codonCount += 1

#assert(len(counts)==64) # Not true because of 'n's

#print(codonCount)
#print(counts.items())
print(counts.most_common(10))
        
    
