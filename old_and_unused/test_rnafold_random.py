from __future__ import print_function
import random
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import RNA
from runningstats import RunningStats


# configuration
windowWidth = 40
sequenceLength = 120
totalSequncesToTry = 10000 + 5
strongFasta = []
weakFasta = []



def randseq(N=100):
    return ''.join([c for c in itertools.starmap(lambda:random.choice('acgt'), itertools.repeat((), N))])

def randSeq_StrongFolding(N=100):
    parts = []

    def newrand():
        return randseq(random.randint(10,20))
    
    while( sum(map(len, parts)) < N ):
        if( len(parts) == 0 ):
            # Add the first part
            parts.append(newrand())
        else:
            # Add another part
            newseq = str(Seq(parts[-1]).reverse_complement())
            parts.append(str(newseq))
            parts.append(newrand())
    
    out = ''.join(parts)
    return out[:N]
    


for i in range(50):
    a = randSeq_StrongFolding(120)
    strct, energy = RNA.fold(a)
    assert(len(a)==120)
    print(a, energy)

import sys
sys.exit()


gstats = RunningStats()
i=0
for i in range(totalSequncesToTry):
    seq = randseq(sequenceLength)

    stats = RunningStats()

    for start in range(len(seq)-windowWidth+1):
        fragment = seq[start:(start+windowWidth)]
        assert(len(fragment)==windowWidth)
        
        # Calculate the RNA folding energy. This is the computation-heavy part.
        strct, energy = RNA.fold(fragment)
        assert(energy <= 0.0)
        stats.push(energy)

    if( stats.mean() > -3.2 and stats.min() > (-6.32-0.50) ):
        weakFasta.append( SeqRecord(Seq(seq, alphabet=generic_dna), id=("rand%d" % i), description=("Random synthetic sequence, mean MFE=%.4g, MFA>=%.4g" % (stats.mean(), stats.min()))))
    elif( stats.mean() < -8.03 and stats.max() < (-6.32+0.50)):
        strongFasta.append( SeqRecord(Seq(seq, alphabet=generic_dna), id=("rand%d" % i), description=("Random synthetic sequence, mean MFE=%.4g, MFA<=%.4g" % (stats.mean(), stats.max()))))

    gstats.push(stats.mean())
    if( i%500==499 ):
        print("weak: %d strong: %d" % (len(weakFasta), len(strongFasta)))
        
print("Totals: weak folding: %d seqs (%.2g%%) strong folding: %d seqs (%.2g%%)" % (len(weakFasta), float(len(weakFasta))/i*100, len(strongFasta), float(len(strongFasta))/i*100))
SeqIO.write( weakFasta, "weak_mfe.fna", "fasta")
SeqIO.write( strongFasta, "strong_mfe.fna", "fasta")

print("Mean-MFE - statistics: Min: %.3g Mean: %.3g Std: %.3g Max: %.3g" % (gstats.min(), gstats.mean(), gstats.stdev(), gstats.max()))

