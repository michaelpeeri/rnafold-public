from itertools import compress
from copy import copy
from math import factorial
import random
from Bio.Seq import Seq
import Bio.Data

def splitCodons(seq):
    assert(len(seq)%3==0)
    return [seq[i*3:i*3+3] for i in range(len(seq)/3)]

"""
Convert a forward codon table to a backwards codon table (AA->codons)
Note: make_back_table from Bio.Data.CodonTable only returns the first codon
Note: the returned table does not include include entries for start/stop codons
"""
def make_full_back_table(forward_table):
    ret = {}
    for codon,aa in forward_table.items():
        codon = codon.lower()
        if aa is None:
            continue
        
        if aa in ret:
            ret[aa].append(codon)
        else:
            ret[aa] = [codon]
    return ret

def applyCodonsPermutation(mutableSeq, newCodons, positions):
    for newCodon, pos in zip(newCodons, positions):
        mutableSeq[pos] = newCodon

class SynonymousCodonPermutingRandomization(object):
    """
    skipIdenticalSequences - refuse returning sequences identical to the source (even though they are valid permutation). This means the sampling from the population of random perumtations is done without replacement; This can be thought of as either distorting the distribution or avoiding distrortion; The difference is only noticeble for sequences which have a very small number of permutations (compared to the number actually generated).
    """
    def __init__(self, geneticCode=1, skipIdenticalSequences=True):
        self._geneticCode = geneticCode
        self._code = Bio.Data.CodonTable.unambiguous_dna_by_id[geneticCode]
        self._back_table = make_full_back_table(self._code.forward_table)
        self._skipIdenticalSequences = skipIdenticalSequences

    def randomize(self, nucleotideSeq):
        seq = Seq(nucleotideSeq.lower())
        codonsSeq = splitCodons(nucleotideSeq.lower())
        origSeqTranslation = seq.translate(table=self._code)

        randomizedCodonsSeq = copy(codonsSeq)
        permutationsCount = 1
        for aa,aaCodons in self._back_table.items(): # Iterate over each AA and its codons
            aaCodonsSet = frozenset(aaCodons)
            isSynonymousPosition = map( lambda x: x in aaCodonsSet, codonsSeq)
            synonymousPositions = filter(lambda x:x != -1, [i if v else -1 for v,i in zip(isSynonymousPosition, range(len(isSynonymousPosition)))])

            if(sum([int(x) for x in isSynonymousPosition]) < 2): # nothing to permute
                continue

            pool = list(compress(codonsSeq, isSynonymousPosition))

            # Calculate the total number of permutations, which is the product
            # of the numbers of permutations for each AA.
            # The permutations for each AA are a multiset permutation:
            # Source: https://en.wikipedia.org/wiki/Permutation#Permutations_of_multisets
            permutationsCountDenom = 1
            for c in aaCodons:
                countOfC = sum( [1 if x==c else 0 for x in codonsSeq] )
                permutationsCountDenom *= factorial(countOfC)
            permutationsCount *= factorial(len(pool)) / permutationsCountDenom

            # Permute synonymous codons for the current aa
            random.shuffle(pool)
            applyCodonsPermutation( randomizedCodonsSeq, pool, synonymousPositions)

        # Final checks
        # refuse providing permutations if the number of possible permutations is very small
        if(permutationsCount < 50):
            print("Warning: sequence only has %d possible permutations" % permutationsCount)
            raise Exception("Sequence only has %d possible permutations" % permutationsCount)

        # calculate nucleotide identity to the original sequence
        identityCount = sum([x==y for x,y in zip(''.join(randomizedCodonsSeq), seq)])
        assert(identityCount >= 0 and identityCount <= len(codonsSeq)*3)
        identity = float(identityCount)/len(seq)

        # refuse returning a sequence identical to the source (even though it is a valid permutation)
        if( self._skipIdenticalSequences and identity > 0.9999 ):
            return self.randomize(nucleotideSeq)

        # warn if the permutation only effects a small number of positions
        # note: identities are normally over 67%
        if( identity > 0.9 ):
            print("Warning: shuffled sequence has %.1f%% identity to original sequence" % (identity*100))


        resultingSeq = ''.join(randomizedCodonsSeq)
        assert(Seq(resultingSeq).translate(table=self._code) == origSeqTranslation)  # translation was not maintained by randomization!
        return (permutationsCount, identity, resultingSeq)


def test():
    c = SynonymousCodonPermutingRandomization()
    s = set()
    totalPermutationsCount = None
    testSeq = 'atcccgcgcccacctaataacacagcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcgatcgcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcgatcgcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcgatcgcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcg'
    #testSeq =  'atcgcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcgatcgcgatcagcgaaccaata'
    #testSeq =  'atcgcggcagcagcg'
    for i in range(1000):
        totalPermutationsCount, identity, perm = c.randomize(testSeq)
        s.add( perm )
    print('-----------------')
    print('Found %d randomizations' % len(s))
    print('Sequence length: %d codons' % (len(testSeq)/3))
    print('Total randomizations: %g' % float(totalPermutationsCount))
    return 0

if __name__=="__main__":
    import sys
    sys.exit(test())
