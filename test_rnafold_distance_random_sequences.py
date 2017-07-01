from __future__ import print_function
import random
import itertools
from collections import Counter
from data_helpers import getAllNativeCDSsForSpecies, decompressNucleicSequence
import _distributed
from rnafold_vienna import RNAfoldWithStructure
from timeit import timeit

scheduler = _distributed.open()


# create a random dna-like sequence (with uniform distribution)
def randseq(N=100):
    return ''.join([c for c in itertools.starmap(lambda:random.choice('acgt'), itertools.repeat((), N))])

def printHist(c):
    maxval = float(max(c.values()))
    for val, count in c.items():
        print("%02d %s" % (val,"="*(int(round(count/maxval*80)))))
        
#c = Counter()
#for i in range(10000):
#    (pairs, _) = RNAfoldWithStructure(randseq(40))
#    c.update([x[1]-x[0] for x in pairs])


#printHist(c)

def test(fraction, taxId, numFractions):
    assert(fraction>=0)
    assert(fraction < numFractions)
    seqs = getAllNativeCDSsForSpecies(taxId, fraction, numFractions)
    return len(seqs)

def testNativeDistances(fraction, taxId, numFractions):
    assert(fraction>=0)
    assert(fraction < numFractions)
    
    numSeqs = 0
    counts = Counter()

    allSeqs = getAllNativeCDSsForSpecies(taxId, fraction, numFractions)
    
    for (seqId, seq) in allSeqs.items():
        cdsSeq = decompressNucleicSequence(seq)
        del seq
        
        for start in range(0, len(cdsSeq)-41, 10):
            end = start+40
            window = cdsSeq[start:end]
            assert(len(window)==40)

            (pairs, _) = RNAfoldWithStructure(window)
            counts.update([x[1]-x[0] for x in pairs])
            
        numSeqs += 1
    return (numSeqs, counts)


def testWithNFractions(numFractions):
    A = scheduler.map(test, range(numFractions), taxId=436017, numFractions=numFractions)
    B = scheduler.submit(sum, A)

    numSequences = B.result()

    print("%d fractions -> %d sequences" % (numFractions, numSequences))

    return numSequences
#newIds = testval.result()


def testDistances(numFractions=500):
    A = scheduler.map(testNativeDistances, range(numFractions), taxId=436017, numFractions=numFractions)
    B = scheduler.gather(A)

    totalCounts = Counter()
    numSeqs = 0

    for result in B:
        totalCounts += result[1]
        numSeqs += result[0]
    
    printHist(totalCounts)
    print("Num of sequences: %d" % numSeqs)
    print("Pairs: %.4g (estimated %.3g per window)" % (sum(totalCounts), float(sum(totalCounts))/(numSeqs*20)))




#tt = (1, 5, 10, 20, 50, 100, 200)
#tt = (1, 5, 10, 20)

#for t in tt:
#    secs = timeit( lambda : testWithNFractions(t), number=1 )
#    print("%.3gs" % secs)


testDistances()



