from __future__ import print_function
from math import log
from time import time
import os
import random
import itertools
from collections import Counter
from timeit import timeit
import sys
from Bio.Data import CodonTable
import logging
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
import pandas as pd
from data_helpers import getAllNativeCDSsForSpecies, decompressNucleicSequence, getSpeciesTranslationTable, countSpeciesCDS, getSpeciesName, allSpeciesSource, nativeSequencesSource
import _distributed
import config
from rnafold_vienna import RNAfoldWithStructure
from codon_randomization import SynonymousCodonPermutingRandomization


taxId = int(sys.argv[1])

scheduler = _distributed.open()

_log2 = log(2.0)
def log2(x):
    return log(x)/_log2

# create a random dna-like sequence (with uniform distribution)
def randseq(N=100):
    return ''.join([c for c in itertools.starmap(lambda:random.choice('acgt'), itertools.repeat((), N))])

def multiRandSeq(Q=10, N=500):
    return dict(list(map(lambda x: (x, randseq(N)), range(Q))))


#print(multiRandSeq())

def printHist(c):
    maxval = float(max(c.values()))
    for val, count in c.items():
        print("%02d %s" % (val,"="*(int(round(count/maxval*80)))))
        
#c = Counter()
#for i in range(10000):
#    (pairs, _) = RNAfoldWithStructure(randseq(40))
#    c.update([x[1]-x[0] for x in pairs])


#printHist(c)

# def test(fraction, taxId, numFractions):
#     assert(fraction>=0)
#     assert(fraction < numFractions)
#     seqs = getAllNativeCDSsForSpecies(taxId, fraction, numFractions)
#     return len(seqs)


class PairedStructureConsensus(object):
    def __init__(self):
        self._consensus = []

    def update(self, pairs):
        if not pairs:
            return

        assert(pairs==sorted(pairs)) # pairs must be sorted
        
        if not self._consensus:
            self._consensus.extend(pairs)
            return
            
        if pairs[0][0] > self._consensus[-1][1]:
            self._consensus.extend(pairs)

    def getConsensus(self):
        return self._consensus

    
class ModulusCounter(object):
    def __init__(self, N):
        self._N = N

        counters = []
        for n in range(N):
            counters.append( Counter() )
        self._counters = counters

    def updateWithIndices(self, instances):
        for n in [x%self._N for x in instances]:
            self._counters[n][1] += 1

    def getResults(self):
        return self._counters

    def __iadd__(self, other):
        assert(other._N == self._N)
        for n in range(self._N):
            self._counters[n] += other._counters[n]
        return self
            

def sequencePairsSource(cdsSeq, windowWidth=40, windowStep=10):
    for start in range(0, len(cdsSeq)-windowWidth-1, 10):
        end = start+windowWidth
        window = cdsSeq[start:end]
        assert(len(window)==windowWidth)

        (pairs, _) = RNAfoldWithStructure(window, cdsOffset=start)
        yield pairs

def sequencePairsSource2(cdsSeq, windowWidth=40, windowStep=10):
    for start in range(0, len(cdsSeq)-windowWidth-1, 10):
        end = start+windowWidth
        window = cdsSeq[start:end]
        assert(len(window)==windowWidth)

        (pairs, energy) = RNAfoldWithStructure(window, cdsOffset=start)
        yield (pairs, energy)
        

def createBackTable(forwardTable):
    out = {}
    for k,v in forwardTable.items():
        if v in out:
            out[v].append(k)
        else:
            out[v] = [k]
    
    return out

def getCodons(seq):
    if(len(seq) % 3 != 0):
        print("Warning: sequence does not contain complete codons, skipping...")
        return
        
    for n in range(0, len(seq)-3, 3):  # Skip the stop codon
        codon = seq[n:(n+3)]
        if( 'N' in codon ):
            continue
        yield codon

class CodonEntropyCalculation(object):
    def __init__(self, geneticCode):
        self._table = CodonTable.unambiguous_dna_by_id[geneticCode]
        self._backTable = createBackTable(self._table.forward_table)
        #print(self._backTable)
        
    def calcCDSCodonEntropy(self, cdsNucleotideSeq):

        H = []
        
        for codon in getCodons(cdsNucleotideSeq.upper()):
            alternatives = self._backTable[self._table.forward_table[codon]]

            for pos in (0,1,2):

                Pc = 1./float(len(alternatives))

                sigmaI = 0.0
                
                for n in ('A','C','G','T'):
                    countn = sum([1 for x in alternatives if x[pos]==n])

                    Pn = countn * Pc

                    if( Pn > 1e-20):
                        sigmaI += Pn * log2(Pn)

                if sigmaI == 0.0:
                    H.append(0.0)  # prevent inserting 0s as -0.0
                else:
                    H.append(-sigmaI)
                        
            
            #print("%s %d %.3g %.3g %.3g" % (codon, len(alternatives), H[-3], H[-2], H[-1]))

        return H
            

        
def roundH(x):
    return int(round(x,1)*10)

class CalcStats(object):
    def __init__(self, taxId):
        self._numSeqs = 0
        self._counts = Counter()
        
        self._mod3 = ModulusCounter(3)
        self._mod4 = ModulusCounter(4)
        self._mod5 = ModulusCounter(5)
        
        self._Hcounts_all    = Counter()
        self._Hcounts_paired = Counter()
        
        geneticCode = getSpeciesTranslationTable(taxId)
        self._entropyCalc = CodonEntropyCalculation(geneticCode)

    def calcSeq(self, cdsSeq):
        
        H = self._entropyCalc.calcCDSCodonEntropy(cdsSeq)

        self._Hcounts_all.update(map(roundH, H))
        
        for pairs in sequencePairsSource(cdsSeq):

            #consensus.update(pairs)
            
            self._mod3.updateWithIndices([x[0] for x in pairs])
            self._mod3.updateWithIndices([x[1] for x in pairs])
            self._mod4.updateWithIndices([x[0] for x in pairs])
            self._mod4.updateWithIndices([x[1] for x in pairs])
            self._mod5.updateWithIndices([x[0] for x in pairs])
            self._mod5.updateWithIndices([x[1] for x in pairs])
            
            # update distances histogram
            self._counts.update([x[1]-x[0] for x in pairs])

            lenH = len(H)

            for x in pairs:
                if x[0] < lenH and x[1] < lenH: # suppress pairs that fall on the stop codon
                    H0 = H[x[0]]
                    H1 = H[x[1]]

                    self._Hcounts_paired.update((roundH(H0), roundH(H1)))



    def getResults(self):
        pass


class CalcStats2(object):
    def __init__(self, taxId):

        self._taxId = taxId
        
        geneticCode = getSpeciesTranslationTable(taxId)
        self._entropyCalc = CodonEntropyCalculation(geneticCode)

        self._countsPaired = {}
        self._countsTotal = {}

        

        #self._energySum = 0.0
        #self._energyCount = 0

    def _upd(self, h, updPaired):
        s = None
        
        if updPaired:
            s = self._countsPaired
        else:
            s = self._countsTotal

            
        if h in s:
            s[h] += 1
        else:
            s[h] = 1

    def totals(self):
        return (sum(self._countsPaired.values()), sum(self._countsTotal.values()))
    
    
    def calcSeq(self, cdsSeq):
        H = self._entropyCalc.calcCDSCodonEntropy(cdsSeq)

        for h in H:
            self._upd(roundH(h), False)
            
        lenH = len(H)

        for (pairs, energy) in sequencePairsSource2(cdsSeq):

            for x in pairs:
                if x[0] < lenH and x[1] < lenH: # suppress pairs that fall on the stop codon
                    H0 = roundH(H[x[0]])
                    H1 = roundH(H[x[1]])

                    self._upd(H0, True)
                    self._upd(H1, True)

            #self._energySum += energy
            #self._energyCount += 1
            
            

    def __add__(self, other):
        assert(self._taxId == other._taxId)
        
        #self._energySum += other._energySum
        #self._energyCount += other._energyCount

        for k,v in other._countsPaired.items():
            if k in self._countsPaired:
                self._countsPaired[k] += v
            else:
                self._countsPaired[k] = v

        for k,v in other._countsTotal.items():
            if k in self._countsTotal:
                self._countsTotal[k] += v
            else:
                self._countsTotal[k] = v

        return self
                
    def __sub__(self, other):
        assert(self._taxId == other._taxId)
        
        #self._energySum -= other._energySum
        #self._energyCount += other._energyCount

        for k,v in other._countsPaired.items():
            if k in self._countsPaired:
                self._countsPaired[k] -= v
            else:
                self._countsPaired[k]  = -v

        for k,v in other._countsTotal.items():
            if k in self._countsTotal:
                self._countsTotal[k] -= v
            else:
                self._countsTotal[k]  = -v

        return self

    def __str__(self):
        return "{ paired: [%s], total: [%s] }" % (self._countsPaired, self._countsTotal)

    def __repr__(self):
        return "{ paired: [%s], total: [%s] }" % (self._countsPaired, self._countsTotal)
    
    
def testNativeDistances(fraction, taxId, numFractions):
    assert(fraction>=0)
    assert(fraction < numFractions)

    assert(type(taxId)==type(0))

    startTime = time()

    shuffler = SynonymousCodonPermutingRandomization(getSpeciesTranslationTable(taxId))

    numShuffles = 1

    numSeqsDone = 0

    diffStats = CalcStats2(taxId)
    allNativeStats = CalcStats2(taxId)
    
    for (seqId, seq) in nativeSequencesSource(taxId, fraction, numFractions):

        if random.randint(0,1)>0:
            continue

        print(seqId)

        nativeStats   = CalcStats2(taxId)
        
        nativeStats.calcSeq(seq)

        allNativeStats += nativeStats
        
        numShufflesIncluded = 0
        numAttempts = 0
        
        totalPermutationsCountForSeq = None

        while True:
            identity = None
            shuffledSeq = None

            if time() - startTime > 300:
                raise Exception("Calculation took to much time!")
            
            try:
                totalPermutationsCountForSeq, identity, shuffledSeq = shuffler.randomize(seq)

                numAttempts += 1
            except Exception as e:
                print(e)
                continue # skip this sequence

            if totalPermutationsCountForSeq < 200:
                break # skip this sequence


            if(identity < 0.95):
                
                shuffledStats = CalcStats2(taxId)
                
                shuffledStats.calcSeq(shuffledSeq)

                diff = nativeStats - shuffledStats

                diffStats += diff

                numShufflesIncluded += 1

                if numShufflesIncluded >= numShuffles:
                    break

            else:
                if numAttempts >= numShuffles * 2:
                    break

        if totalPermutationsCountForSeq < 200:
            continue # skip this sequence
        
        if numShufflesIncluded < numShuffles:
            continue
        
        numSeqsDone += 1
                
    #logging.warning(mod3.getResults())
    #logging.warning(mod4.getResults())
    #logging.warning(mod5.getResults())
    return (allNativeStats, diffStats, taxId, fraction)



# def testWithNFractions(numFractions):
#     A = scheduler.map(test, range(numFractions), taxId=436017, numFractions=numFractions)
#     B = scheduler.submit(sum, A)

#     numSequences = B.result()

#     print("%d fractions -> %d sequences" % (numFractions, numSequences))

#     return numSequences
# #newIds = testval.result()

def sign(x):
    return x/abs(x)

def plotResults(taxId, results):
    print("Plotting %d..." % taxId)

    df = pd.DataFrame({
        'entropy':pd.Categorical([]),
        'fraction':pd.Series(dtype='float')
    })

    total = sum(results._countsPaired.values())
    for k,v in results._countsPaired.items():
        v = float(v)
        df = df.append(pd.DataFrame({
            'entropy':pd.Categorical([k/10.0]),
            'fraction':pd.Series([sign(v)*v/total], dtype='float')
        }))
        
    
    fig, ax = plt.subplots()    

    df.plot.bar(x='entropy', y='fraction', ax=ax)

    plt.title(getSpeciesName(taxId))

    outputFile1 = "excess_pairs_taxid_%d.pdf" % taxId
    outputFile2 = "excess_pairs_taxid_%d.svg" % taxId
    
    plt.savefig(outputFile1)
    plt.savefig(outputFile2)
    plt.close(fig)

    


def testDistances():
    import dask

    delayedCalls = []

    fractionSize = 10

    for taxId in allSpeciesSource():

        if random.randint(0, 20) > 0:
            continue

        #if not getSpeciesProperty(taxId, 'paired-mRNA-fraction')[0] is None:
        #    continue
        outputFile1 = "excess_pairs_taxid_%d.pdf" % taxId
        outputFile2 = "excess_pairs_taxid_%d.svg" % taxId
        
        if os.path.exists(outputFile1):
            continue
        
        size = countSpeciesCDS(taxId)

        numFractions = size/fractionSize
        for i in range(numFractions):
            call = dask.delayed( testNativeDistances )(i, taxId, numFractions)
            delayedCalls.append( call )
            #taxids.append(taxId)

    print("Starting %d calls..." % len(delayedCalls))

    futures = scheduler.compute(delayedCalls) # submit all delayed calculations; obtain futures immediately

    try:
        _distributed.progress(futures) # wait for all calculations to complete
    except Exception as e:
        print(E)
    print("\n")

    print("Waiting for all tasks to complete...")
    _distributed.wait(futures)

    errorsCount = 0
    
    results = {}
    for f in futures:
        try:
            (allNativeStats, diffStats, taxId, fraction) = scheduler.gather(f)

            current = None
            if taxId in results:
                current = results[taxId]
            else:
                current = CalcStats2(taxId)

            current += diffStats

            results[taxId] = current
            
        except Exception as e:
            print(e)
            errorsCount += 1

    if errorsCount:
        print("=="* 20)
        print("Finished with %d errors!" % errorsCount)
        print("=="* 20)

    print(results)
    for taxId, result in results.items():
        plotResults(taxId, result)
    

def testDistances2():
    numFractions = 500
    
    A = scheduler.map(testNativeDistances, range(numFractions), taxId=taxId, numFractions=numFractions)
    _distributed.progress(A)
    B = scheduler.gather(A)

    #totalCounts = Counter()
    numSeqs = 0


    diffStats = CalcStats2(taxId)
    allNativeStats = CalcStats2(taxId)
    
    #mod3 = ModulusCounter(3)
    #mod4 = ModulusCounter(4)
    #mod5 = ModulusCounter(5)
    #Hcounts_all = Counter()
    #Hcounts_paired = Counter()

    for result in B:
        #totalCounts += result[1]
        #numSeqs += result[0]
        #mod3 += result[2]
        #mod4 += result[3]
        #mod5 += result[4]
        #Hcounts_all += result[5]
        #Hcounts_paired += result[6]

        allNativeStats += result[0]
        diffStats      += result[1]
    
    #printHist(totalCounts)
    #print("Num of sequences: %d" % numSeqs)
    #print("Pairs: %.4g (estimated %.3g per window)" % (sum(totalCounts), float(sum(totalCounts))/(numSeqs*20)))

    #print(mod3.getResults())
    #print(mod4.getResults())
    #print(mod5.getResults())
    #print("Hcounts - all")
    #print(Hcounts_all)
    #print("Hcounts - paired")
    #print(Hcounts_paired)
    print("All:")
    print(allNativeStats)
    print("Diffs:")
    print(diffStats)




#tt = (1, 5, 10, 20, 50, 100, 200)
#tt = (1, 5, 10, 20)

#for t in tt:
#    secs = timeit( lambda : testWithNFractions(t), number=1 )
#    print("%.3gs" % secs)

Ecoli_PolA = """ATGGTTCAGATCCCCCAAAATCCACTTATCCTTGTAGATGGTTCATCTTATCTTTATCGCGC\
ATATCACGCGTTTCCCCCGCTGACTAACAGCGCAGGCGAGCCGACCGGTGCGATGTATGGTGTCCTCAACATGCTGCG\
CAGTCTGATCATGCAATATAAACCGACGCATGCAGCGGTGGTCTTTGACGCCAAGGGAAAAACCTTTCGTGATGAACT\
GTTTGAACATTACAAATCACATCGCCCGCCAATGCCGGACGATCTGCGTGCACAAATCGAACCCTTGCACGCGATGGT\
TAAAGCGATGGGACTGCCGCTGCTGGCGGTTTCTGGCGTAGAAGCGGACGACGTTATCGGTACTCTGGCGCGCGAAGC\
CGAAAAAGCCGGGCGTCCGGTGCTGATCAGCACTGGCGATAAAGATATGGCGCAGCTGGTGACGCCAAATATTACGCT\
TATCAATACCATGACGAATACCATCCTCGGACCGGAAGAGGTGGTGAATAAGTACGGCGTGCCGCCAGAACTGATCAT\
CGATTTCCTGGCGCTGATGGGTGACTCCTCTGATAACATTCCTGGCGTACCGGGCGTCGGTGAAAAAACCGCGCAGGC\
ATTGCTGCAAGGTCTTGGCGGACTGGATACGCTGTATGCCGAGCCAGAAAAAATTGCTGGGTTGAGCTTCCGTGGCGC\
GAAAACAATGGCAGCGAAGCTCGAGCAAAACAAAGAAGTTGCTTATCTCTCATACCAGCTGGCGACGATTAAAACCGA\
CGTTGAACTGGAGCTGACCTGTGAACAACTGGAAGTGCAGCAACCGGCAGCGGAAGAGTTGTTGGGGCTGTTCAAAAA\
GTATGAGTTCAAACGCTGGACTGCTGATGTCGAAGCGGGCAAATGGTTACAGGCCAAAGGGGCAAAACCAGCCGCGAA\
GCCACAGGAAACCAGTGTTGCAGACGAAGCACCAGAAGTGACGGCAACGGTGATTTCTTATGACAACTACGTCACCAT\
CCTTGATGAAGAAACACTGAAAGCGTGGATTGCGAAGCTGGAAAAAGCGCCGGTATTTGCATTTGATACCGAAACCGA\
CAGCCTTGATAACATCTCTGCTAACCTGGTCGGGCTTTCTTTTGCTATCGAGCCAGGCGTAGCGGCATATATTCCGGT\
TGCTCATGATTATCTTGATGCGCCCGATCAAATCTCTCGCGAGCGTGCACTCGAGTTGCTAAAACCGCTGCTGGAAGA\
TGAAAAGGCGCTGAAGGTCGGGCAAAACCTGAAATACGATCGCGGTATTCTGGCGAACTACGGCATTGAACTGCGTGG\
GATTGCGTTTGATACCATGCTGGAGTCCTACATTCTCAATAGCGTTGCCGGGCGTCACGATATGGACAGCCTCGCGGA\
ACGTTGGTTGAAGCACAAAACCATCACTTTTGAAGAGATTGCTGGTAAAGGCAAAAATCAACTGACCTTTAACCAGAT\
TGCCCTCGAAGAAGCCGGACGTTACGCCGCCGAAGATGCAGATGTCACCTTGCAGTTGCATCTGAAAATGTGGCCGGA\
TCTGCAAAAACACAAAGGGCCGTTGAACGTCTTCGAGAATATCGAAATGCCGCTGGTGCCGGTGCTTTCACGCATTGA\
ACGTAACGGTGTGAAGATCGATCCGAAAGTGCTGCACAATCATTCTGAAGAGCTCACCCTTCGTCTGGCTGAGCTGGA\
AAAGAAAGCGCATGAAATTGCAGGTGAGGAATTTAACCTTTCTTCCACCAAGCAGTTACAAACCATTCTCTTTGAAAA\
ACAGGGCATTAAACCGCTGAAGAAAACGCCGGGTGGCGCGCCGTCAACGTCGGAAGAGGTACTGGAAGAACTGGCGCT\
GGACTATCCGTTGCCAAAAGTGATTCTGGAGTATCGTGGTCTGGCGAAGCTGAAATCGACCTACACCGACAAGCTGCC\
GCTGATGATCAACCCGAAAACCGGGCGTGTGCATACCTCTTATCACCAGGCAGTAACTGCAACGGGACGTTTATCGTC\
AACCGATCCTAACCTGCAAAACATTCCGGTGCGTAACGAAGAAGGTCGTCGTATCCGCCAGGCGTTTATTGCGCCAGA\
GGATTATGTGATTGTCTCAGCGGACTACTCGCAGATTGAACTGCGCATTATGGCGCATCTTTCGCGTGACAAAGGCTT\
GCTGACCGCATTCGCGGAAGGAAAAGATATCCACCGGGCAACGGCGGCAGAAGTGTTTGGTTTGCCACTGGAAACCGT\
CACCAGCGAGCAACGCCGTAGCGCGAAAGCGATCAACTTTGGTCTGATTTATGGCATGAGTGCTTTCGGTCTGGCGCG\
GCAATTGAACATTCCACGTAAAGAAGCGCAGAAGTACATGGACCTTTACTTCGAACGCTACCCTGGCGTGCTGGAGTA\
TATGGAACGCACCCGTGCTCAGGCGAAAGAGCAGGGCTACGTTGAAACGCTGGACGGACGCCGTCTGTATCTGCCGGA\
TATCAAATCCAGCAATGGTGCTCGTCGTGCAGCGGCTGAACGTGCAGCCATTAACGCGCCAATGCAGGGAACCGCCGC\
CGACATTATCAAACGGGCGATGATTGCCGTTGATGCGTGGTTACAGGCTGAGCAACCGCGTGTACGTATGATCATGCA\
GGTACACGATGAACTGGTATTTGAAGTTCATAAAGATGATGTTGATGCCGTCGCGAAGCAGATTCATCAACTGATGGA\
AAACTGTACCCGTCTGGATGTGCCGTTGCTGGTGGAAGTGGGGAGTGGCGAAAACTGGGATCAGGCGCACTAA"""

#entropyCalc = CodonEntropyCalculation(11)
#print(entropyCalc.calcCDSCodonEntropy(Ecoli_PolA))

#print(testNativeDistances(15, 243230, 300))
#print(testNativeDistances(15, taxId, 300))
testDistances()




