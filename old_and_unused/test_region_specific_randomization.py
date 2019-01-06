"""
Compare distance between paired nucleotides in native and random sequences
"""
from __future__ import print_function
from math import log10
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
import seaborn as sns; sns.set()
import pandas as pd
from data_helpers import getAllNativeCDSsForSpecies, decompressNucleicSequence, getSpeciesTranslationTable, countSpeciesCDS, getSpeciesName, allSpeciesSource, nativeSequencesSource
import _distributed
import config
from rnafold_vienna import RNAfoldWithStructure
from codon_randomization import SynonymousCodonPermutingRandomization


scheduler = _distributed.open()


def getCodonMaskForSeq(ntseq, _from, _to):
    lengthInCodons = len(ntseq)/3
    assert( len(ntseq) == lengthInCodons*3 )
    mask = [False]*lengthInCodons
    mask[_from:_to] = [True]*(_to-_from)
    return mask
    
    
def testRegionSpecificRandomization(fraction, taxId, numFractions):
    assert(fraction>=0)
    assert(fraction < numFractions)

    assert(type(taxId)==type(0))

    startTime = time()

    shuffler = SynonymousCodonPermutingRandomization(getSpeciesTranslationTable(taxId))

    numShuffles = 1

    numSeqsDone = 0

    #diffStats = CalcStats2(taxId)
    #allNativeStats = CalcStats2(taxId)

    data = []
    
    for (seqId, seq) in nativeSequencesSource(taxId, fraction, numFractions):

        if random.randint(0,1)>0:
            continue

        #print(seqId)

        #nativeStats   = CalcStats2(taxId)
        
        #nativeStats.calcSeq(seq)

        #allNativeStats += nativeStats
        
        numShufflesIncluded = 0
        numAttempts = 0
        
        totalPermutationsCountForSeq = None

        while True:
            identity = None
            shuffledSeq = None

            #if time() - startTime > 300:
            #    raise Exception("Calculation took to much time!")
            
            try:
                numAttempts += 1
                totalPermutationsCountForSeq, identity, shuffledSeq = shuffler.randomizeWithMask(seq, getCodonMaskForSeq(seq, 0, 22) )

            except Exception as e:
                print(e)
                #continue # skip this sequence
                raise e

            if numAttempts >= 3:
                    break

        data.append( (len(seq), totalPermutationsCountForSeq) )
        
        numSeqsDone += 1
                
    #logging.warning(mod3.getResults())
    #logging.warning(mod4.getResults())
    #logging.warning(mod5.getResults())
    return (taxId, fraction, numSeqsDone, data)


def testRegionSpecificRandomization2(codon, taxId):
    #assert(fraction>=0)
    #assert(fraction < numFractions)

    assert(type(taxId)==type(0))

    startTime = time()

    shuffler = SynonymousCodonPermutingRandomization(getSpeciesTranslationTable(taxId))

    numShuffles = 1

    numSeqsDone = 0

    #diffStats = CalcStats2(taxId)
    #allNativeStats = CalcStats2(taxId)

    data = []

    poolCodons = []
    
    for (seqId, seq) in nativeSequencesSource(taxId, 0, 1):
        if len(seq) >= (codon+1)*3-1:
            codon = seq[codon*3:(codon+1)*3]
            assert(len(codon)==3)
            poolCodons.append( codon )

    pool = ''.join(poolCodons)

    totalPermutationsCountForSeq = None

    numAttempts = 0

    while True:
        identity = None
        shuffledSeq = None
        
        #if time() - startTime > 300:
        #    raise Exception("Calculation took to much time!")
        
        try:
            numAttempts += 1
            totalPermutationsCountForSeq, identity, shuffledSeq = shuffler.randomize(pool)
            
        except Exception as e:
            print(e)
            #continue # skip this sequence
            raise e

        if numAttempts >= 3:
            break

    #data.append( (len(seq), totalPermutationsCountForSeq) )
        
    #numSeqsDone += 1
                
    #logging.warning(mod3.getResults())
    #logging.warning(mod4.getResults())
    #logging.warning(mod5.getResults())
    return (taxId, codon, len(pool), totalPermutationsCountForSeq)


def testRegionSpecificRandomization3(taxId):
    #assert(fraction>=0)
    #assert(fraction < numFractions)

    assert(type(taxId)==type(0))

    startTime = time()

    shuffler = SynonymousCodonPermutingRandomization(getSpeciesTranslationTable(taxId))

    numShuffles = 1

    numSeqsDone = 0

    #diffStats = CalcStats2(taxId)
    #allNativeStats = CalcStats2(taxId)

    data = []

    poolCodons = []
    allCdss = []

    anyWarnings = False

    for (seqId, seq) in nativeSequencesSource(taxId, 0, 1):
        if len(seq)%3!=0:
            print("WARNING: Bad CDS length for taxId: {} seqId: {} length: {}".format(taxId, seqId, len(seq)))
            #return(taxId, 0, [])
            #break
            anyWarnings = True
            continue

        allCdss.append(seq)
    print(len(allCdss))

    
    #(ret, identities) = shuffler.verticalPermutation(allCdss)
    return (taxId
            
    return (taxId, len(ret), identities)



# def testWithNFractions(numFractions):
#     A = scheduler.map(test, range(numFractions), taxId=436017, numFractions=numFractions)
#     B = scheduler.submit(sum, A)

#     numSequences = B.result()

#     print("%d fraNctions -> %d sequences" % (numFractions, numSequences))

#     return numSequences
# #newIds = testval.result()

def sign(x):
    return x/abs(x)

def plotResults(taxId, results):
    print("Plotting %d..." % taxId)
    #assert(len(results) > 0)
    if( len(results)==0 ): return

    df = pd.DataFrame({
        'length_nt':pd.Series(dtype='int'),
        'log_count_rand':pd.Series(dtype='float')
    })

    print(results[:20])

    for n, x in enumerate(results):
        #assert(len(x)==2)
        #assert(x[0] > 40)
        #assert(x[1] >= 1)
        df = df.append(pd.DataFrame({
            'length_nt':pd.Series([n]),
            #'log_count_rand':pd.Series([log10(x[1])], dtype='float')
            'log_count_rand':pd.Series([x], dtype='float')
        }))
        
    
    fig, ax = plt.subplots()    

    sns.jointplot( x='length_nt', y='log_count_rand', data=df )

    plt.title(getSpeciesName(taxId))

    outputFile1 = "test_region_specific_randomization_taxid_%d.pdf" % taxId
    outputFile2 = "test_region_specific_randomization_taxid_%d.svg" % taxId
    
    plt.savefig(outputFile1)
    plt.savefig(outputFile2)
    plt.close(fig)

    


def testRandomization():
    import dask

    delayedCalls = []

    fractionSize = 100000

    for taxId in allSpeciesSource():

        #if random.randint(0, 1) > 0:
        #    continue

        #if not getSpeciesProperty(taxId, 'paired-mRNA-fraction')[0] is None:
        #    continue
        #outputFile1 = "test_region_specific_randomization_taxid_%d.pdf" % taxId
        #outputFile2 = "test_region_specific_randomization_taxid_%d.svg" % taxId
        #
        #if os.path.exists(outputFile1):
        #    continue

        #testRegionSpecificRandomization3(taxId)
        #return ## DEBUG ONLY ####
        
        size = countSpeciesCDS(taxId)

        call = dask.delayed( testRegionSpecificRandomization3 )(taxId)
        delayedCalls.append( call )
        #taxids.append(taxId)

        #numFractions = size/fractionSize
        #numFractions = 100
        #for i in range(numFractions):
        #    #call = dask.delayed( testRegionSpecificRandomization )(i, taxId, numFractions)
        #    call = dask.delayed( testRegionSpecificRandomization2 )(i, taxId)
        #    delayedCalls.append( call )
        #    #taxids.append(taxId)
            
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
            #(taxId, fraction, numSeqsDone, data) = scheduler.gather(f)
            #(taxId, codon, numSeqsDone, totalPermutationsCountForSeq) = scheduler.gather(f)
            (taxId, numSeqsDone, identities) = scheduler.gather(f)
            print(numSeqsDone)
            #print(len(data))
            #data = (codon, totalPermutationsCountForSeq)

            #if taxId in results:
            #    current = results[taxId]
            #else:
            #    current = []

            #current.append( data )
            #print(len(current))

            results[taxId] = identities
            
        except Exception as e:
            print(e)
            errorsCount += 1

    if errorsCount:
        print("=="* 20)
        print("Finished with %d errors!" % errorsCount)
        print("=="* 20)

    print(results.keys())
    for taxId, result in results.items():
        print("{} :: {}".format(taxId, result))
        plotResults(taxId, result)
    


testRandomization()




