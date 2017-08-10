from random import randint
from data_helpers import getAllNativeCDSsForSpecies, decompressNucleicSequence, countSpeciesCDS, setSpeciesProperty, allSpeciesSource, getSpeciesProperty
import _distributed
from rnafold_vienna import RNAfoldWithStructure

# configuration
windowWidth = 40
windowStep = 10


def nativeSequencesSource(fraction, taxId, numFractions):
    for (seqId, seqData) in getAllNativeCDSsForSpecies(taxId, fraction, numFractions).items():
        cdsSeq = decompressNucleicSequence(seqData)
        del seqData
        yield (seqId, cdsSeq)


def sequencePairsSource(cdsSeq, windowWidth=windowWidth, windowStep=windowStep):
    for start in range(0, len(cdsSeq)-windowWidth-1, windowStep):
        end = start+windowWidth
        window = cdsSeq[start:end]
        assert(len(window)==windowWidth)

        (pairs, _) = RNAfoldWithStructure(window, cdsOffset=start)
        yield pairs



def calcNativePairedFraction(taxId, fraction, numFractions):
    
    countPairedNucleotides = 0
    countTotalNucleotides  = 0
    cdsCount = 0
    
    for seqId, seq in nativeSequencesSource(fraction, taxId, numFractions):
        for windowPairs in sequencePairsSource(seq):
            countTotalNucleotides += windowWidth
            countPairedNucleotides += len(windowPairs)*2
            
        cdsCount += 1


    #print("Total:  %d" % countTotalNucleotides)
    #print("Paired: %d (%.3g%%)" % (countPairedNucleotides, float(countPairedNucleotides)/countTotalNucleotides*100))

    return (taxId, fraction, cdsCount, countPairedNucleotides, countTotalNucleotides)




def runDistributed():
    import _distributed
    import dask

    scheduler = _distributed.open()

    results = {}

    #taxids = []
    delayedCalls = []

    fractionSize = 20

    for taxId in allSpeciesSource():

        if randint(0, 20) > 0:
            continue

        if not getSpeciesProperty(taxId, 'paired-mRNA-fraction')[0] is None:
            continue
        
        size = countSpeciesCDS(taxId)

        numFractions = size/fractionSize
        for i in range(numFractions):
            call = dask.delayed( calcNativePairedFraction )(taxId, i, numFractions)
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

    results = {}
    errorsCount = 0
    for f in futures:
        try:
            (taxId, fraction, cdsCount, countPairedNucleotides, countTotalNucleotides) = scheduler.gather(f)

            current = None
            if taxId in results:
                current = results[taxId]
            else:
                current = (0, 0, 0, set())

            current = (current[0] + cdsCount,
                       current[1] + countPairedNucleotides,
                       current[2] + countTotalNucleotides,
                       current[3].union( set((fraction,)) ) )

            results[taxId] = current
            
        except Exception as e:
            print(e)
            errorsCount += 1


    for taxId, result in results.items():
        if len(result[3]) != max(result[3])+1:
            #raise Exception("Found invalid number of items for taxId=%d" % taxId)
            print("Found invalid number of items for taxId=%d" % taxId)
            continue

        fraction = float(result[1]) / result[2]

        setSpeciesProperty( taxId, "paired-mRNA-fraction", "%.4g"%fraction, "computed (v3)", overwrite=False )
        
        print("TaxId: %d\t\tFraction: %.4g" % (taxId, fraction))
        

    print("Finished %d species with %d errors" % (len(results), errorsCount))
    return results
    
#test()
runDistributed()
