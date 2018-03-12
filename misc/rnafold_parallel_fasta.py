# Use RNAfold to calculate LFE profile for sequences read from fasta file
# Use dask to parallelize.
# Output results as csv file
from __future__ import print_function
import random
import sys
import csv
import _distributed
import config
from rnafold_vienna import RNAfold_direct
from Bio import SeqIO


fastaInput = sys.argv[1]

scheduler = _distributed.open()

   
    
def calculateNativeProfileForSequence(seqId, seq, windowWidth, windowStep):
    out = []
    for start in range(0, len(seq)-windowWidth+1, windowStep):
        end = start+windowWidth
        window = seq[start:end]
        assert(len(window)==windowWidth)

        # strip missing nucleotides
        window = str(window).replace("-", "")

        if len(window) < 1:
            out.append(0)
        else:
            energy = RNAfold_direct(window)
            out.append(energy)
        
    return( seqId, out )



    
def fastaSequencesSource(fastaIn):
    with open(fastaIn, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            yield (record.id, record.seq)


def testSequencesSource():
    yield (0, "actgcagcagtcgagtctagctagctgatcgatcgtagctagcgccggctatcgatcgtagccggcgc")
    yield (1, "actgcagcagtcgagtctagctagctgatcgatcgtagctagcgccggctatcgatcgtagccggcgc")
    yield (2, "actgcagcagtcgagtctagctagctgatcgatcgtagctagcgccggctatcgatcgtagccggcgc")


    

def calcFastaProfiles(fastaInput):
    import dask

    delayedCalls = []


    for (seqId, seq) in fastaSequencesSource(fastaInput):

        #if random.randint(0, 100) > 0:
        #    continue

        call = dask.delayed( calculateNativeProfileForSequence )(seqId, seq, 31, 1)
        delayedCalls.append( call )


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
    
    results = []
        
    for f in futures:
        try:
            (seqId, profile) = scheduler.gather(f)

            results.append((seqId, profile))

        except Exception as e:
            print(e)
            errorsCount += 1

    with open("{}.profiles.csv".format(fastaInput), "wb") as csvfile:
        csvout = csv.writer(csvfile)
        for (seqId, profile) in results:
            csvout.writerow( [seqId] + profile )
                
    if errorsCount:
        print("=="* 20)
        print("Finished with %d errors!" % errorsCount)
        print("=="* 20)

    #for taxId, result in results.items():
    #    plotResults(taxId, result)
    


calcFastaProfiles(fastaInput)




