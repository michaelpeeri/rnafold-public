import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from data_helpers import SpeciesCDSSource, CDSHelper, countSpeciesCDS, calcCrc, DropSequenceWithResults
import data_helpers # REMOVE THIS!!!
from rate_limit import RateLimit
from runningstats import RunningStats

# Configuration
taxId = int(sys.argv[1])

statsLength = RunningStats()
statsShuffles = RunningStats()

recordsCount = 0
warningsCount = 0

rl = RateLimit(30)


total = countSpeciesCDS(taxId)

seqDeleter = DropSequenceWithResults()

for protId in SpeciesCDSSource(taxId):
    cds = CDSHelper( taxId, protId )
    recordsCount += 1

    print(protId)

    statsLength.push( cds.length())

    if( len(cds.sequence()) != cds.length() ):
        print("WARNING: incorrect sequence length detected for record (taxid=%d, protId=%s); real-length=%d, recorded-length=%d." % (taxId, protId, len(cds.sequence()), cds.length()))
        warningsCount += 1

    recomputedCrc = calcCrc( cds.sequence() )
    annotatedCrc = cds.crc()
    assert( recomputedCrc == annotatedCrc )
    print(cds.sequence()[:15])

    seq1trans = Seq(cds.sequence(), generic_dna).translate()
    crc1 = calcCrc(seq1trans)
    

    print(data_helpers.r.lrange(data_helpers.shuffledSeqIdsKey % (taxId, protId), 0, 1))
    shuffles = cds.shuffledSeqIds()
    unique = len(frozenset(shuffles))

    if( unique != len(shuffles)):
        print("WARNING: duplicate shuffles found in record (taxid=%d, protId=%s); count=%d, unique=%d." % (taxId, protId, len(shuffles), unique))
        warningsCount += 1
    statsShuffles.push( len(shuffles))

    
    shuffledCrcs = set()
    shuffledTransCrcs = set()

    shuffledCrcsToIds = {}


    for shuffledSeqId in shuffles:
        shufSeq = cds.getShuffledSeq2(shuffledSeqId)

        shuffledCrc = calcCrc(shufSeq)
        
        shuffledCrcs.add( shuffledCrc )

        if shuffledCrc in shuffledCrcsToIds:
            shuffledCrcsToIds[shuffledCrc].append( shuffledSeqId )
        else:
            shuffledCrcsToIds[shuffledCrc] = [ shuffledSeqId ]
        
        seq1shufftrans = Seq(shufSeq, generic_dna).translate()
        shuffledTransCrcs.add( calcCrc( seq1shufftrans ) )

    duplicateSeqs = 0
    for shuffledCrc, shuffleIds in shuffledCrcsToIds.items():
        if( len(shuffleIds) > 1 ):  # redundant sequence found
            #print("%s: Ids [%s] have the same Crc %x" % (cds.getProtId(), shuffleIds, shuffledCrc))

            # Drop all but the first of the redundant sequences
            for seqId in shuffleIds[1:]:
                seqDeleter.markSequenceForDropping(taxId, protId, seqId)
                duplicateSeqs += len(shuffleIds)-1
                
    if( duplicateSeqs ):
        print("(debug) %d duplicates found for %s" % (duplicateSeqs, protId))
    else:
        print("(debug) No duplicates found for %s" % protId)

    print("Shuffled distinct: %d" % (len(shuffledCrcs)))
    print("Shuffled trans distinct: %d" % len(shuffledTransCrcs))
        

    #for s in range(len(shuffles)):
    #    shuf = cds.getShuffledSeq(s)
    #    if( len(cds.sequence()) != cds.length() ):
    #        print("WARNING: incorrect shuffled sequence length detected for record (taxid=%d, protId=%s, seqId=%d); real-length=%d, recorded-length=%d." % (taxId, protId, shuffles[s], len(shuff.sequence()), cds.length()))
    #        warningsCount += 1
        

    if(rl()):
        print("processed %d records (%.2g%%)" % (recordsCount, float(recordsCount)/total*100))

seqDeleter.performDroppingNow()  # Drop any remaining redundant sequences
        
print(statsLength.count())
print("%.3g %.3g +-%.3g %.3g" % (statsLength.min(), statsLength.mean(), 2*statsLength.stdev(), statsLength.max()))

print(statsShuffles.count())
print("%.3g %.3g +-%.3g %.3g" % (statsShuffles.min(), statsShuffles.mean(), 2*statsShuffles.stdev(), statsShuffles.max()))
    
print("Done - Processed %d records, found %d warnings" % (recordsCount, warningsCount))
if( warningsCount > 0 ):
    print("%d warnings found!" % (warningsCount,))

