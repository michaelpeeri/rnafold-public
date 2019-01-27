from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from data_helpers import SpeciesCDSSource, CDSHelper, countSpeciesCDS, calcCrc
from runningstats import RunningStats
from rate_limit import RateLimit

# Configuration
taxId = 9066

statsLength = RunningStats()
statsShuffles = RunningStats()

recordsCount = 0
warningsCount = 0

rl = RateLimit(30)


total = countSpeciesCDS(taxId)

for protId in SpeciesCDSSource(taxId):
    cds = CDSHelper( taxId, protId )
    recordsCount += 1

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


    shuffles = cds.shuffledSeqIds()
    unique = len(frozenset(shuffles))

    if( unique != len(shuffles)):
        print("WARNING: duplicate shuffles found in record (taxid=%d, protId=%s); count=%d, unique=%d." % (taxId, protId, len(shuffles), unique))
        warningsCount += 1
    statsShuffles.push( len(shuffles))

    
    shuffledCrcs = set()
    shuffledTransCrcs = set()

    for shuffledSeqId in shuffles:
        shufSeq = cds.getShuffledSeq2(shuffledSeqId)
        
        shuffledCrcs.add( calcCrc(shufSeq) )
        seq1shufftrans = Seq(shufSeq, generic_dna).translate()
        shuffledTransCrcs.add( calcCrc( seq1shufftrans ) )

    print("Shuffled distinct: %d" % (len(shuffledCrcs)))
    print("Shuffled trans distinct: %d" % len(shuffledTransCrcs))
        
        
        
        

    #for s in range(len(shuffles)):
    #    shuf = cds.getShuffledSeq(s)
    #    if( len(cds.sequence()) != cds.length() ):
    #        print("WARNING: incorrect shuffled sequence length detected for record (taxid=%d, protId=%s, seqId=%d); real-length=%d, recorded-length=%d." % (taxId, protId, shuffles[s], len(shuff.sequence()), cds.length()))
    #        warningsCount += 1
        

    if(rl()):
        print("processed %d records (%.2g%%)" % (recordsCount, float(recordsCount)/total*100))

print(statsLength.count())
print("%.3g %.3g +-%.3g %.3g" % (statsLength.min(), statsLength.mean(), 2*statsLength.stdev(), statsLength.max()))

print(statsShuffles.count())
print("%.3g %.3g +-%.3g %.3g" % (statsShuffles.min(), statsShuffles.mean(), 2*statsShuffles.stdev(), statsShuffles.max()))
    
print("Done - Processed %d records, found %d warnings" % (recordsCount, warningsCount))
if( warningsCount > 0 ):
    print("%d warnings found!" % (warningsCount,))

