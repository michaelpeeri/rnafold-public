# Scan all seqs for a CDS; delete all associated shuffled sequences from redis+mysql
#
import argparse
from data_helpers import SpeciesCDSSource, CDSHelper, countSpeciesCDS, calcCrc
from rate_limit import RateLimit
#from runningstats import RunningStats
from runningstats import OfflineStats



argsParser = argparse.ArgumentParser()
argsParser.add_argument("--taxId", type=int, required=True)
argsParser.add_argument("--keep-first-n-shuffles", type=int, default=None)
args = argsParser.parse_args()

# Configuration
taxId = args.taxId

#statsShuffles = RunningStats()
statsShuffles = OfflineStats()

recordsCount = 0
warningsCount = 0

rl = RateLimit(30)

total = countSpeciesCDS(taxId)

for protId in SpeciesCDSSource(taxId):
    cds = CDSHelper( taxId, protId )

    statsShuffles.push( cds.dropShuffledSeqs(lastItemToKeep=args.keep_first_n_shuffles) )

    recordsCount += 1

    if(rl()):
        print("processed %d records (%.2g%%)" % (recordsCount, float(recordsCount)/total*100))

    # DEBUG ONLY # DEBUG ONLY # DEBUG ONLY # DEBUG ONLY # DEBUG ONLY # DEBUG ONLY #
    #if( recordsCount > 20 ):
    #    break 
    # DEBUG ONLY # DEBUG ONLY # DEBUG ONLY # DEBUG ONLY # DEBUG ONLY # DEBUG ONLY #

print(statsShuffles.count())
print("%.3g %.3g +-%.3g %.3g" % (statsShuffles.min(), statsShuffles.mean(), 2*statsShuffles.stdev(), statsShuffles.max()))
    
print("Done - Processed %d records, found %d warnings" % (recordsCount, warningsCount))

