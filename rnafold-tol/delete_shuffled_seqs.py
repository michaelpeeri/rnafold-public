# RNAFold - Analyze mRNA folding bias (Local Folding Energy) based on randomizations.
# Copyright (C) 2016-2020 Michael Peeri
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

