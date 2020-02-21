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
# Enqueue failed entries (i.e. entries missing the calculation results) for repeated processing.
# Scan all CDS entries from a given taxid; Find entries missing a computation result; Re-insert them into the queue for repeated processing
# Input - taxid1,taxid2,taxid3
#       - computationTag
#       - randomFraction
# Example:
# python2 requeue_sequences_missing_energies_for_sliding_window.py 3055,556484 rna-fold-window-40-0 10
# 
# TODO: This script will requeue sequences that have already in the queue but haven't been completed yet.
# TODO: Add support for step-size >1
import sys
import codecs
from random import randint
#import redis
import config
import mysql_rnafold as db
#from sqlalchemy import sql
#from sqlalchemy.sql.expression import func
from data_helpers import CDSHelper, countSpeciesCDS, getSpeciesName, SpeciesCDSSource, numItemsInQueue, RateLimit, getAllComputedSeqs


# command-line arguments
species = map(int, sys.argv[1].split(","))
#computationTag = sys.argv[2]
#if( computationTag.find(':') != -1 ): raise Exception("computation tag cannot contain ':' (should be compatible with redis key names)")
# e.g. rna-fold-0
#randomFraction = int(sys.argv[3])

# Configuration
#queueKey = "queue:tag:awaiting-%s:members" % computationTag
#nativeCdsSeqIdKey = "CDS:taxid:%d:protid:%s:seq-id"
#seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"
windowWidth = 40
calculationWidth = 150
seriesSourceNumber = db.Sources.RNAfoldEnergy_SlidingWindow40
computationTag = "rna-fold-window-40-0"
# TODO: Add support for step-size >1

# Establish DB connections
#r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
#session = db.Session()

skipped = 0
selected = 0
alreadyCompleted = 0


rl = RateLimit(30)

computed = getAllComputedSeqs(seriesSourceNumber)
print("Collecting data from %d computation results..." % len(computed))


match = {}

for taxIdForProcessing in species:
    print("Procesing %d sequences for tax-id %d (%s)..."
          % (countSpeciesCDS(taxIdForProcessing),
             taxIdForProcessing,
             getSpeciesName(taxIdForProcessing)))

    # Iterate over all CDS entries for this species
    for protId in SpeciesCDSSource(taxIdForProcessing):
        #protId = codecs.decode(protId)
        # Filtering


        # Skip sequences with partial CDS annotations
        #if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
        #    skipped += 1
        #    continue

        #if( not r.exists(nativeCdsSeqIdKey % (taxIdForProcessing, protId)) ):
        #    skipped +=1
        #    continue

        if(rl()):
            print("%d %d" % (selected, skipped))

        cds = CDSHelper(taxIdForProcessing, protId)

        seqLength = cds.length()
        if( not seqLength is None ):
            # Skip sequences with length <40nt (window width)
            if(seqLength < calculationWidth + windowWidth - 1 ):
                skipped += 1
                continue
        else:
            print("Warning: Could not find CDS length entry for taxid=%d, protid=%s" % (taxIdForProcessing, protId) )
            skipped += 1
            continue

        #requiredNumWindows = seqLength - windowWidth + 1

        cdsSeqId = cds.seqId()

        match[cdsSeqId] = (protId,-1)

        shuffledIds = cds.shuffledSeqIds()

        computedShufflesCount = 0
        missingShuffles = []

        for n in range(50):
            match[shuffledIds[n]] = (protId,n)
            if shuffledIds[n] in computed:
                computedShufflesCount += 1
            else:
                missingShuffles.append(n)

        if( computedShufflesCount<50 ):
            print("%s - found only %d computed results" % (protId,computedShufflesCount))

            #print("Enqueueing...")
            #for shuf in missingShuffles:
            #    cds.enqueueForProcessing(computationTag, shuf)

            continue

        if( computedShufflesCount == 50 and cdsSeqId in computed ):
            print("%s - OK" % (protId,))
            #print(cds.getCalculationResult( seriesSourceNumber, -1))
        
        selected += 1
        alreadyCompleted += 1

print(len(match))


print("%d selected, %d skipped, %d already completed (%d total)" % (selected, skipped, alreadyCompleted, selected+skipped))
#print("queue contains %d items" % numItemsInQueue(computationTag))
