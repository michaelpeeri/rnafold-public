# Transfer sliding-window rna-fold energies from redis to mysql
# Input: 
#
from __future__ import print_function
import sys
import redis
import numpy
from math import sqrt
import numpy as np
from scipy import stats
from scipy.stats import norm
#import matplotlib
#matplotlib.use("cairo")
#import matplotlib.pyplot as plt
import config
import mysql_rnafold as db
from sqlalchemy import sql
from sqlalchemy.sql.expression import func
import runningstats


# Configuration
calculationWidth = 200
species = (3055, 556484)
windowWidth = 40
seriesSourceNumber = db.Sources.RNAfoldEnergy_SlidingWindow40
seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"
cdsSeqIdKey = "CDS:taxid:%d:protid:%s:seq-id"
shuffledSeqIdsKey = "CDS:taxid:%d:protid:%s:shuffled-seq-ids-v2"
partialCDSKey = "CDS:taxid:%d:protid:%s:partial"
debugMode = True

# Connect to redis (source DB)
r = redis.StrictRedis(host=config.host, port=config.port, db=config.db, password=config.password)
# Connect to MySQL (destination DB)
session = db.Session()
#db.db.echo = True



skipped = 0
selected = 0


#        prevWindowCount = db.connection.execute( sql.select(( sql.func.count('*'),)).select_from(db.sequence_series).where(
#                sql.and_(
#                    db.sequence_series.c.sequence_id==seqId,
#                    db.sequence_series.c.source==seriesSourceNumber,
#                    db.sequence_series.c.ext_index==0
#                    )) ).scalar()


def checkWindows(windows):
    if( not debugMode ): return

    lastIndex = -1
    for vals in windows:
        currIndex = vals[1]
        assert( currIndex == lastIndex + 1 )
        lastIndex = currIndex
    assert(lastIndex <= calculationWidth-1)

#xaxis = np.arange(0,200,1)
#def plotCdsMFE(identifier, windows):
#
#    fig = plt.figure()
#    plt.plot(xaxis, windows)
#    plt.xlabel('window start (nt, from cds)')
#    plt.ylabel('MFE')
#    plt.grid(True)
#    plt.savefig("mfe_40nt_cds_%s.png" % identifier )
#    plt.close(fig)

statsLines = []

for taxIdForProcessing in species:
    # Iterate over all CDS entries for this taxId

    cdsStats = []
    shuffledStats = []
    for i in range(calculationWidth):
        cdsStats.append( runningstats.RunningStats(debug=True) )
        shuffledStats.append( runningstats.RunningStats(debug=True) )
        
    for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
        # Filters

        # 
        if( r.exists(partialCDSKey % (taxIdForProcessing, protId) ) ):
            skipped += 1
            continue

        seqLength = r.get(seqLengthKey % (taxIdForProcessing, protId))
        if( not seqLength is None ):
            # Skip sequences too short to have a sufficient number of windows
            seqLength = int(seqLength)
            if(seqLength < calculationWidth + windowWidth):
                skipped +=1
                continue
        else:
            print("Warning: Could not find CDS length entry for taxid=%d, protid=%s" % (taxIdForProcessing, protId) )
            skipped += 1
            continue

        requiredNumWindows = seqLength - windowWidth + 1


        straightCdsSeqId = r.get(cdsSeqIdKey % (taxIdForProcessing, protId))
        if( straightCdsSeqId is None ):
            continue
        straightCdsSeqId = int(straightCdsSeqId)


        # get main seq
        straightWindows = db.connection.execute( sql.select(( db.sequence_series.c.value, db.sequence_series.c.index)).select_from(db.sequence_series).where(
                sql.and_(
                    db.sequence_series.c.sequence_id==straightCdsSeqId,
                    db.sequence_series.c.source==seriesSourceNumber,
                    db.sequence_series.c.ext_index==0,
                    db.sequence_series.c.index < calculationWidth
                    )).order_by(sql.asc(db.sequence_series.c.index))).fetchall()
        # returns array of 1-tupples
        if(len(straightWindows) < calculationWidth):
            skipped +=1
            continue

        # Save the first N running windows calculated for this protein
        straightWindows = straightWindows[:calculationWidth]
        # straightWindows is a 200-array of 1-tuples (i.e. the MFE values for each window)
        #
        # Make sure the number
        checkWindows(straightWindows)
        print(straightWindows)


        # get shuffled seqs
        shuffledSeqIds = map(int, r.lrange(shuffledSeqIdsKey % (taxIdForProcessing, protId), 0, -1))

        fullSeries = []

        for idx, currentSeqId in enumerate(shuffledSeqIds):
            # Count the number of window results stored for this series (in the current sequence)

            shuffledWindows = db.connection.execute( sql.select(( db.sequence_series.c.value, db.sequence_series.c.index)).select_from(db.sequence_series).where(
                    sql.and_(
                        db.sequence_series.c.sequence_id==currentSeqId,
                        db.sequence_series.c.source==seriesSourceNumber,
                        db.sequence_series.c.ext_index==0,
                        db.sequence_series.c.index < calculationWidth
                        )).order_by(sql.asc(db.sequence_series.c.index))).fetchall()
            # returns array of 1-tupples
            if(len(shuffledWindows) < calculationWidth):
                continue

            # Save the first N running windows calculated for this protein
            shuffledWindows = shuffledWindows[:calculationWidth]
            checkWindows(shuffledWindows)

            fullSeries.append( shuffledWindows )



        # Disable plotting
        #plotCdsMFE( protId, straightWindows )

        if( len(fullSeries) < 10 ):
            skipped += 1
            continue

        a = []

        # 
        for i in range(calculationWidth):
            # array 1-tuple
            vs = straightWindows[i][0]
            vt = []
            
            for n in range(10):
                # 10-array calculationWidth-array 1-tuple 
                vt.append( fullSeries[n][i][0] ) 
            
            n = numpy.array(vt)
            sMean = n.mean()
            sStd = n.std()
            if( sStd < 1e-10 ):
                # TODO: HANDLE THIS CASE BETTER
                break

            z = (vs - sMean) / sStd ################################
            #a.append(z)
            # Testing - store the original window energy (not z-score)
            a.append(vs)

        if( len(a) != calculationWidth):
            skipped += 1
            continue


        # CDS selected - No more skipping after this point...
        selected += 1

        # Update the cumulative stats for this protein
        for i in range(calculationWidth):
            # array 1-tuple
            cdsStats[i].push(straightWindows[i][0])

            for n in range(10):
                # 10-array calculationWidth-array 1-tuple 
                shuffledStats[i].push( fullSeries[n][i][0] )

        # Emit the values for downstream processing
        print("%d,%s,%s" % (taxIdForProcessing, protId, ",".join(map(lambda x: "%.4g"%x, a))))

        # Display status messages
        #if( selected % 500 == 499 ):
        #    print("#Processed %d sequences (%d skipped)" % (selected, skipped) )

        # DEBUG ONLY
        #if( selected > 1005 ):
        #    break

    # Finished processing all CDSs for this species

    # Make sure cdsStats for all windows is based on the same number of sequences
    prev = -1
    curr = -1
    for i in range(calculationWidth):
        curr = cdsStats[i].count()
        if( prev > 0 ):
            assert(curr==prev)
        else:
            prev = curr
    print("#cdsStats elements contain %d items" % curr)

    # Make sure cdsStats for all windows is based on the same number of sequences
    prev = -1
    for i in range(calculationWidth):
        curr = shuffledStats[i].count()
        if( prev > 0 ):
            assert(curr==prev)
        else:
            prev = curr
    print("#shuffledStats elements contain %d items" % curr)

    # Generate stats line (unless there are too few results)
    if( cdsStats[0].count() >= 1000 and shuffledStats[0].count() >= 1000 ):
        out = ["#" + str(taxIdForProcessing)]
        for i in range(calculationWidth):
            compositeZ = (cdsStats[i].mean() - shuffledStats[i].mean()) / shuffledStats[i].stdev()
            out.append("%.3g" % compositeZ)
        statsLines.append( ",".join(out) )
    else:
        print("Warning: Skipping stats-line for species %d, because there aren't enough result..." % taxIdForProcessing)


print("Done - processed %d sequences (skipped %d sequences)." % (selected, skipped))
for l in statsLines:
    print(l)
