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
# Transfer sliding-window rna-fold energies from redis to mysql
# Input: 
#
from __future__ import print_function
import sys
import redis
import config
import mysql_rnafold as db
from sqlalchemy import sql
from sqlalchemy.sql.expression import func


# Configuration
species = (3055, 556484)
windowWidth = 40
windowKey = "CDS:taxid:%d:protid:%s:computed:rna-fold-window-40-0:value"
seqIdKey = "CDS:taxid:%d:protid:%s:seq-id"
seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"

# Connect to redis (source DB)
r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
# Connect to MySQL (destination DB)
session = db.Session()


skipped = 0
selected = 0

for taxIdForProcessing in species:
    # Iterate over all CDS entries for this taxId
    for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
        # Filters

        # Skip sequences with no sliding-window energies 
        if( not r.exists(windowKey % (taxIdForProcessing, protId) ) ):
            skipped += 1
            continue

        # Skip sequences with no sequence-id
        currentSequenceId = r.get( seqIdKey % (taxIdForProcessing, protId) )
        if( currentSequenceId is None):
            print("Warning: skipping sequence (taxid=%d, prot-id=%s) because it doesn't have a sequence-id" % (taxIdForProcessing, protId))
            skipped += 1
            continue
        currentSequenceId = int(currentSequenceId)

        # Get the length of the source series
        origLength = r.llen( windowKey % (taxIdForProcessing, protId))
        if( origLength is None or origLength==0 ):
            skipped += 1
            continue;

        # get the length of the CDS sequence
        sequenceLength = db.connection.execute( sql.select(( sql.func.length(db.sequences.c.sequence),)).where(db.sequences.c.id==currentSequenceId) ).scalar()

        expectedLength = sequenceLength - windowWidth + 1
        if( origLength != expectedLength ):
            print("Warning: skipping sequence (taxid=%d, prot-id=%s) because it has the wrong number of elements (%d, expected %d)" % (taxIdForProcessing, protId, origLength, expectedLength))
            skipped += 1
            continue

        selected += 1

        # Move each one of the values to the destination
        nextValueId = 0
        while(True):
            origValue = r.lpop( windowKey % (taxIdForProcessing, protId))
            if( origValue == None): break

            # Convert the value.
            # NOTE: The values should be negative, but where stored without the '-' to save space...
            numericVal = -float(origValue)
            assert(numericVal <= 0.0) # sanity test for current value

            f  = db.SequenceSeries(
                taxid=taxIdForProcessing,
                source=db.Sources.RNAfoldEnergy_SlidingWindow40,
                index=nextValueId,
                sequence_id=currentSequenceId,
                value=numericVal )
            session.add(f)
            nextValueId += 1
        session.commit() # add the new items in the series
        

        # Side task: store the length of the sequence with the metadata (as it's no longer there)
        assert(not r.exists( seqLengthKey % (taxIdForProcessing, protId) ) )
        r.set( seqLengthKey % (taxIdForProcessing, protId), sequenceLength )

        assert( nextValueId == origLength )
        r.delete( windowKey % (taxIdForProcessing, protId), windowKey % (taxIdForProcessing, protId) )

        if( selected % 500 == 499 ):
            print("Processed %d sequences (%d skipped)" % (selected, skipped) )


print("Done - moved %d sequences." % selected)
