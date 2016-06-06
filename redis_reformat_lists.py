# Convert 'rna-fold-window' fields in redis (containing lists of floats) into compact format (to save space)
# Source format example:  "-8.800000190734863" (who is in responsible for this insane format?)
# Target format example:  "8.8"  (the '-' is removed!)
#
from __future__ import print_function
import sys
import redis
import config

# configuration
windowKey = "CDS:taxid:%d:protid:%s:computed:rna-fold-window-40-0:value"
tempWindowKey = "CDS:taxid:%d:protid:%s:computed:rna-fold-window-40-0:temp-value"

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)

species = (3055, 556484)

skipped = 0
selected = 0

bytesBefore = 0
bytesAfter = 0

for taxIdForProcessing in species:
    # Iterate over all CDS entries for this taxId
    for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
        # Filters
        # Skip partial CDSs
        if(not r.exists(windowKey % (taxIdForProcessing, protId))):
            skipped += 1
            continue

        selected += 1
        # Extract the CDS sequence; write in fasta format to standard output

        origLength = r.llen( windowKey % (taxIdForProcessing, protId))

        while(True):
            origValue = r.lpop( windowKey % (taxIdForProcessing, protId))
            if( origValue == None): break

            numericVal = float(origValue)

            assert(numericVal <= 0.0)
            newValue = "%.3g" % (-numericVal)  # Store values as positive! (to save space on the "-"...)

            bytesBefore += len(origValue)
            bytesAfter += len(newValue)
            if( len(origValue) < len(newValue) ):
                print("Warning: length increased for '%s' -> '%s'" % (origValue, newValue) )
            assert(len(origValue) >= len(newValue) )

            r.rpush( tempWindowKey % (taxIdForProcessing, protId), newValue )

        assert( r.llen( tempWindowKey % (taxIdForProcessing, protId)) == origLength )
        r.rename( tempWindowKey % (taxIdForProcessing, protId), windowKey % (taxIdForProcessing, protId) )

if( bytesBefore>0 ):
    print("Converted %d keys; Saved %d bytes (%.2f%%)" % (selected, bytesBefore-bytesAfter, float(bytesBefore-bytesAfter)/bytesBefore*100 ) )
