# Examine data in the "pending updates" table (sequence_series2_updates)
# Try to identify errors, inconsistencies and duplications
# 
from collections import Counter
import argparse
from hashlib import md5
from data_helpers import seriesUpdatesSource, splitLongSequenceIdentifier
from rate_limit import RateLimit
import mysql_rnafold as db
from sqlalchemy import sql
from sqlalchemy.sql.expression import func


argsParser = argparse.ArgumentParser()
argsParser.add_argument("--verbose", type=int, default=0)
argsParser.add_argument("--groupSize", type=int, default=800)
argsParser.add_argument("--find-duplicates", type=bool, default=False)
argsParser.add_argument("--perform-delete", type=bool, default=False)
args = argsParser.parse_args()


rl = RateLimit(30)


total = 0
class ErrorTypes:
    MergeIdentityErrorCount = 1
    UpdateRecordIdentityMismatch = 2
    OriginalRecordIdentityMismatch = 3
    UpdateRecordFormatError = 4
    UpdateProfileTooShort = 5
    ExistingProfileValueCorrupted = 6
    UpdateProfileIsIdentical = 7
    SeqCRCchanged = 8
    IdentifiterInsideRecordChanged = 9
    ConflictingOrDuplicateUpdateRecords = 10

err = Counter()
recordsByTaxId = Counter()
badUpdateRecords = set()
badOriginalRecords = set()

def extractRecordId(longIdentifier):
    return longIdentifier.split(":")

def allUpdatesSeriesSource():
    A = db.sequence_series2_updates.alias(name="A")

    series = db.connection.execute( sql.select((func.distinct(A.c.source),)) ).fetchall()

    # yield 505 # test only
    
    for s in [x[0] for x in series]:
        yield s
    

def allUpdatesSource(series=None):

    A = db.sequence_series2_updates.alias(name="A")
    #B = db.sequence_series2_updates.alias(name="B")

    maxIndex = None
    
    if series is None:
        maxIndex = db.connection.execute( sql.select((func.max(A.c.dummy_id),)) ).fetchone()[0]
    else:
        maxIndex = db.connection.execute( sql.select((func.max(A.c.dummy_id),)).where(A.c.source==series) ).fetchone()[0]
        
    if maxIndex is None:
        print("No updates found")
        return
    
    print("Total update records: %d" % maxIndex)
    #maxIndex = 4949
    groupSize = 5000


    total = 0
    for fromId in range(0, maxIndex, groupSize):
        q = None
        if series is None:
            q = sql.select(( A.c.dummy_id, A.c.sequence_id, func.length(A.c.content))).where(
                sql.and_(
                    #                A.c.dummy_id < B.c.dummy_id,
                    #                A.c.sequence_id == B.c.sequence_id,
                    #                A.c.source == B.c.source,
                    sql.between(A.c.dummy_id, fromId, min(fromId+groupSize-1, maxIndex))
                )).order_by(A.c.dummy_id)
        else:
            q = sql.select(( A.c.dummy_id, A.c.sequence_id, func.length(A.c.content))).where(
                sql.and_(
                    #                A.c.dummy_id < B.c.dummy_id,
                    #                A.c.sequence_id == B.c.sequence_id,
                    A.c.source == series,
                    sql.between(A.c.dummy_id, fromId, min(fromId+groupSize-1, maxIndex))
                )).order_by(A.c.dummy_id)

        #print("Executing query for range(%d, %d)..." % (fromId, fromId+groupSize))
        result = db.connection.execute( q )
        pairRecords = result.fetchall()

        for r in pairRecords:
            yield r
        
        total += len(pairRecords)

    print("Total: %d" % total)
    

    
"""
Return redundant update records, i.e., those updating the same series (computation) for the same sequence
"""
def findDuplicateUpdates():
    recordsToDelete = set()

    prevRecords = {}

    # Iterate of all duplicate pairs, i.e., pairs of update records that apply to the same (sequence_id, source)
    for series in allUpdatesSeriesSource():

        print("Processing all updates for series %d..." % series)
        
        for (idB, sequenceId, lengthB) in allUpdatesSource(series):

            if not sequenceId in prevRecords:
                prevRecords[sequenceId] = (idB, lengthB)
                continue

            (idA, lengthA) = prevRecords[sequenceId]
            assert(idA < idB)


            recordToDelete = None

            # Try to delete the smaller record (by compressed byte size) - since we are about to destroy data, try to destroy the least
            # if both records have the same size, delete the newest
            # 
            # Note: since any pair may be part of a larger group of duplicate records, this rule must define a strong ordering on the members of the
            #       group. That way, only one member (the "largest" one) will be maintained after a single iteration.
            #
            if lengthA is None:
                if lengthB is None:
                    recordToDelete = idB
                else: # lengthA is None and lengthB is not None:
                    recordToDelete = idA
            elif lengthB is None:
                recordToDelete = idB
            elif lengthA is None and lengthB is None:
                recordToDelete = idB
            elif lengthA > lengthB:
                recordToDelete = idB
            elif lengthB > lengthA:
                recordToDelete = idA
            else:
                recordToDelete = idB

            if( recordToDelete==idA ):
                print("Found duplicate pair: X%d %d" % (idA, idB))
            elif( recordToDelete==idB):
                print("Found duplicate pair: %d X%d" % (idA, idB))
            else:
                assert(False)

            recordsToDelete.add(recordToDelete)

    return recordsToDelete
        

        
        
if( args.find_duplicates ):
    duplicates = findDuplicateUpdates()
    print(len(duplicates))
    duplicatesList = sorted(list(duplicates))
    print(md5(str(duplicatesList)).hexdigest())
    badUpdateRecords.update(duplicates)
    if len(duplicates) > 0:
        err[ErrorTypes.ConflictingOrDuplicateUpdateRecords] += len(duplicates)


        
    
for updates in seriesUpdatesSource(102, 500):
    total += len(updates)

    for updateRecord, originalRecord in updates:
        assert(len(updateRecord)==3)
        assert((originalRecord is None) or (len(originalRecord)==2))

        if( args.verbose>2 ):
            print("================================================================")
            print("<<< %s" % str(originalRecord))
            print("---")
            print(">>> %s" % str(updateRecord))

        # Update identity test
        externalId = updateRecord[0]
        internalLongId = updateRecord[1]['id']
        try:
            internalId = splitLongSequenceIdentifier(internalLongId)
            internalSeqId = internalId[2]

            if(internalSeqId != externalId ):
                err[ErrorTypes.UpdateRecordIdentityMismatch] += 1
                badUpdateRecords.add(updateRecord[2])
                if( args.verbose>2 ):
                    print("UpdateRecordIdentityMismatch  (internalSeqId '%d' != externalId '%d')" % (internalSeqId, externalId))

            recordsByTaxId[internalId[0]] += 1

        except Exception as e:
            err[ErrorTypes.UpdateRecordFormatError] += 1
            badUpdateRecords.add(updateRecord[2])
            print(e)
            

        if not originalRecord is None:
            # Tests requiring and original record

            # Original identity test
            if(originalRecord[0] != splitLongSequenceIdentifier(originalRecord[1]['id'])[2]):
                # Note: This is an error in the original record. How to treat this?
                err[ErrorTypes.OriginalRecordIdentityMismatch] += 1
                badOriginalRecords.add(originalRecord[0])
                if( args.verbose>2 ):
                    print("OriginalRecordIdentityMismatch (originalRecord[0] '%d' != internalId '%d')" % (originalRecord[0], splitLongSequenceIdentifier(originalRecord[1]['id'])[2]))
            
            # records match test
            if( updateRecord[0] != originalRecord[0] ):
                #print(originalRecord)
                #print(updateRecord)
                err[ErrorTypes.MergeIdentityErrorCount] += 1
                badUpdateRecords.add(updateRecord[2])
                if( args.verbose>2 ):
                    print("MergeIdentityErrorCount")

            if( updateRecord[1]['seq-crc'] != originalRecord[1]['seq-crc'] ):
                err[ErrorTypes.SeqCRCchanged] += 1
                badUpdateRecords.add(updateRecord[2])
                if( args.verbose>2 ):
                    print("SeqCRCchanged")

            if( updateRecord[1]['id'] != originalRecord[1]['id'] ):
                err[ErrorTypes.IdentifiterInsideRecordChanged] += 1
                badUpdateRecords.add(updateRecord[2])
                if( args.verbose>2 ):
                    print("IdentifiterInsideRecordChanged")

            profile0 = originalRecord[1]['MFE-profile']
            profile1 = updateRecord[1]['MFE-profile']
            if( len(profile1) < len(profile0) ):
                err[ErrorTypes.UpdateProfileTooShort] += 1
                badUpdateRecords.add(updateRecord[2])
                if( args.verbose>2 ):
                    print("UpdateProfileTooShort")
            else:
                #print('-------- 0 --------')
                #print(len(profile0))
                #print(profile0)
                #print('-------- 1 --------')
                #print(len(profile1))
                #print(profile1)
                profile0.extend([None]*(len(profile1)-len(profile0))) # Add 'None's at the end (to allow comparison of new values)
                hasNewWindows = False
                
                for v0, v1 in zip( profile0, profile1 ):
                    # Note: this check rejects changes to existing windows (though those might be needed at some point)
                    if( (not v0 is None) and
                        ( (v1 is None) or
                          (v1 > 0.0) or
                          (abs(v0-v1) >= 1e-8) )
                        ):
                        err[ErrorTypes.ExistingProfileValueCorrupted] += 1
                        badUpdateRecords.add(updateRecord[2])
                        if( args.verbose>2 ):
                            print("ExistingProfileValueCorrupted")
                            
                        break
                    if( (v0 is None) and (not v1 is None) ): # Is this a new value
                        hasNewWindows = True
                        
                if( not hasNewWindows ):  # Did we find any new values?
                    err[ErrorTypes.UpdateProfileIsIdentical] += 1
                    badUpdateRecords.add(updateRecord[2])
                    if( args.verbose>2 ):
                        print("UpdateProfileIsIdentical")
                    
                    
        
    if(rl()):
        print(total, err, recordsByTaxId)
        
    #if( total > 90000):
    #    print("DEBUG: Stopped before end!")
    #    break



    
print(total)
print(err)
print(recordsByTaxId)

print("-- To delete %d invalid records: --" % len(badUpdateRecords))
print("use rnafold;")
todelete = list(badUpdateRecords)
firstItem = 0
while(True):
    lastItem = min(firstItem+500, len(todelete))
    if( lastItem <= firstItem ):
        break
    print("-- items %d - %d: --" % (firstItem, lastItem))
    print("delete from sequence_series2_updates")
    print("where dummy_id in (%s);" % (", ".join(map(str, todelete[firstItem:lastItem]))))
    print("go")
    firstItem = lastItem
    
if(badOriginalRecords):
    print("WARNING: Found error in original records.")
    print("WARNING: Original records: ")
    print(badOriginalRecords)

print(total, err, recordsByTaxId)

if( badUpdateRecords and args.perform_delete ):
    response = raw_input("Enter 'yes' to perform deletion of %d updates (or anything else to skip): " % len(badUpdateRecords))
    if( response=="yes" ):
        todelete = list(badUpdateRecords)
        firstItem = 0
        rowsDeleted = 0

        bulkSize = 1000
        
        while(True):
            lastItem = min(firstItem+bulkSize, len(todelete))
            if( lastItem <= firstItem ):
                break
            
            stmt = db.sequence_series2_updates.delete().where( db.sequence_series2_updates.c.dummy_id.in_( todelete[firstItem:lastItem] ) )
            #print(stmt)
            
            result = db.connection.execute( stmt )

            count = result.rowcount
            result.close()
            #count = lastItem-firstItem-1  # test only

            #print("%d <= %d-%d" % (count, lastItem, firstItem))
            assert( count <= lastItem-firstItem )
            rowsDeleted += count
            
            firstItem = lastItem


        print("Rows deleted: %d (expected: %d)" % (rowsDeleted, len(badUpdateRecords)))

