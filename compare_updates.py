# Examine data in the "pending updates" table (sequence_series2_updates)
# Try to identify errors, inconsistencies and duplications
# 
from collections import Counter
import argparse
from data_helpers import seriesUpdatesSource, splitLongSequenceIdentifier
from rate_limit import RateLimit


argsParser = argparse.ArgumentParser()
argsParser.add_argument("--verbose", type=int, default=0)
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

err = Counter()
recordsByTaxId = Counter()
badUpdateRecords = set()
badOriginalRecords = set()

def extractRecordId(longIdentifier):
    return longIdentifier.split(":")

for updates in seriesUpdatesSource(102, 500):
    total += len(updates)

    for updateRecord, originalRecord in updates:
        assert(len(updateRecord)==3)
        assert(len(originalRecord)==2)

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
