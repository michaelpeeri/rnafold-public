#
#
from __future__ import print_function
import sys
import time
import datetime
import gzip
import json
import logging
from math import floor
from socket import gethostname
from collections import Iterable, Set
from os import getpid
from cStringIO import StringIO
import redis
from binascii import crc32
from Bio import SeqIO
import config
import mysql_rnafold as db
from mysql_rnafold import Sequence2, SequenceSeries2, SequenceSeries2Updates, SequenceFloats2, SequenceIntegers2 # expose ORM objects through this namespace
from sqlalchemy import sql
from sqlalchemy.sql.expression import func
from sqlalchemy.exc import IntegrityError # expose errors through this namespace
from contextlib import contextmanager
import nucleic_compress


# Configuration
speciesCDSList = "species:taxid:%d:CDS"
speciesNameKey = "species:taxid:%d:name"
queueItemsKey = "queue:tag:awaiting-%s:members"
seqLengthKey = "CDS:taxid:%d:protid:%s:length-nt"
cdsSeqIdKey = "CDS:taxid:%d:protid:%s:seq-id"
cdsSeqChecksumKey = "CDS:taxid:%d:protid:%s:cds-seq-checksum"
shuffledSeqIdsKey = "CDS:taxid:%d:protid:%s:shuffled-seq-ids-v2"
workerKeepAliveKey = "status:worker-keep-alive"
jobStatusKey = "status:job:%s:progress"
taxGroupKey = "species:taxid:%d:tax-group"
speciesTranslationTableKey = "species:taxid:%d:genomic-transl-table"
speciesPropertyValueKey = "species:taxid:%d:properties:%s"
speciesPropertySourceKey = "species:taxid:%d:properties:%s:source"
speciesTaxIdKey = "species:name:%s:taxid"


# establish connections
# metadata server (redis)
r = redis.StrictRedis(host=config.host, port=config.port, db=config.db, password=config.password)
# sequences server (mysql)
session = db.Session()

"""
Enuemrate all CDSs from a given species
"""
def SpeciesCDSSource(taxId):
    # Make sure taxid is valid
    assert(r.exists(speciesNameKey % taxId))

    for protId in r.sscan_iter(speciesCDSList % taxId):
        yield protId


"""
Enumerate (and remove) queued items
"""
def QueueSource(queueTag):
    queueKey = queueItemsKey % queueTag
    while(True):
        itemToProcess = r.lpop(queueKey)
        if( not itemToProcess is None):
            yield itemToProcess
        else:
            return

"""
Yield all existing taxIds (in arbitrary order)
"""
def allSpeciesSource():
    keysFilter = speciesCDSList.replace("%d", "*")

    for key in r.keys(pattern=keysFilter):
        fields = key.split(":")
        yield int(fields[2])
        
"""
Return the number of CDS records for a given taxid
"""
def countSpeciesCDS(taxId):
    # Make sure taxid is valid
    assert(r.exists(speciesNameKey % taxId))

    return r.scard(speciesCDSList % taxId)

"""
Look-up the species name, by taxid (e.g., for display)
"""
def getSpeciesName(taxId):
    return r.get("species:taxid:%d:name" % taxId)

"""
"""
def getSpeciesFileName(taxId):
    p = r.get("species:taxid:%d:name" % taxId).split(" ")
    return p[0][0] + p[1]

"""
Get a taxonomic group containing this taxId (the groups were chosen, rather arbitrarily, for plotting)

TODO: Transition this to use NCBI data instead of the manually annotated values
"""
def getSpeciesTaxonomicGroup(taxId):
    return r.get(taxGroupKey % taxId)

"""
Get the nuclear translation table used for this genome
"""
def getSpeciesTranslationTable(taxId):
    return int(r.get(speciesTranslationTableKey % taxId))

"""
Return True if all specied taxIds are valid.
"""
def checkSpeciesExist(taxIds):
    if isinstance(taxIds, Iterable):  # usually, species will be a sequence of numeric taxid values
        if isinstance(taxIds, basestring):
            raise Exception("species cannot be string")
        # all set - proceed...
    else:
        taxIds = (taxIds,) # assume we got a single (numeric) taxid value
        
    for t in taxIds:
        if not (r.exists(speciesNameKey % t)):
            raise Exception("Species not found matching taxid=%s" % t)
    
    return True


def decompressSeriesRecord(compressed):
    if( compressed[:4] != "\x1f\x8b\x08\x00" ):
        raise Exception("Compressed format not recognized")
        
    f = gzip.GzipFile("", "rb", 9, StringIO(compressed))
    decoded = f.read()
    f.close()
    return decoded

def decodeJsonSeriesRecord(jsonstr):
    if jsonstr is None:
        return None
    else:
        return json.loads(jsonstr.replace('id=', '"id":').replace('seq-crc=', '"seq-crc":').replace('MFE-profile=','"MFE-profile":').replace('MeanMFE=','"Mean-MFE":'))


def decompressNucleicSequence(compressed):
    return nucleic_compress.decode(compressed)

def version():
    return "1.0"

"""
Wrapper for common CDS-related operations
"""
class CDSHelper(object):
    def __init__(self, taxId, protId):
        self._taxId = taxId
        self._protId = protId
        self._cache = {}
        self.updatescount = 0

    def _getScalarRedisProperty(self, cacheTag, redisKey, convertFunc = None):
        cachedVal = self._cache.get(cacheTag)
        if( cachedVal != None ):
            return cachedVal

        newVal = r.get( redisKey )
        if( convertFunc != None ):
            newVal = convertFunc(newVal)

        self._cache[cacheTag] = newVal
        return newVal

    def _getListRedisProperty(self, cacheTag, redisKey, convertFunc=None):
        cachedVal = self._cache.get(cacheTag)
        if( cachedVal != None ):
            return cachedVal

        newVal = r.lrange( redisKey, 0, -1 )
        if( convertFunc != None ):
            newVal = convertFunc(newVal)

        self._cache[cacheTag] = newVal
        return newVal


    def _getListRedisPropertyItem(self, redisKey, index, convertFunc=None):
        # No cache for this function

        newVal = r.lindex( redisKey, index )
        if( newVal is None ):
            raise Exception("CDS %s does not have an item in position %d (key=%s)" % (self._protId, index, redisKey))
        if( convertFunc != None ):
            newVal = convertFunc(newVal)

        return newVal


    def _fetchSequence(self, seqIds):
        if( not isinstance(seqIds, Iterable)):
            seqIds = (seqIds,)

        logging.info("Fetching sequences: %s" % seqIds)
        ret = {}
        notFoundInCache = set()
        # Try finding the sequence(s) in the cache
        for sid in seqIds:
            cacheTag = "%d:seq" % sid
            cachedVal = self._cache.get(cacheTag)
            if( cachedVal != None ):
                ret[sid] = cachedVal
            else:
                notFoundInCache.add(sid)

        # Is there something we did not find in the cache?
        if( notFoundInCache ):

            #raise Exception("Debug: fetching %s" % notFoundInCache)

            logging.info("Sequences not found in cache: %s" % notFoundInCache)

            records = None
            with db.connection.begin() as xact:
            
                # Get the sequence for this entry
                results = db.connection.execute( sql.select( (db.sequences2.c.id, db.sequences2.c.sequence)).select_from(db.sequences2).where(
                                db.sequences2.c.id.in_(notFoundInCache)
                                ) )  # Note: order_by not needed, because id is used
                if( results.rowcount < len(notFoundInCache) ):
                    logging.warning("Some results were not found for taxid=%d, protid=%s." % (self._taxId, self._protId))
                records = results.fetchall()
                del results

            # Decode all results, and store in cache
            for sid, encoded in records:
                assert(sid in notFoundInCache)
                logging.info("%d %d" % (sid, len(encoded)))
                logging.info(encoded[:10])
                    
                seq = nucleic_compress.decode(encoded)
                del encoded

                ret[sid] = seq

                # Store the sequence in the cache
                cacheTag = "%d:seq" % sid
                self._cache[cacheTag] = seq
                del seq
            del records

        assert(len(ret) <= len(seqIds))
        
        if len(ret)==1:
            return ret.popitem()[1]  # return the only item
        else:
            return ret


    def length(self):
        return self._getScalarRedisProperty( "cds-length-nt", seqLengthKey % (self._taxId, self._protId), int)
        
    def crc(self):
        return self._getScalarRedisProperty( "cds-crc", cdsSeqChecksumKey % (self._taxId, self._protId), int)

    def seqId(self):
        return self._getScalarRedisProperty( "cds-seq-id", cdsSeqIdKey % (self._taxId, self._protId), int)
        
    def shuffledSeqIds(self):
        return self._getListRedisProperty( "shuffled-seq-ids", shuffledSeqIdsKey % (self._taxId, self._protId), lambda x: map(int, x) )

    def getProtId(self):
        return self._protId

    def getTaxId(self):
        return self._taxId

    def getTranslationTable(self):
        return getSpeciesTranslationTable(self._taxId)

    def sequence(self):
        # Get the sequence-id for the CDS sequence
        seqId = self.seqId()
        # Fetch that sequence from the DB
        return self._fetchSequence(seqId)

    def getShuffledSeq(self, pos):
        assert(pos>=0)
        shuffleSeqId = self.getShuffledSeqId(pos)
        
        return self._fetchSequence(shuffleSeqId)

    def getShuffledSeq2(self, seqId):
        return self._fetchSequence(seqId)

    def getShuffledSeqId(self, pos):
        assert(pos>=0)
        # TODO - Cache this...
        return self._getListRedisPropertyItem( shuffledSeqIdsKey % (self._taxId, self._protId), pos, int )


    def clearCalculationResult(self, calculationId, shuffleId=-1):
        # TODO - impl this...
        pass

    def saveCalculationResult(self, calculationId, results, shuffleId=-1):
        seqId = None
        if( shuffleId < 0 ):
            seqId = self.seqId()
        else:
            seqId = self.getShuffledSeqId( shuffleId )

        gzBuffer = StringIO()
        f = gzip.GzipFile("", "wb", 9, gzBuffer)
        f.write(results)
        f.close()
        print(len(results), len(gzBuffer.getvalue()))

        # TODO - Add local session (instead of reusing the global one)
        
        ss2 = db.SequenceSeries2( sequence_id=seqId,
                                  content=gzBuffer.getvalue(),
                                  source=calculationId,
                                  ext_index=0)

        session.add(ss2)

        try:
            session.commit()
        except IntegrityError as e:
            print(e)
            # Ignore and continue...
            # TODO - improve this...
        gzBuffer.close()

        #r.hset( key, field, value )


    def saveCalculationResult2(self, calculationId, results, seqId, commit=True):

        gzBuffer = StringIO()
        f = gzip.GzipFile("", "wb", 9, gzBuffer)
        f.write(results)
        f.close()
        #print(len(results), len(gzBuffer.getvalue()))

        #print(repr(gzBuffer.getvalue())[:100])

        ss2 = db.SequenceSeries2Updates( sequence_id=seqId,
                                         content=gzBuffer.getvalue(),
                                         source=calculationId,
                                         ext_index=0)

        session.add(ss2)
        if commit:
            session.commit()
        else:
            self.updatescount += 1
            
        gzBuffer.close()

    def commitChanges(self):
        print(self.updatescount)
        self.updatescount = 0
        session.commit()

    def getDebugInfo(self):
        out = []
        out.append("CDS id: %d" % self.seqId())
        shuffids = self.shuffledSeqIds()
        out.append("Shuffle-ids: %s" % ([(i,v) for i,v in enumerate(shuffids)]))
        return out

    def isCalculationDone(self, calculationId, shuffleId=-1):
        seqId = None
        if( shuffleId < 0 ):
            seqId = self.seqId()
        else:
            seqId = self.getShuffledSeqId( shuffleId )

        count = db.connection.execute( sql.select(( sql.func.count('*'),)).select_from(db.sequence_series2).where(
                sql.and_(
                    db.sequence_series2.c.sequence_id==seqId,
                    db.sequence_series2.c.source==calculationId,
                    )) ).scalar()
        return count>0

    def getCalculationResult(self, calculationId, shuffleId=-1):
        seqId = None
        if( shuffleId < 0 ):
            seqId = self.seqId()
        else:
            seqId = self.getShuffledSeqId( shuffleId )

        compressed = db.connection.execute( sql.select(( db.sequence_series2.c.content,)).select_from(db.sequence_series2).where(
                sql.and_(
                    db.sequence_series2.c.sequence_id==seqId,
                    db.sequence_series2.c.source==calculationId,
                    )) ).scalar()

        try:
            return decompressSeriesRecord(compressed)
        except Exception as e:
            raise Exception("gzip format not recognized for sequence %d (taxId=%d, protId=%s)"  % (seqId, self._taxId, self._protId))

    """
    Return existing calculation results, for the specified calculationId, belonging to all specified shuffleIds
    Output: dict with results, indexed by calculation-Id.
            Missing results will have value None.

            Results will be uncompressed. In addition, if parseAsJson==True, results will be decoded as JSON.
    """
    def getCalculationResult2(self, calculationId, shuffleIds, parseAsJson=False):
        assert(shuffleIds >= -1)
        cdsSeqId = self.seqId()
        allSeqIds = self.shuffledSeqIds()

        print("requested calcid is: ", calculationId)
        print("requested shuffleIds: ", shuffleIds)

        def shuffleIdToSeqId(i):
            if i<0:
                return cdsSeqId
            elif i >= len(allSeqIds):
                return None
            else:
                return allSeqIds[i]

        def seqIdToShuffleId(i):
            if i==cdsSeqId:
                return -1
            else:
                return allSeqIds.index(i)

        # Create a list of sequence-ids for the requested shuffle-ids
        #requestedSeqIds = filter( lambda x: not x is None, map(shuffleIdToSeqId, shuffleIds) )
        requestedSeqIds = [y for y in [ shuffleIdToSeqId(x) for x in shuffleIds ] if not y is None]
        print("requested SeqIds (first 20): ", requestedSeqIds[:20])

        # Get the all result records
        results = db.connection.execute( sql.select(( db.sequence_series2.c.content, db.sequence_series2.c.sequence_id)).select_from(db.sequence_series2).where(
                sql.and_(
                    db.sequence_series2.c.source==calculationId,
                    db.sequence_series2.c.sequence_id.in_(requestedSeqIds)
                    )).order_by(db.sequence_series2.c.sequence_id) )
        if( results.rowcount < len(shuffleIds) ):
            print("Warning: some results were not found for taxid=%d, protid=%s." % (self._taxId, self._protId))
        records = results.fetchall()

        print("Got %d results from db" % (len(records)))

        out = dict( zip( shuffleIds, [None]*len(requestedSeqIds) ) )

        for record in records:
            compressed = record[0]
            currSeqId = record[1]

            decoded = None
            try:
                decoded = decompressSeriesRecord(compressed)
            except Exception as e:
                raise Exception("gzip format not recognized for sequence %d (taxId=%d, protId=%s)"  % (currSeqId, self._taxId, self._protId))

            out[ seqIdToShuffleId(currSeqId) ] = decoded

        if( parseAsJson ):
            out2 = {}
            for shuffleId, encoded in out.items():
                if( encoded is None ):
                    #print("Warning: Missing data for protein %s, shuffle-id %d" % (protId, shuffleId))
                    decoded = None
                else:
                    # Cover for silly JSON formatting error in an early version
                    decoded = decodeJsonSeriesRecord(encoded)
                    del encoded
                out2[shuffleId] = decoded
            out = out2

        assert(len(out) == len(shuffleIds))
        return out


    """
    Get lists of CDS positions, for each requested shuffleId, that are missing calculated MFE values.

    Params:
    * calculationId
    * shuffleIds -- list of shuffle-ids (contained in the range 0..n-1). the arbitrary list order determines the output order.
    * windowsToCheck -- a collection of window offsets, relative to CDS start

    Return value:
      returned value is a list with items matching those requested in shuffleIds.

      each item in the list can be:
    * if some positions are missing, the item will be a collection of CDS positions found in windowsToCheck
      but for which the results are missing.
    * if no values are missing, the item will be an empty list ([])
    * if the shuffle-id does not exist (in the results table), item will be None
    """
    def checkCalculationResultWithWindows(self, calculationId, shuffleIds, windowsToCheck):

        # Fetch the sequence-ids for the native and existing shuffles
        cdsSeqId = self.seqId()
        allSeqIds = self.shuffledSeqIds()

        #print("found %d shuffles" % len(allSeqIds))

        if not isinstance(windowsToCheck, Set):
            windowsToCheck = frozenset(windowsToCheck)

        # helper function to get the sequence-id (DB identifier) for each shuffle-id
        # note: existing shuffled seqs are assigned consecutive shuffle-ids (0..n-1)
        def shuffleIdToSeqId(i):
            if i<0:
                return cdsSeqId
            elif i<len(allSeqIds):
                return allSeqIds[i]
            else:
                return None

        #print("shuffleIds= %s" % shuffleIds)
        requestedSeqIds = filter( lambda x: not x is None, map(shuffleIdToSeqId, shuffleIds))

        # get the existing results for all sequences
        results = db.connection.execute( sql.select(( db.sequence_series2.c.sequence_id, db.sequence_series2.c.content)).select_from(db.sequence_series2).where(
                sql.and_(
                    db.sequence_series2.c.source==calculationId,
                    db.sequence_series2.c.sequence_id.in_(requestedSeqIds)
                    )) )
        records = results.fetchall()

        # check all sequence records to make sure they exist
        results = db.connection.execute( sql.select(( db.sequences2.c.id, )).select_from(db.sequences2).where(
                sql.and_(
                    db.sequences2.c.id.in_(requestedSeqIds),
                    db.sequences2.c.alphabet == db.Alphabets.RNA_Huff,
                    db.sequences2.c.source.in_( (db.Sources.ShuffleCDSv2_matlab, db.Sources.ShuffleCDSv2_python) ),
                    db.sequences2.c.sequence.isnot(None)
                    )) )
        allExistingSequenceIds = frozenset([x[0] for x in results.fetchall()])
        assert(len(allExistingSequenceIds) <= len(requestedSeqIds))

        # returned value is an array, with items matching those requested in shuffleIds (which are not required to be in consecutive or in ascending order)
        out = [None]*len(shuffleIds)

        # map each sequence-id to its position in shuffleIds
        sequenceIdToPos = dict(zip(map(shuffleIdToSeqId, shuffleIds), range(len(shuffleIds))))

        # iterate over each returned result, i.e., the existing results for all requested shuffle-ids in this gene
        for n, record in enumerate(records):
            #logging.warning("%d %s" % (n, record))
            currSeqId = record[0]
            compressed = record[1]

            # results are stored as gzip-compressed json.
            # first, decompress the result
            encoded = None
            try:
                encoded = decompressSeriesRecord(compressed)
            except Exception as e:
                raise Exception("gzip format not recognized for sequence %d (taxId=%d, protId=%s)"  % (currSeqId, self._taxId, self._protId))

            # next, parse the json string
            decoded = decodeJsonSeriesRecord(encoded)
            del encoded

            # extract the indices fo which profiles were already calculated
            # (the indices also correspond to positions in CDS coordinates; missing values are stored as json nulls)
            calculatedWindows = frozenset([ n for n,v in enumerate(decoded['MFE-profile']) if not v is None ])
            #print(sorted(list(calculatedWindows)))

            # calculate the missing windows (using set difference)
            missingWindows = windowsToCheck - calculatedWindows

            #print("Existing: %d  Needed: %d   Missing: %d" % (len(calculatedWindows), len(windowsToCheck), len(missingWindows)))
            #if( len(calculatedWindows)>151 ):
            #    print("Missing: %s" % str(missingWindows))

            #print(currSeqId)
            #print(sequenceIdToPos[currSeqId])
            
            if( missingWindows ):
                out[ sequenceIdToPos[currSeqId] ] = missingWindows
            else:
                out[ sequenceIdToPos[currSeqId] ] = []
            
            #print(out[ sequenceIdToPos[currSeqId] ])
            

        print("-----------")
        #logging.warning("---------------------")
        # For items where the output value (matching item in out) is None, the are two cases:
        # 1) The random sequence might not exist; in this case, we will return None, and the sequence will need to be created
        # 2) The random seqeunce might already exists; in this case, we will return a list of all windows, since they are all missing and ready to be computed
        itemPositionsMissingResults = [i for i,x in enumerate(out) if x is None]
        #logging.warning("shuffleIds: %s" % shuffleIds)
        #logging.warning("itemPositions: %s" % itemPositionsMissingResults)
        #logging.warning("allSeqIds: %s" % allSeqIds)
        
        for pos in itemPositionsMissingResults:
            #logging.warning("pos: %d" % pos)
            itemShuffleId = shuffleIds[pos]
            if itemShuffleId == -1:
                logging.warning("Skipping checking native CDS sequence")
                out[pos] = frozenset(windowsToCheck)
                continue   # we assume the native CDS sequence must already exist...
            
            #logging.warning("itemShuffleId: %d" % itemShuffleId)
            if( itemShuffleId < len(allSeqIds) ):
                #logging.warning("*********(-1)************")
                shuffleSequenceId = allSeqIds[itemShuffleId]
                #logging.warning("shuffleSequenceId: %d" % shuffleSequenceId)
            else:
                #logging.warning("*********(0)************")
                continue # case 1) - nothing else needs to be done

            #logging.warning("*********(1)************")
            if( shuffleSequenceId in allExistingSequenceIds ): # the randomized sequence for this shuffle-id already exists (although we found no results for it)
                #logging.warning("*********(2)************")
                out[pos] = frozenset(windowsToCheck) # case 2) - re-calculate all windows for this sequence
            else:
                pass # case 1) - nothing else needs to be done

        assert(len(out) == len(shuffleIds))
        return out

    
    """
    Enqueue the following items, as a single work item.
    shuffleId may be a single shuffleId (in the range -1..(n-1)), or a sequence of shuffleIds
    """
    def enqueueForProcessing(self, computationTag, shuffleId=-1, lastWindowStart=None, windowStep=None):
        queueKey = queueItemsKey % computationTag

        if( not isinstance(shuffleId, Iterable)):
            shuffleId = (shuffleId,)

        if( not shuffleId): # No shuffle-ids to enqueue
            print("Warning: enqueueForProcessing() called with no shuffle-ids to add...")
            return

        seqIds = []
        for sid in shuffleId:
            if( sid < 0 ):
                seqIds.append( self.seqId() )
            else:
                seqIds.append( self.getShuffledSeqId( sid ) )

        optionalParams = ""
        if( lastWindowStart is not None and windowStep is not None):
            optionalParams = "/%d/%d" % (lastWindowStart, windowStep)

        queueItem = "%d/%s/%s/%s%s" % (self._taxId, self._protId, ",".join(map(str, seqIds)), ",".join(map(str, shuffleId)), optionalParams) # Separator changed to '/' as a work-around for protein-ids containing ':'...
        #print(queueItem)
        
        r.rpush(queueKey, queueItem)

    def dropShuffledSeqs(self, lastItemToKeep=0):
        shuffledSeqIds = self.shuffledSeqIds()

        cdsRecordsBeforeDeleting = len(shuffledSeqIds)
        shuffledSeqIdsToDelete = shuffledSeqIds[lastItemToKeep:None]  # delete the range lastItemToKeep:end

        expectedDeletedRecords = len(shuffledSeqIdsToDelete)
        
        if not expectedDeletedRecords:
            return 0

        print("Debug: Found %d records in redis; expecting to delete %d records from mysql" % (cdsRecordsBeforeDeleting, expectedDeletedRecords))

        result = db.connection.execute( db.sequences2.delete().where(
                    db.sequences2.c.id.in_(shuffledSeqIdsToDelete)
                    ) )

        print("Debug: deleted %d records from mysql" % result.rowcount)

        assert(result.rowcount <= expectedDeletedRecords )  # can't delete more records than we specified!
        
        if(result.rowcount < expectedDeletedRecords):
            print("Warning: Deleted less records (%d) than expected (%d)" % (result.rowcount, expectedDeletedRecords))
            
        shuffleCountBeforeDeleting = r.llen( shuffledSeqIdsKey % (self._taxId, self._protId) )
        print("Debug: Before deleting from redis: %d items" % shuffleCountBeforeDeleting )
        #r.delete(shuffledSeqIdsKey % (self._taxId, self._protId))
        r.ltrim( shuffledSeqIdsKey % (self._taxId, self._protId), 0, lastItemToKeep-1 )

        shuffleCountAfterDeleting = r.llen( shuffledSeqIdsKey % (self._taxId, self._protId) )
        print("Debug: After deleting from redis: %d items" % shuffleCountAfterDeleting )

        assert( shuffleCountAfterDeleting <= lastItemToKeep )

        return result.rowcount

    
    def dropNativeSeq(self):
        nativeSeqId = self.seqId()

        if( nativeSeqId is None ):
            raise Exception("No native seq defined")

        result = db.connection.execute( db.sequences2.delete().where(
                    db.sequences2.c.id == nativeSeqId
                    ) )
        
        if(result.rowcount >= 1):
            r.delete( cdsSeqIdKey % (self._taxId, self._protId) )
        else:
            raise Exception("Failed to delete native seq with id %d (prot-id=%s)" % (nativeSeqId, self._protId))

    def dropRecord(self, dropCDSKeys=False):
        keysToDelete = set()
        
        # Collect the list of all values for this protein
        # Note: This is *really* slow
        #       As an alternative, use matchCDSKeyNamesSource (see drop_species.py)
        if( dropCDSKeys ):
            for k in r.scan_iter(match='CDS:taxid:%d:protid:%s:*'%(self._taxId, self._protId) ):
                keysToDelete.add(k)
            # Delete all keys for this protein
            for k in keysToDelete:
                r.delete(k)

        # Delete the protId entry from the proteins set
        r.srem(speciesCDSList % self._taxId, self._protId)


def matchCDSKeyNamesSource(taxId):
    for k in r.scan_iter(match='CDS:taxid:%d:protid:*' % taxId ):
        yield k




"""
Get the number of items queued
"""
def numItemsInQueue(computationTag):
    queueKey = queueItemsKey % computationTag
    return r.llen(queueKey)


"""
Update worker's heartbeat time; optionally, update job's progress counter
"""
SecsPerDay = float(24*60*60)
def updateJobStatus_ItemCompleted( workerKey, itemCompleted=None ):
    # Send job keep-alive
    timeSinceEpoch = datetime.datetime.utcnow() - datetime.datetime(2016, 06, 01)
    timeIndex = timeSinceEpoch.days + float(timeSinceEpoch.seconds)/SecsPerDay

    r.zadd( workerKeepAliveKey, timeIndex, workerKey)

    if( not itemCompleted is None):
        # Update job progress counter
        r.incr( jobStatusKey % itemCompleted )

"""
Create a string (redis key compatible) to identify a running worker process
"""
def createWorkerKey(computationTag):
    return "%s-%d-%s" % (gethostname(), getpid(), computationTag)

"""
Calculate CRC for a given sequence
"""
def calcCrc(seq):
    return crc32(str(seq).lower()) & 0xffffffff


def splitLongSequenceIdentifier(longIdentifier):
    recordIdentifier = ()
    if( longIdentifier.find('/') == -1 ):
        recordIdentifier = longIdentifier.split(":")
    else:
        recordIdentifier = longIdentifier.split("/")

    if( len(recordIdentifier) != 4 ):
        raise Exception("Invalid record identifier format: '%s'" % longIdentifier)
    

    return (int(recordIdentifier[0]), recordIdentifier[1], int(recordIdentifier[2]), int(recordIdentifier[3]))
    

"""
Get a merged list of matching records showing an update and the original record, i.e., a record from sequence_series2_updates, and the corresponding original record in sequence_series.
Updates are allowed for non-existent records (i.e., a new version was computed but the original recored is missing).

Note: There is no treatment of multiple updates for the same record (i.e., both will be returned).
"""
def seriesUpdatesSource(calculationId, bulkSize=500):

    lastId = 0
    offset = 0 

    while(True):
        # Approach 1: Use sequence-id (not primary-key!)
        #results = db.connection.execute( sql.select(( db.sequence_series2_updates.c.sequence_id, )).select_from(db.sequence_series2_updates).where(
        #    sql.and_(
        #        db.sequence_series2_updates.c.source==calculationId,
        #        db.sequence_series2_updates.c.sequence_id > lastId
        #    )).order_by(db.sequence_series2_updates.c.sequence_id).limit(bulkSize) ).fetchall()

        # Approach 2: Use native DB paging impl. (LIMIT / OFFSET)
        #results = db.connection.execute( sql.select(( db.sequence_series2_updates.c.sequence_id, )).select_from(db.sequence_series2_updates).where(
        #    sql.and_(
        #        db.sequence_series2_updates.c.source==calculationId,
        #    )).limit(bulkSize).offset(offset) ).fetchall()

        # Approach 3: Use primary-key (dummy_id)
        # Note: This seems to be the most efficient
        updates = db.connection.execute( sql.select(( db.sequence_series2_updates.c.dummy_id, db.sequence_series2_updates.c.sequence_id, db.sequence_series2_updates.c.content )).select_from(db.sequence_series2_updates).where(
            sql.and_(
                db.sequence_series2_updates.c.source==calculationId,
                db.sequence_series2_updates.c.dummy_id > lastId
            )).order_by(db.sequence_series2_updates.c.dummy_id).limit(bulkSize) ).fetchall()

        if( len(updates) == 0 ):
            return

        # Fetch the corresponding original records
        updatedIds = [x[1] for x in updates]

        originals = db.connection.execute( sql.select(( db.sequence_series2.c.sequence_id, db.sequence_series2.c.content)).select_from(db.sequence_series2).where(
            sql.and_(
                db.sequence_series2.c.source==calculationId,
                db.sequence_series2.c.sequence_id.in_(updatedIds)
            )) ).fetchall()
        
        lastId = updates[-1][0]
        #offset += len(updates)

        if(len(updates) != len(originals)):
            print("Warning: updates returned different number of results than original rows (updates=%d, originals=%d)" % (len(updates), len(originals)))
            assert(len(updates) >= len(originals)) # There cannot be more matching original records than updates for them

        out = []
        # Merge original and updated results on sequence_id (use left merge)
        
        originalsById = dict( zip( [x[0] for x in originals], originals ) )
        assert(len(originalsById)==len(originals))
        
        for updateRecord in updates:
            # Find the matching original record (if any)
            updatedId = updateRecord[1]
            matchingOriginal = originalsById.get(updatedId)  # May be None
            
            # Decode the update record
            # Decoded format: (sequence_id, content_record, dummy_id)
            convertedUpdateRecord = (updateRecord[1], decodeJsonSeriesRecord( decompressSeriesRecord( updateRecord[2] )), updateRecord[0] )

            # Decode the original record
            # Decoded format: (sequence_id, content_record)
            convertedOriginalRecord = None
            if( not matchingOriginal is None ):
                convertedOriginalRecord = (matchingOriginal[0], decodeJsonSeriesRecord( decompressSeriesRecord( matchingOriginal[1] )) )
                
            out.append( (convertedUpdateRecord, convertedOriginalRecord) )
        
        yield(out)



"""
Return a dict mapping sequence-ids to sequences (in compressed format), for all native CDS sequences belonging to a species
This is much faster than fetching the sequences separately.

If fraction and modulus are set, only return the sequences for which the sequenceId % modulus == fraction
(used to parallelize)
"""
def getAllNativeCDSsForSpecies(taxId, fraction=None, modulus=None):
    # First, collect all sequence-ids for the given species. This manual filtering by species is much faster than simply fetching all computation results...
    #print("Collecting sequence ids for taxid=%d..." % taxId)
    sequenceIdsForTaxid = set()
    for protId in SpeciesCDSSource(taxId):
        cds = CDSHelper(taxId, protId)
        sequenceIdsForTaxid.add(cds.seqId())

    if( not fraction is None ):
        assert(fraction >= 0)
        assert(modulus > 0)
        assert(fraction < modulus)
        sequenceIdsForTaxid = set( [x for x in sequenceIdsForTaxid if x % modulus == fraction] )
    
    print("Fetching results for %d sequences..." % len(sequenceIdsForTaxid))
    calculated = db.connection.execute( sql.select(( db.sequences2.c.id, db.sequences2.c.sequence)).select_from(db.sequences2).where(
                sql.and_(
                    db.sequences2.c.id.in_(sequenceIdsForTaxid)
                    )) ).fetchall()
    out = dict(calculated)

    print("Got %d results" % len(out))

    return out


def setSpeciesProperty(taxId, propName, propVal, source, overwrite=True):
    if( propName.find(":") != -1 ):
        raise Exception("Invalid property name '%s'" % propName)

    if (not overwrite) and r.exists(speciesPropertyValueKey  % (taxId, propName)):
        return
    
    r.set(speciesPropertyValueKey  % (taxId, propName), propVal)
    r.set(speciesPropertySourceKey % (taxId, propName), source)

def getSpeciesProperty(taxId, propName):
    propVal    = r.get(speciesPropertyValueKey % (taxId, propName))
    propSource = r.get(speciesPropertySourceKey % (taxId, propName))
    return (propVal, propSource)

def getSpeciesTemperatureInfo(taxId):
    tempCategoricalProp = getSpeciesProperty(taxId, 'temperature-range')
    tempNumericalProp   = getSpeciesProperty(taxId, 'optimum-temperature')

    syntheticCategory = None

    if tempCategoricalProp[0] is None:
        if not tempNumericalProp[0] is None:
            tempVal = float(tempNumericalProp[0])


            #Counter({'Mesophilic': 79, 'Hyperthermophilic': 20, 'Thermophilic': 13, 'Psychrophilic': 4, 'Unknown': 1})

            if tempVal < 25.0:
                assert(tempVal > -30.0)
                syntheticCategory = 'Psychrophilic'

            elif tempVal < 42.0:
                syntheticCategory = 'Mesophilic'

            elif tempVal < 70.0:
                syntheticCategory = 'Thermophilic'

            else:
                assert(tempVal < 120.0)
                syntheticCategory = 'Hyperthermophilic'

            tempCategoricalProp = (syntheticCategory, 'FromOptimumTemperature')
        else:
            tempCategoricalProp = ('Unknown', 'Unknown')

    return( tempNumericalProp, tempCategoricalProp) 
            

    

def getGenomicGCContent(taxId):
    (propVal, _)    = getSpeciesProperty(taxId, "gc-content") # annotated by ncbi_entrez.py
    
    if not propVal is None:
        return float(propVal)
    else:
        return None
    

def getAllComputedSeqsForSpecies(calculationIds, taxId, maxShuffleId=100):
    # First, collect all sequence-ids for the given species. This manual filtering by species is much faster than simply fetching all computation results...
    #print("Collecting sequence ids for taxid=%d..." % taxId)

    if type(calculationIds)==type(0):
        calculationIds = (calculationIds,)
        
    sequenceIdsForTaxid = set()
    for protId in SpeciesCDSSource(taxId):
        cds = CDSHelper(taxId, protId)

        newIds = set()
        newIds.add(cds.seqId())

        allShuffleIds = cds.shuffledSeqIds()
        if( len(allShuffleIds) > maxShuffleId ):  # limit the included shuffled ids to the first N=maxShuffleId
            newIds.update(cds.shuffledSeqIds()[:maxShuffleId] )
        else:
            newIds.update(cds.shuffledSeqIds() )
            
        sequenceIdsForTaxid |= newIds
    
    print("Fetching results for %d sequences..." % len(sequenceIdsForTaxid))
    calculated = db.connection.execute( sql.select(( db.sequence_series2.c.sequence_id, db.sequence_series2.c.content)).select_from(db.sequence_series2).where(
                sql.and_(
                    db.sequence_series2.c.source.in_(calculationIds),
                    db.sequence_series2.c.sequence_id.in_(sequenceIdsForTaxid)
                    )) ).fetchall()
    #print("Converting data...")
    #out = set([x[0] for x in calculated])
    out = dict(calculated)

    print("Got %d results" % len(out))

    return out

class DropSequenceWithResults(object):
    def __init__(self, bulkSizeTarget=1000):
        self._bulkSize = bulkSizeTarget
        self.reset()

    """
    Mark a single sequence for dropping.
    This includes:
    * deleting the sequence from sequences2
    * deleting related results from sequences_series2
    * updating the redis keys to remove the deleted sequences
    """
    def markSequenceForDropping(self, taxId, protId, sequenceId):
        self._sequenceIdsToBeDropped.append((taxId, protId, sequenceId))
        ownerId = (taxId, protId)
        
        if( ownerId in self._owners):
            self._owners[ownerId].append(sequenceId)
        else:
            self._owners[ownerId] = [sequenceId]

        if(len(self._sequenceIdsToBeDropped) >= self._bulkSize ):
            self.performDroppingNow()

    def reset(self):
        self._sequenceIdsToBeDropped = []
        self._owners = {}

    def performDroppingNow(self):
        if( not self._sequenceIdsToBeDropped ):
            return

        droppedSeqsFromRedis = 0
        droppedSeqsFromMysql = 0
        droppedSeqResults = 0
        
        # 1) Delete pointers to sequences being dropped from redis, by updating all owner sequences
        for owner, sequenceIdsToDelete in self._owners.items():
            taxId  = owner[0]
            protId = owner[1]
            keyId = shuffledSeqIdsKey % (taxId, protId)
            originalItems = list( map(int, r.lrange( keyId, 0, -1 ) ) )
            assert(len(originalItems) >= len(sequenceIdsToDelete))

            updatedItems = set(originalItems) - set(sequenceIdsToDelete)  # may not be 'stable' (i.e., may change order of remaining items)

            assert(len(updatedItems) < len(originalItems) )

            print("(debug) deleting key %s  -- %s" % (keyId, originalItems) )

            r.delete( keyId )
            r.rpush( keyId, *updatedItems )

            assert( r.llen(keyId) == len(updatedItems) )
            droppedSeqsFromRedis += len(sequenceIdsToDelete)
 

        # 2) Delete the sequence records
        seqIdsForDeletion = [x[2] for x in self._sequenceIdsToBeDropped]
        result = db.connection.execute( db.sequences2.delete().where(
            db.sequences2.c.id.in_( seqIdsForDeletion )
        ) )
        droppedSeqsFromMysql = result.rowcount
        #if(result.rowcount < len(seqIdsForDeletion)):
        #    raise Exception("Only deleted %d records" % result.rowcount)

        # 3) Delete calculation results for the sequences
        # Note: this will delete all calculations for these sequence-ids (i.e., calc-id is not used)
        result = db.connection.execute( db.sequence_series2.delete().where(
            db.sequence_series2.c.sequence_id.in_( seqIdsForDeletion )
        ) )
        droppedSeqResults = result.rowcount

        # 4) Clear the items-to-be-deleted container
        self.reset()

        print("INFO: Deletion counts: Redis %d, Mysql %d seqs, %d results" % (droppedSeqsFromRedis, droppedSeqsFromMysql, droppedSeqResults) )

@contextmanager
def session_scope():
    """
    Source: http://docs.sqlalchemy.org/en/latest/orm/session_basics.html
    """
    session = db.Session()
    try:
        yield session
        session.commit()
    except:
        session.rollback()
        raise
    finally:
        session.close()
    
    
        
class AddShuffledSequences(object):
    def __init__(self, taxId, protId, sourceTag):
        self._taxId = taxId
        self._protId = protId
        self._sourceTag = sourceTag
        
        assert( r.sismember(speciesCDSList % (taxId,),  protId))

    def addSequence(self, shuffleId, seq):
        # Store the shuffled CDS sequence

        logging.info('compress seq: %s...' % seq[:15])

        # Compress the CDS sequence
        encodedCds = nucleic_compress.encode(seq)

        newSequenceId = None
        with session_scope() as s:
            s1 = db.Sequence2(
                sequence=encodedCds,
                alphabet=db.Alphabets.RNA_Huff,
                source=self._sourceTag)

            s.add(s1)
            s.commit()    # perform commit now, so we can obtain the generated sequence id
            
            newSequenceId = s1.id
            
        assert(not newSequenceId is None)

        # mysql work done; record the new sequence id in redis
        shufflesKey = shuffledSeqIdsKey % (self._taxId, self._protId)

        with r.pipeline() as pipe: # will automatically call pipe.reset()
            attempt = 1
            while True:
                try:
                    pipe.watch( shufflesKey )

                    listItems = pipe.llen( shufflesKey )

                    # note: redis list indices are 0-based (e.g., for lset, lrange)
        
                    if( shuffleId < listItems ):   # we are replacing an existing item
                        pipe.lset( shufflesKey, shuffleId, newSequenceId )
                    elif( shuffleId == listItems ): # we are appending the next item
                        pipe.rpush( shufflesKey, newSequenceId )
                    else:
                        pipe.reset() # not required (will be called by context manager)
                        raise Exception("Can't append items past end of list (key=%s; length=%d; new index=%d)" % (shufflesKey, listItems, shuffleId))

                    pipe.execute()

                    break # All done!
                except redis.WatchError:
                    pass
                    
                attempt += 1
                if( attempt > 10 ):
                    raise Exception("Failed too many times while trying to update key %s" % shufflesKey)
                else:
                    continue

        testValue = r.lindex( shufflesKey, shuffleId )
        if( testValue is None):
            raise Exception("Updating sequence taxId=%d, protId=%s, shuffleId=%d seems to have failed to update the key (key=%s): Actual value is missing" % (self._taxId, self._protId, shuffleId, shufflesKey) )
            
        if( int(testValue) != newSequenceId ): # verify we successfully updated the key
            raise Exception("Updating sequence taxId=%d, protId=%s, shuffleId=%d seems to have failed to update the key (key=%s): Expected value=%d, actual value: %s" % (self._taxId, self._protId, shuffleId, shufflesKey, newSequenceId, testValue) )
        # Note: this is not part of the transaction and can theoretically fail if someone rewrites the key (or otherwise updates the list in a way they shouldn't)
        
        return newSequenceId
