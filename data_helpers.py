#
#
from __future__ import print_function
import sys
import time
import datetime
import gzip
from socket import gethostname
from os import getpid
from cStringIO import StringIO
#import re
import redis
from binascii import crc32
from Bio import SeqIO
import config
import mysql_rnafold as db
from sqlalchemy import sql
from sqlalchemy.sql.expression import func
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

def getSpeciesTaxonomicGroup(taxId):
    return r.get(taxGroupKey % taxId)


"""
Wrapper for common CDS-related operations
"""
class CDSHelper(object):
    def __init__(self, taxId, protId):
        self._taxId = taxId
        self._protId = protId
        self._cache = {}

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


    def _fetchSequence(self, seqId):
        cacheTag = "%d:seq" % seqId
        cachedVal = self._cache.get(cacheTag)
        if( cachedVal != None ):
            return cachedVal

        # Get the sequence for this entry
        encoded = db.connection.execute( sql.select( (db.sequences2.c.sequence,)).select_from(db.sequences2).where(
                        db.sequences2.c.id==seqId )
                        ).scalar()
        seq = nucleic_compress.decode(encoded)

        self._cache[cacheTag] = seq
        return seq
        

    def length(self):
        return self._getScalarRedisProperty( "cds-length-nt", seqLengthKey % (self._taxId, self._protId), int)
        
    def crc(self):
        return self._getScalarRedisProperty( "cds-crc", cdsSeqChecksumKey % (self._taxId, self._protId), int)

    def seqId(self):
        return self._getScalarRedisProperty( "cds-seq-id", cdsSeqIdKey % (self._taxId, self._protId), int)
        
    def shuffledSeqIds(self):
        return self._getListRedisProperty( "shuffled-seq-ids", shuffledSeqIdsKey % (self._taxId, self._protId), lambda x: map(int, x) )

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

        ss2 = db.SequenceSeries2( sequence_id=seqId,
                                  content=gzBuffer.getvalue(),
                                  source=calculationId,
                                  ext_index=0)

        session.add(ss2)

        try:
            session.commit()
        except sqlalchemy.exc.IntegrityError as e:
            print(e)
            # Ignore and continue...
            # TODO - improve this...
        gzBuffer.close()

        #r.hset( key, field, value )


    def saveCalculationResult2(self, calculationId, results, seqId):

        gzBuffer = StringIO()
        f = gzip.GzipFile("", "wb", 9, gzBuffer)
        f.write(results)
        f.close()
        print(len(results), len(gzBuffer.getvalue()))

        print(repr(gzBuffer.getvalue())[:100])

        ss2 = db.SequenceSeries2( sequence_id=seqId,
                                  content=gzBuffer.getvalue(),
                                  source=calculationId,
                                  ext_index=0)

        session.add(ss2)
        session.commit()
        gzBuffer.close()

        #r.hset( key, field, value )


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

        if( compressed[:4] != "\x1f\x8b\x08\x00" ):
            print("Warning: gzip format not recognized for sequence %d (taxId=%d, protId=%s)"  % (seqId, self._taxId, self._protId))

        f = gzip.GzipFile("", "rb", 9, StringIO(compressed))
        decoded = f.read()
        f.close()
        return decoded


    def getCalculationResult2(self, calculationId, shuffleIds):

        cdsSeqId = self.seqId()
        allSeqIds = self.shuffledSeqIds()

        def shuffleIdToSeqId(i):
            if i<0:
                return cdsSeqId
            else:
                return allSeqIds[i]
            
        requestedSeqIds = map(shuffleIdToSeqId, shuffleIds)

        results = db.connection.execute( sql.select(( db.sequence_series2.c.content, db.sequence_series2.c.sequence_id)).select_from(db.sequence_series2).where(
                sql.and_(
                    db.sequence_series2.c.source==calculationId,
                    db.sequence_series2.c.sequence_id.in_(requestedSeqIds)
                    )) )
        if( results.rowcount < len(shuffleIds) ):
            print("Warning: some results were not found for taxid=%d, protid=%s." % (self._taxId, self._protId))
            #return None

        if( (results.rowcount < len(shuffleIds) - 5) or (results.rowcount < 10) ):
            print("Warning: Too many results were missing for taxid=%d, protid=%s." % (self._taxId, self._protId))
            return None

        records = results.fetchall()

        out = [None]*len(shuffleIds)

        sequenceIdToPos = dict( zip( requestedSeqIds, range(len(shuffleIds)) ) )
        
        for n, record in enumerate(records):
            compressed = record[0]
            currSeqId = record[1]

            if( compressed[:4] != "\x1f\x8b\x08\x00" ):
                print("Warning: gzip format not recognized for sequence %d (taxId=%d, protId=%s)"  % (currSeqId, self._taxId, self._protId))

            f = gzip.GzipFile("", "rb", 9, StringIO(compressed))
            decoded = f.read()
            f.close()

            out[ sequenceIdToPos[currSeqId] ] = decoded

        return out


    def checkCalculationResult(self, calculationId, shuffleIds):

        cdsSeqId = self.seqId()
        allSeqIds = self.shuffledSeqIds()

        def shuffleIdToSeqId(i):
            if i<0:
                return cdsSeqId
            else:
                return allSeqIds[i]
            
        requestedSeqIds = map(shuffleIdToSeqId, shuffleIds)

        results = db.connection.execute( sql.select(( db.sequence_series2.c.sequence_id,)).select_from(db.sequence_series2).where(
                sql.and_(
                    db.sequence_series2.c.source==calculationId,
                    db.sequence_series2.c.sequence_id.in_(requestedSeqIds)
                    )) )
        records = results.fetchall()

        out = [False]*len(shuffleIds)

        sequenceIdToPos = dict( zip( requestedSeqIds, range(len(shuffleIds)) ) )
        
        for n, record in enumerate(records):
            currSeqId = record[0]
            out[ sequenceIdToPos[currSeqId] ] = True

        return out
            
            

    def enqueueForProcessing(self, computationTag, shuffleId=-1):
        queueKey = queueItemsKey % computationTag
        seqId = None
        if( shuffleId < 0 ):
            seqId = self.seqId()
        else:
            seqId = self.getShuffledSeqId( shuffleId )

        queueItem = "%d:%s:%d:%d" % (self._taxId, self._protId, seqId, shuffleId)
        #print(queueItem)
        
        r.rpush(queueKey, queueItem)

    def dropShuffledSeqs(self):
        shuffledSeqIds = self.shuffledSeqIds()

        #print(shuffledSeqIds[:10])

        result = db.connection.execute( db.sequences2.delete().where(
                    db.sequences2.c.id.in_(shuffledSeqIds)
                    ) )
        
        if(result.rowcount == len(shuffledSeqIds)):
            r.delete(shuffledSeqIdsKey % (self._taxId, self._protId))
            return result.rowcount
        else:
            raise Exception("Only deleted %d records" % result.rowcount)
        
    def dropNativeSeq(self):
        nativeSeqId = self.seqId()

        result = db.connection.execute( db.sequences2.delete().where(
                    db.sequences2.c.id == nativeSeqId
                    ) )
        
        if(result.rowcount >= 1):
            r.delete( cdsSeqIdKey % (self._taxId, self._protId) )
        else:
            raise Exception("Failed to delete native seq with id %d (prot-id=%s)" % (nativeSeqId, self._protId))

    def dropRecord(self):
        keysToDelete = set()
        # Collect the list of all values for this protein
        for k in r.scan_iter(match='CDS:taxid:%d:protid:%s:*'%(self._taxId, self._protId) ):
            keysToDelete.add(k)
        # Delete all keys for this protein
        for k in keysToDelete:
            r.delete(k)
        # Delete the protId entry from the proteins set
        r.srem(speciesCDSList % self._taxId, self._protId)

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


def getAllComputedSeqs(calculationId):
    calculated = db.connection.execute( sql.select(( db.sequence_series2.c.sequence_id,)).select_from(db.sequence_series2).where(
                sql.and_(
                    db.sequence_series2.c.source==calculationId,
                    )) ).fetchall()
    out = set()
    for x in calculated:
        out.add(x[0])

    return out

    

