# Upload sliced archive files to glacier using multipart upload
# Example: python glacier_multipart_upload.py --parts "test_archive.part_*" --account-id 123456789 --vault-name my-vault-name --archive-description "Test"
# Reference:
# https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/glacier.html
#
import os
import hashlib
import codecs
import logging
import argparse
import boto3


# Configuration
logging.basicConfig(filename='glacier_multipart_upload.log', level=logging.DEBUG, format='%(asctime)s %(module)s %(message)s')

bufferSize       =   32*1024
treeLeafHashSize = 1024*1024 # 1MB - Reference: https://docs.aws.amazon.com/amazonglacier/latest/dev/checksum-calculations.html
assert( treeLeafHashSize % bufferSize == 0 )

logging.info("Buffer size: {}".format(bufferSize))
logging.info("treeLeafHashSize: {}".format(treeLeafHashSize))

def findParts(args):
    from glob import glob
    return sorted([x for x in glob(args.parts) if os.path.exists(x)])

def getPartSizes(partFiles):
    return [os.stat(p).st_size for p in partFiles]
"""
Determine the part size of the multipart-archive (i.e. the size of all parts except the last)
"""
def findPartSize(sizes):
    assert(all([s==sizes[0] for s in sizes[:-1]]))  # all sizes except the last size must be the same
    assert(sizes[-1] <= sizes[0]) # the last size can be smaller
    return sizes[0]

"""
Calculate sha256 hash for file stored on disk
"""
def calcHashForFile(partFile):
    m = hashlib.sha256()
    with open( partFile, "rb" ) as f:
        while True:
            chunk = f.read( bufferSize )
            if len( chunk ) > 0:
                m.update( chunk )
                
            else:  # EOF
                break
            
    return m.digest()

def getPairs(items):
    i = iter(items)

    try:
        while True:
            c = 0
            a0 = next(i)
            c = 1
            a1 = next(i)
            c = 2
            yield (a0,a1)
    except StopIteration:
        if c == 1:
            yield (a0,)

#print(list(getPairs(range(13))))
#print(list(getPairs(range(12))))
#print(list(getPairs(range(11))))
#print(list(getPairs(range(2))))
#print(list(getPairs(range(1))))
#print(list(getPairs([])))
#import sys
#sys.exit()
    

# Join hashes:
def joinHashesAsTree(hashes):
    
    def joinHashesSingleLevel(hashes):
        out = []
        for pair in getPairs(hashes):
            if len(pair)==1:
                out.append(pair[0])

            elif len(pair)==2:
                m = hashlib.sha256()
                m.update(pair[0])
                m.update(pair[1])

                out.append( m.digest() )

            else:
                assert(False)

        return out

    while len(hashes) > 1:
        hashes = joinHashesSingleLevel(hashes)
        
    return hashes[0]

"""
Calculate tree hash for file stored on disk
Reference: https://docs.aws.amazon.com/amazonglacier/latest/dev/checksum-calculations.html
"""
def calcTreeHashForFile(partFile, treeLeafHashSize=1024*1024):
    
    m = hashlib.sha256()
    leafHashes = []
    with open( partFile, "rb" ) as f:
        sizeSoFar = 0
        while True:
            chunk = f.read( bufferSize )
            if len( chunk ) > 0:
                m.update( chunk )
                sizeSoFar += len(chunk)
                assert( sizeSoFar <= treeLeafHashSize )

                if sizeSoFar == treeLeafHashSize:
                    leafHashes.append( m.digest() )
                    m = hashlib.sha256() # reset
                    #print(sizeSoFar)
                    sizeSoFar = 0
                
            else:  # EOF
                break

    if sizeSoFar > 0:
        leafHashes.append( m.digest() )
        m = hashlib.sha256() # reset
        #print(sizeSoFar)
        sizeSoFar = 0

        

    return joinHashesAsTree(leafHashes)


def calcAllHashes(partFiles):
    out = {}
    partNum = 1
    totalParts = len(partFiles)
    for part in partFiles:
        logging.warning("Calculating hash for part {}/{}...".format(partNum, totalParts))
        partHash = calcTreeHashForFile( part )
        out[part] = partHash
        partNum += 1

    return out

def up(args):

    partFiles = findParts(args)
    #partFiles = partFiles[:2] + [partFiles[-1]]  # DEBUG ONLY
    partSizes = getPartSizes(partFiles)
    partSize  = findPartSize(partSizes)
    logging.info("Command-line args:")
    logging.info(args)

    glacier = boto3.resource('glacier')
    
    vault = glacier.Vault(args.account_id, args.vault_name)
    logging.info("Opened vault with arn: {}".format(vault.vault_arn))

    hashes = calcAllHashes( partFiles )
    logging.debug("All part hashes:")
    logging.debug(hashes)
    fullHash = joinHashesAsTree( [hashes[x] for x in partFiles] )
    logging.debug("Root hash:")
    logging.debug(fullHash)
    

    logging.debug("Starting multipart upload...")
    multipart_upload = vault.initiate_multipart_upload(archiveDescription=args.archive_description, partSize=str(partSize) )
    logging.info("multipart-id: {}".format(multipart_upload.multipart_upload_id))

    try:
        startPos = 0
        for n in range(len(partFiles)):
            logging.info("-----"*40)
            print("Uploading part {}/{}".format(n+1, len(partFiles) ))
            logging.info("Uploading part {}/{}".format(n+1, len(partFiles) ))
            logging.debug("Part file: {}".format( partFiles[n]) )
            logging.debug("File range: {}-{}".format( startPos, startPos + partSizes[n] - 1 ))
            logging.debug("Part tree hash (calculated internally): {}".format( hashes[partFiles[n]]))

            with open( partFiles[n], "rb" ) as content:
                ret = multipart_upload.upload_part(
                    #range="Content-Range:bytes {}-{}/*".format(startPos, startPos + partSizes[n] - 1),
                    range = "bytes {}-{}/*".format(startPos, startPos + partSizes[n] - 1),
                    body  = content )
                logging.debug("Multipart response:")
                logging.debug(ret)
                assert(ret['ResponseMetadata']['HTTPStatusCode'] == 204)
                externalChecksum = codecs.decode(ret['checksum'], encoding='hex_codec')
                logging.debug("Got external checksum: {}".format(externalChecksum))
                assert(externalChecksum == hashes[partFiles[n]])
            
            startPos += partSizes[n]
    
        if not args.dry_run:
            logging.info("Completing upload")
            ret = multipart_upload.complete(
                archiveSize = str(sum(partSizes)),
                checksum = codecs.decode( codecs.encode( fullHash, encoding='hex_codec' ), 'ascii' )
            )
            assert(ret['ResponseMetadata']['HTTPStatusCode'] == 201)
            logging.info(ret)
            
        else:
            logging.info("Dry-run requested, aborting upload")
            ret = multipart_upload.abort()
            logging.info(ret)
            
    except Exception as e:
        logging.error("Caught exception, aborting upload")
        
        ret = multipart_upload.abort()
        logging.info(ret)

        logging.error(e)
        logging.error(str(e))
        logging.error(e.args)
        print(e)
        print(str(e))
        print(e.args)
        
        raise e
        
    
    #multipart_upload = glacier.MultipartUpload('account_id','vault_name','id')

def standalone():
    import argparse

    argsParser = argparse.ArgumentParser()
    argsParser.add_argument("--parts", type=str, required=True)
    argsParser.add_argument("--account-id", type=str, required=True)
    argsParser.add_argument("--vault-name", type=str, required=True)
    argsParser.add_argument("--archive-description", type=str, required=True)
    argsParser.add_argument("--dry-run", action="store_true", default=False)
    args = argsParser.parse_args()

    up( args )

    

if __name__=="__main__":
    import sys
    sys.exit(standalone())
