import sys
from time import sleep
from data_helpers import SpeciesCDSSource, CDSHelper, getSpeciesName, countSpeciesCDS

taxId = int(sys.argv[1])

assert(countSpeciesCDS(taxId)>0)
print("Species %d (%s) has %d proteins stored." % (taxId, getSpeciesName(taxId), countSpeciesCDS(taxId)))
print("Will delete it in 10 seconds...")
sleep(10)

for protId in SpeciesCDSSource(taxId):
    print(protId)
    cds = CDSHelper(taxId, protId)
    try:
        cds.dropShuffledSeqs()
    except Exception as e:
        print(e)

    try:
        cds.dropNativeSeq()
    except Exception as e:
        print(e)
        
    cds.dropRecord()
    

print("Species %d (%s) now has %d proteins stored." % (taxId, getSpeciesName(taxId), countSpeciesCDS(taxId)))
