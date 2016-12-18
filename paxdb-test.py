import sys
import csv
import requests
from time import sleep
from data_helpers import SpeciesCDSSource


f = sys.argv[1]
taxid = int(sys.argv[2])

#http://rest.ensemblgenomes.org/xrefs/id/CAB12275?content-type=application/json;all_levels=1

server = "http://rest.ensemblgenomes.org"
ext = "/xrefs/id/%s?content-type=application/json;all_levels=1"


def getRecord(protid):
        r = requests.get(
            server + ext % protid)
            #headers={ "Content-Type" : "application/json"})
 
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        # Rate-limit requests
        sleep(0.5)
 
        decoded = r.json()
        return decoded

for protid in SpeciesCDSSource(taxid):
    record = getRecord(protid)
    print(record)
        

#with open(f, 'r') as csvfile:
#    for row in csv.reader(csvfile, delimiter='\t'):
#        #['3706992', '224308.Bsubs1_010100004063', '4.46']

#        paxId = row[1].split(".")[1]
#        pa = float(row[2])

        
        

