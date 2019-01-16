# Extract the CDS sequences for a given species into fasta format (write to standard output)
# Input: taxid
#
from __future__ import print_function
import sys
import redis
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import NucleotideAlphabet
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)


taxIdForProcessing = int(sys.argv[1])


skipped = 0
selected = 0
# Iterate over all CDS entries for this taxId
for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
    # Filters
    # Skip partial CDSs
    if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
        skipped += 1
        continue

    selected += 1
    # Extract the CDS sequence; write in fasta format to standard output
    cds = r.get("CDS:taxid:%d:protid:%s:seq" % (taxIdForProcessing, protId))
    record = SeqRecord(Seq(cds, NucleotideAlphabet), id=protId)
    print(record.format("fasta"), end='')
