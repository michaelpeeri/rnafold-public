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
import sys
import redis
from Bio import SeqIO
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)


#f = open('/home/michael/rnafold/Creinherdtii.json', 'r')
#f = open('/home/michael/rnafold/Ptricornutum.json', 'r')
taxId = int(sys.argv[1])
f = open(sys.argv[2], 'r')

visitedProteinIds = set()


# Store a species entry for this species
#r.set('species:taxid:%d:name' % (taxId,), species)
#r.set('species:name:%s:taxid' % (species,), taxId)

cdsCount = 0
notFoundCount = 0
skippedCount = 0
for record in SeqIO.parse(f, "fasta"):
    proteinId = record.id

    # verify there are no duplicates entries
    assert(proteinId not in visitedProteinIds)
    visitedProteinIds.add(proteinId)

    if( not r.exists("CDS:taxid:%d:protid:%s:seq" % (taxId, proteinId))):
        notFoundCount += 1
        print("Warning: Couldn't find existing entry for protein %s" % proteinId)
        continue

    cdsLength = r.strlen("CDS:taxid:%d:protid:%s:seq" % (taxId, proteinId))
    if( len(record.seq) < cdsLength ): # The cDNA must be at least as long as the CDS
        skippedCount += 1
        print("Warning: Found entry for protein %d, but the cDNA length (%d) is shorter than the CDS length (%d)" % (proteinId, len(record.seq), cdsLength))
        continue


    # Store the cDNA sequence
    r.set('CDS:taxid:%d:protid:%s:cdna-seq' % (taxId, proteinId), record.seq)
    cdsCount += 1

if( notFoundCount + skippedCount > 0):
    print("Warning: %d entries skipped and %d entries not found" % (skippedCount, notFoundCount))

print("Processed %d CDS entries" % (cdsCount,))
print("(out of %d CDS entries for this specied)" % (r.scard("species:taxid:%d:CDS" % (taxId,))))
