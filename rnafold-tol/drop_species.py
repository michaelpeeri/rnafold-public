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
from time import sleep
from data_helpers import SpeciesCDSSource, CDSHelper, getSpeciesName, countSpeciesCDS, matchCDSKeyNamesSource, r
from rate_limit import RateLimit

taxId = int(sys.argv[1])

rl = RateLimit(10)

if(countSpeciesCDS(taxId)==0):
    print("Species %d (%s) doesn't have any proteins..." % (taxId, getSpeciesName(taxId)))
    print("Nothing left to do...")
    sys.exit(0)

print("Species %d (%s) has %d proteins stored." % (taxId, getSpeciesName(taxId), countSpeciesCDS(taxId)))
print("Will delete it in 10 seconds...")
sleep(10)

count = 0

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

    count += 1

    if( rl()):
        print("Done processing %d records" % count)

    #if count>100:  # DEBUG ONLY
    #    break      # DEBUG ONLY

keysToDelete = []
print("Enumerating keys...")
c = 0
for k in matchCDSKeyNamesSource(taxId):
    keysToDelete.append(k)
    c += 1
    if( rl()):
        print("%d done" % c)
print("Found %d keys for deletion. Deleting..." % c)

c = 0
for k in keysToDelete:
    r.delete(k)
    c += 1
    if( rl()):
        print("%d done" % c)

print("Species %d (%s) now has %d proteins stored." % (taxId, getSpeciesName(taxId), countSpeciesCDS(taxId)))
