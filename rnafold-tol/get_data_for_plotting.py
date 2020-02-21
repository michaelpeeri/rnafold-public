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
# Input: taxid
# Enqueue all complete CDSs for a particular organism for a given processing task

import sys
import redis
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)

taxIdForProcessing = int(sys.argv[1])
print("#Procesing %d sequences for tax-id %d (%s)..."
    % (r.scard("species:taxid:%d:CDS" % taxIdForProcessing),
    taxIdForProcessing,
    r.get("species:taxid:%d:name" % taxIdForProcessing)))

skipped = 0
selected = 0
for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
    # Filters
    # Skip partial CDSs
    if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
        skipped += 1
        continue

    if(not r.exists("CDS:taxid:%d:protid:%s:computed:rna-fold-0:energy" % (taxIdForProcessing, protId))):
        skipped += 1
        continue

    selected += 1
    # Insert into the processing queue
    cdsLength = int(r.strlen("CDS:taxid:%d:protid:%s:seq" % (taxIdForProcessing, protId)))
    foldEnergy = float(r.get("CDS:taxid:%d:protid:%s:computed:rna-fold-0:energy" % (taxIdForProcessing, protId)))
    print("%s,%d,%f" % (protId, cdsLength, foldEnergy))

print("#%d selected, %d skipped (%d total)" % (selected, skipped, selected+skipped))
