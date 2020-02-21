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
# Enqueue failed entries (i.e. entries missing the calculation results) for repeated processing.
# Scan all CDS entries from a given taxid; Find entries missing a computation result; Re-insert them into the queue for repeated processing
# Input - taxid
import sys
import redis
import config

# Configuration
queueKey = "queue:tag:awaiting-rna-fold-for-shuffled-0:members"
calculationResultKey = "CDS:taxid:%d:protid:%s:computed:rna-fold-for-shuffled-0:energy"

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)

taxIdForProcessing = int(sys.argv[1])
print("Procesing %d sequences for tax-id %d (%s)..."
    % (r.scard("species:taxid:%d:CDS" % taxIdForProcessing),
    taxIdForProcessing,
    r.get("species:taxid:%d:name" % taxIdForProcessing)))

skipped = 0
selected = 0

# Iterate over all CDS entries for this species
for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
    # Filtering
    # Skip sequences with partial CDS annotations
    if(r.exists("CDS:taxid:%d:protid:%s:partial" % (taxIdForProcessing, protId))):
        skipped += 1
        continue

    # Is the computed folding energy missing (or not a non-positive float?)
    existingEnergy = r.get(calculationResultKey % (taxIdForProcessing, protId))
    if(existingEnergy is None or not( float(existingEnergy) <= 0.0) ):
        selected += 1
        # Enqueue this protein for processing
        r.rpush(queueKey, "%s:%s" % (taxIdForProcessing, protId))

    skipped += 1

print("%d selected, %d skipped (%d total)" % (selected, skipped, selected+skipped))
print("queue contains %d items" % r.llen(queueKey))
