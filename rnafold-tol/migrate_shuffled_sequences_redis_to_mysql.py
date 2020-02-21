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
# Transfer main sequences from redis to MySQL
# Input: 
#
from __future__ import print_function
import sys
import redis
import config
import mysql_rnafold as db

# Configuration
species = (3055, 556484)

# Connect to redis (source DB)
r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)
# Connect to MySQL (destination DB)
session = db.Session()


skipped = 0
selected = 0

for taxIdForProcessing in species:
    # Iterate over all CDS entries for this taxId
    for protId in r.sscan_iter("species:taxid:%d:CDS" % taxIdForProcessing):
        # Filters
        # (no filters defined)

        # Extract the CDS sequence
        cds = r.get("CDS:taxid:%d:protid:%s:computed:cds-shuffled-seq" % (taxIdForProcessing, protId))
        if( cds is None):
            continue
        if( len(cds) == 0 ):
            continue
            

        selected += 1
        # Store the CDS sequence in destination db; get id of new record
        s1 = db.Sequence(sequence=cds, taxid=taxIdForProcessing, alphabet=db.Alphabets.RNA, source=db.Sources.Computed)
        session.add(s1)
        session.commit()

        # Store record-id in redis
        r.set("CDS:taxid:%d:protid:%s:computed:cds-shuffled-seq-id" % (taxIdForProcessing, protId), s1.id)
        # Remove the sequence itself from redis
        r.delete("CDS:taxid:%d:protid:%s:computed:cds-shuffled-seq" % (taxIdForProcessing, protId))

print("Moved %d sequences." % selected)
