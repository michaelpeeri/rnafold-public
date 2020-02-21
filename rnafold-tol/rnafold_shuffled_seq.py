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
# Worker process for the shuffled sequences queue
# Read sequences queued for processing, do the processing for each, and store the result in the sequence's entry
import redis
import RNA
import config

# Configuration
queueTag                = "queue:tag:awaiting-rna-fold-for-shuffled-0:members"
sequenceTag             = "CDS:taxid:%d:protid:%s:computed:cds-shuffled-seq"
computationResultTag    = "CDS:taxid:%d:protid:%s:computed:rna-fold-for-shuffled-0:energy"


r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)


class XException(object):
    pass

while(True):
    try:
        itemToProcess = None
        # Remove an item from the queue.
        # Note: Failed sequences (i.e. those removed from the queue but not process successfully) will be requeued using
        # the requeue_sequences scripts.
        itemToProcess = r.lpop(queueTag)
        if( itemToProcess is None):
            print("No more items to process...")
            break

        # Each entry in the queue is in the format "taxid:protid"
        print("Processing item %s..." % itemToProcess)
        taxId, protId = itemToProcess.split(":")
        taxId = int(taxId)

        # Get the shuffled sequence for this entry
        seq = r.get(sequenceTag % (taxId, protId))
        # Calculate the RNA folding energy. This is the computation-heavy part.
        strct, energy = RNA.fold(seq)

        # Store the calculation result
        print("%d:%s --> %f" % (taxId, protId, energy))
        r.set(computationResultTag % (taxId, protId), energy)

    except XException:
        print(str(Exception))
        break

print("Done!")
# No more items in the queue; exit...
