import sys
from datetime import datetime
import data_helpers

# config
queueTag = "rna-fold-window-40-0"

command = sys.argv[1]

rl = data_helpers.RateLimit(30)


active_queue = data_helpers.queueItemsKey % queueTag
suspended_queue = "queue:tag:suspended-%s:members" % queueTag
source_queue = None
dest_queue = None

if( command == "suspend" ):
    source_queue = active_queue
    dest_queue = suspended_queue
elif( command == "resume" ):
    source_queue = suspended_queue
    dest_queue = active_queue
else:
    assert(False)

itemsMoved = 0
while( True):
    item = data_helpers.r.rpoplpush( source_queue, dest_queue )
    if( item is None ):
        break

    itemsMoved += 1
    if(rl()):
        print("%s - moved %d items, %d items remaining" % (datetime.now(), itemsMoved, data_helpers.r.llen( source_queue ) ))

print("Moved %d items" % itemsMoved)
print("Active queue contains %d items" % data_helpers.r.llen( active_queue ) )
print("Suspended queue contains %d items" % data_helpers.r.llen( suspended_queue ) )
    
