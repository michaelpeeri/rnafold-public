import argparse
import string
from random import choice, randint
from datetime import datetime
from data_helpers import session, SequenceSeries2Updates, IntegrityError
from rate_limit import RateLimit

argsParser = argparse.ArgumentParser()
argsParser.add_argument("--count", type=int, default=1000000)
argsParser.add_argument("--bulk", type=int, default=200)
args = argsParser.parse_args()

rl = RateLimit(30)

def randstr():
    return ''.join(map(lambda _: choice(string.lowercase), range(50)))

recordsDone = 0
while( recordsDone < args.count):
    print("g")
    for i in range(args.bulk):
        rec = SequenceSeries2Updates()
        rec.sequence_id = randint(1, 999999)
        rec.source = randint(1, 49)
        rec.content = randstr()
        rec.ext_index = 0
        recordsDone += 1
        if recordsDone >= args.count: break
        session.add(rec)

    print("c")
    try:
        session.commit()
    except IntegrityError as e:
        print(e)
        # Ignore and continue...
        # TODO - improve this...
        
    if( rl()):
        print("%s - Written %d records..." % (datetime.now(), recordsDone))

print("Written %d records..." % recordsDone)
    
