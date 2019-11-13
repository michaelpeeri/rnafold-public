from __future__ import print_function
import random
import itertools
from sqlalchemy import sql
from sqlalchemy.sql.expression import func
import mysql_rnafold as db


session = db.Session()


def randseq():
    return ''.join([c for c in itertools.starmap(lambda:random.choice('acgt'), itertools.repeat((), 50000))])

def insertSequence():
    s = db.Sequence(sequence=randseq(), taxid=0, alphabet=0, source=0)
    session.add(s )
    session.commit()
    return s.id


firstId = insertSequence()

for i in range(8):
    s1 = db.Sequence(sequence=randseq(), taxid=0, alphabet=0, source=0)
    session.add(s1)
session.commit()
    #print(s1.id)

lastId = insertSequence()

#q = session.query(db.Sequence.id).filter(db.Sequence.id==lastId)
#print(q)
#print(q.one())

s = db.connection.execute( sql.select(( sql.func.length(db.sequences.c.sequence),)).where(db.sequences.c.id==lastId) ).scalar()
print(s)

s = db.connection.execute( sql.select(( sql.func.count('*'),)).select_from(db.sequences).where(db.sequences.c.id==lastId) ).scalar()
print(s)


for i in range(firstId, lastId+1):
    for n in range(100):
        f  = db.SequenceSeries( source=0, ext_index=0, index=n, sequence_id=i, value=random.uniform(-100.0, 0.0) )
        session.add(f)
    session.commit()
