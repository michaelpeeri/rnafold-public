import sys
from sqlalchemy import create_engine, Table, Column, MetaData, sql, String
from sqlalchemy.dialects.sqlite import INTEGER, FLOAT, VARCHAR
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import IntegrityError
import config

_debug = True

def get_db(filename):
    db = create_engine(config.make_sqlite_host_connection(filename), encoding='ascii', echo=_debug)
    connection = db.connect()
    return db, connection


md = MetaData()

sequences_floats = Table("sequences_floats2", md,
                         Column("gene_id", String, primary_key=True),
                         Column("value", FLOAT),
                         Column("source", INTEGER, primary_key=True))

Base = declarative_base()


class SequencesFloats(Base):
    __tablename__ = "sequences_floats2"
    gene_id = Column(String, primary_key=True)
    value = Column(FLOAT)
    source = Column(INTEGER, primary_key=True)



class Store(object):
    def __init__(self, filename):
        self._db, self._connection = get_db(filename)
        self._filename = filename
        md.create_all(self._db)
        self._session = sessionmaker(bind=self._db)()

    def get_float(self, gene_id, source):
        results = self._connection.execute( sql.select( (sequences_floats.c.value, )).select_from(sequences_floats).where(
            sql.and_(
                sequences_floats.c.gene_id==gene_id,
                sequences_floats.c.source==source
                )
        ) )  # Note: order_by not needed, because id is used

        results = self._connection.execute( sql.select( (1,) ).select_from(sequences_floats) )

        results = self._connection.execute( "select value from sequences_floats" )

        print(self._connection)

        print(results.rowcount)
        if( results.rowcount < 1 ):
            print("Error: some results were not found for file=%s, gene_id=%s, source=%d." % (self._filename, gene_id, source))
            return None
        
        return results.fetchall()


    def set_float(self, gene_id, source, float_val):
        sf = SequencesFloats( gene_id=gene_id,
                              value=float_val,
                              source=source)
        self._session.add(sf)
        self._session.commit() # may throw IntegrityError


    
print("*************************")
_db = create_engine(config.make_sqlite_host_connection("test"), echo=_debug)
_connection = _db.connect()
results = _connection.execute( "select * from sequences_floats2")
print(results.rowcount)
results.close()
_connection.close()
print("*************************")
