import sys
from sqlalchemy import create_engine, Table, Column, Integer, Text, String, SmallInteger, Float, MetaData
from sqlalchemy.dialects.mysql import BLOB
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import config

db = create_engine(config.mysql_host_connection, encoding='ascii', echo=False)
connection = None


def testHostAvailable(host, port):
    import socket
    s = None
    try:
        s = socket.create_connection((host, 3306), 2.0)
        s.close()
        return True
    except socket.error as e:
        pass
    except socket.timeout as e:
        pass

    return False

def getMysqlHostFromConnectionString(connString):
    from re import match
    host = match("mysql://\w+:\w+@([a-zA-z0-9-._]+)/\w+", connString).group(1)
    port = 3306
    return (host, port)

connectionInfo = getMysqlHostFromConnectionString(config.mysql_host_connection)
if not testHostAvailable(*connectionInfo):
    print("Connection test to mysql server %s failed." % (str(connectionInfo)))
    sys.exit(-1)
    

try:
    connection = db.connect()
except Exception as e:
    print("Failed to connect to mysql database on %s" % config.mysql_host_connection)
    print(e)
    sys.exit(-1)


#######################################################
# Table definitions for non-ORM use
#######################################################

md = MetaData()
sequences = Table("sequences", md,
                  Column("id", Integer, primary_key=True),
                  Column("sequence", Text),
                  Column("alphabet", SmallInteger),
                  Column("taxid", Integer),
                  Column("source", Integer))

sequences2 = Table("sequences2", md,
                  Column("id", Integer, primary_key=True),
                  Column("alphabet", SmallInteger),
                  Column("sequence", BLOB),
                  Column("source", Integer))

sequence_series = Table("sequence_series", md,
                        Column("sequence_id", Integer, primary_key=True),
                        Column("value", Float),
                        Column("source", Integer, primary_key=True),
                        Column("ext_index", SmallInteger, primary_key=True),
                        Column("index", Integer, primary_key=True))

sequence_series2 = Table("sequence_series2", md,
                        Column("sequence_id", Integer, primary_key=True),
                        Column("content", BLOB),
                        Column("source", Integer, primary_key=True),
                        Column("ext_index", Integer))


#######################################################
# Table definitions for ORM use
#######################################################

Base = declarative_base()

class Sequence(Base):
    __tablename__ = "sequences"
    id = Column(Integer, primary_key=True)
    sequence = Column(Text)
    alphabet = Column(SmallInteger)
    taxid = Column(Integer)
    source = Column(Integer)

class Sequence2(Base):
    __tablename__ = "sequences2"
    id = Column(Integer, primary_key=True)
    alphabet = Column(SmallInteger)
    sequence = Column(BLOB)
    source = Column(Integer)

class SequenceSeries(Base):
    __tablename__ = "sequence_series"
    sequence_id = Column(Integer, primary_key=True)
    value = Column(Float)
    source = Column(Integer, primary_key=True)
    ext_index = Column(SmallInteger, primary_key=True)
    index = Column(Integer, primary_key=True)

class SequenceSeries2(Base):
    __tablename__ = "sequence_series2"
    sequence_id = Column(Integer, primary_key=True)
    content = Column(BLOB)
    source = Column(Integer, primary_key=True)
    ext_index = Column(Integer)



Session = sessionmaker(bind=db)


class Alphabets:
    DNA = 1
    RNA = 2
    RNA_Huff = 3

# Note: this mixes source numbers for the sequences and sequence_series tables
class Sources:
    External = 1 # imported sequence
    Computed = 2
    ShuffleCDSv2 = 10
    RNAfoldEnergy_SlidingWindow40 = 102
    RNAfoldEnergy_SlidingWindow40_v2 = 103
    
