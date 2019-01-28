from builtins import str
from builtins import object
import sys
from sqlalchemy import create_engine, Table, Column, Integer, Text, String, SmallInteger, Float, MetaData
from sqlalchemy.dialects.mysql import BLOB
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import config

if config.run_without_mysql_server:
    db = None
else:
    db = create_engine( config.mysql_host_connection, encoding='latin1', echo=False, echo_pool=True, pool_recycle=120 )
connection = None


def testHostAvailable(host, port):
    import socket
    s = None
    try:
        s = socket.create_connection((host, port), 2.0)
        s.close()
        return True
    except socket.error as e:
        pass
    except socket.timeout as e:
        pass

    return False

def getMysqlHostFromConnectionString(connString):
    from re import match
    result = match("mysql([^:]+)?://\w+:\w+@([a-zA-z0-9-._]+)(:[0-9]+)?/\w+", connString)
    host = result.group(2)
    portMatch = result.group(3)
    
    if (not portMatch):
        port = 3306
    else:
        assert(portMatch[0]==":")
        port = int(portMatch[1:])
        
    return (host, port)

connectionInfo = getMysqlHostFromConnectionString(config.mysql_host_connection)

if not config.run_without_mysql_server:
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
# sequences = Table("sequences", md,
#                   Column("id", Integer, primary_key=True),
#                   Column("sequence", Text),
#                   Column("alphabet", SmallInteger),
#                   Column("taxid", Integer),
#                   Column("source", Integer))

sequences2 = Table("sequences2", md,
                  Column("id", Integer, primary_key=True),
                  Column("alphabet", SmallInteger),
                  Column("sequence", BLOB),
                  Column("source", Integer))

# sequence_series = Table("sequence_series", md,
#                         Column("sequence_id", Integer, primary_key=True),
#                         Column("value", Float),
#                         Column("source", Integer, primary_key=True),
#                         Column("ext_index", SmallInteger, primary_key=True),
#                         Column("index", Integer, primary_key=True))

sequence_series2 = Table("sequence_series2", md,
                        Column("sequence_id", Integer, primary_key=True),
                        Column("content", BLOB),
                        Column("source", Integer, primary_key=True),
                        Column("ext_index", Integer))

sequence_series2_updates = Table("sequence_series2_updates", md,
                         Column("dummy_id", Integer, primary_key=True),
                         Column("sequence_id", Integer),
                         Column("content", BLOB),
                         Column("source", Integer),
                         Column("ext_index", Integer))


sequence_floats2 = Table("sequence_floats2", md,
                         Column("sequence_id", Integer, primary_key=True),
                         Column("value", Float),
                         Column("source", Integer, primary_key=True))


#######################################################
# Table definitions for ORM use
#######################################################

Base = declarative_base()

# class Sequence(Base):
#     __tablename__ = "sequences"
#     id = Column(Integer, primary_key=True)
#     sequence = Column(Text)
#     alphabet = Column(SmallInteger)
#     taxid = Column(Integer)
#     source = Column(Integer)

class Sequence2(Base):
    __tablename__ = "sequences2"
    id = Column(Integer, primary_key=True)
    alphabet = Column(SmallInteger)
    sequence = Column(BLOB)
    source = Column(Integer)

# class SequenceSeries(Base):
#     __tablename__ = "sequence_series"
#     sequence_id = Column(Integer, primary_key=True)
#     value = Column(Float)
#     source = Column(Integer, primary_key=True)
#     ext_index = Column(SmallInteger, primary_key=True)
#     index = Column(Integer, primary_key=True)

class SequenceSeries2(Base):
    __tablename__ = "sequence_series2"
    sequence_id = Column(Integer, primary_key=True)
    content = Column(BLOB)
    source = Column(Integer, primary_key=True)
    ext_index = Column(Integer)

class SequenceSeries2Updates(Base):
    __tablename__ = "sequence_series2_updates"
    dummy_id = Column(Integer, primary_key=True)
    sequence_id = Column(Integer)
    content = Column(BLOB)
    source = Column(Integer)
    ext_index = Column(Integer)


class SequenceFloats2(Base):
    __tablename__ = "sequence_floats2"
    sequence_id = Column(Integer, primary_key=True)
    value = Column(Float)
    source = Column(Integer, primary_key=True)


# TESTING ONLY #### TESTING ONLY #### TESTING ONLY #### TESTING ONLY #
Session = sessionmaker(bind=db)
#Session = None


class Alphabets(object):
    DNA = 1
    RNA = 2
    RNA_Huff = 3

# Note: this mixes source numbers for the sequences and sequence_series tables
class Sources(object):
    External = 1 # imported sequence
    Computed = 2
    ShuffleCDSv2_matlab = 10
    ShuffleCDSv2_python = 11
    ShuffleCDS_vertical_permutation_1nt = 12
    RNAfoldEnergy_SlidingWindow40_v2 = 102
    RNAfoldEnergy_SlidingWindow40_v2_alt = 103
    RNAfoldEnergy_SlidingWindow40_v2_native_temp = 110
    CDS_length_nt = 201
    PA_paxdb_single_assay_or_weighted_average = 202
    GC_content_all_CDS = 203
    GC_content_codon_pos_1 = 204
    GC_content_codon_pos_2 = 205
    GC_content_codon_pos_3 = 206
    MFE_mean_window_40nt_estimated = 207
    TEST_StepFunction_BeginReferenced = 801
    TEST_StepFunction_EndReferenced   = 802
    
    
