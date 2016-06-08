from sqlalchemy import create_engine, Table, Column, Integer, Text, String, SmallInteger, Float, MetaData
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

db = create_engine('mysql://rnafold:rnafold@localhost/rnafold', encoding='ascii', echo=False)
connection = db.connect()

Base = declarative_base()

class Sequence(Base):
    __tablename__ = "sequences"
    id = Column(Integer, primary_key=True)
    sequence = Column(Text)
    alphabet = Column(SmallInteger)
    taxid = Column(Integer)
    source = Column(Integer)

md = MetaData()
sequences = Table("sequences", md,
                  Column("id", Integer, primary_key=True),
                  Column("sequence", Text),
                  Column("alphabet", SmallInteger),
                  Column("taxid", Integer),
                  Column("source", Integer))
sequence_series = Table("sequence_series", md,
                        Column("sequence_id", Integer, primary_key=True),
                        Column("value", Float),
                        Column("taxid", Integer),
                        Column("source", Integer, primary_key=True),
                        Column("index", Integer, primary_key=True))

class SequenceSeries(Base):
    __tablename__ = "sequence_series"
    sequence_id = Column(Integer, primary_key=True)
    value = Column(Float)
    taxid = Column(Integer)
    source = Column(Integer, primary_key=True)
    index = Column(Integer, primary_key=True)


Session = sessionmaker(bind=db)


class Alphabets:
    DNA = 1
    RNA = 2

class Sources:
    External = 1
    Computed = 2
    RNAfoldEnergy_SlidingWindow40 = 102
    
