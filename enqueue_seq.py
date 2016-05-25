import sys
import redis
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)

if( len(sys.argv)==2):
    seqid = ""
    seq = sys.argv[1]
elif( len(sys.argv)>=3):
    seqid = sys.argv[1]
    seq = sys.argv[2]

r.set("sequences/")
r.lpush("queue/rnafold/items", seq)
