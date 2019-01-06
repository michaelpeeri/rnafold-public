from __future__ import print_function
import sys
import codecs
import redis
import config
import cPickle as pickle

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)

value = -3.1419265

a1 = str(value)
a2 = pickle.dumps(value, pickle.HIGHEST_PROTOCOL)

print("Storing as string '%s': %d bytes" % (str(a1), len(a1)))
print("Storing as pickle '%s': %d bytes" % (str(a2), len(a2)))
print(type(a2))


