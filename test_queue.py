import redis
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)

queueTag = "queue:tag:awaiting-rna-fold-0:members"
print("Queue contains %d items" % r.llen(queueTag))


with r.pipeline() as pipe:
  try:
    pipe.watch(queueTag)
    pipe.multi()
    v = pipe.lpop(queueTag)
    print("Got %s" % v)

    raise Exception("oops!")
    pipe.execute()
  except Exception:
    pass
