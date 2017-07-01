import _distributed
from dask import delayed
import time
import logging

def f(x):
    logging.warning("f(%d) starting" % x)
    time.sleep(10)
    logging.warning("f(%d) returning" % x)
    return x*10


scheduler = _distributed.open()  # create client, connect to the Dask scheduler

# represent functions as delayed calls
delayed_funcs = []
for x in range(10):
    delayed_funcs.append( delayed(f)(x) )

logging.warning("submitting to scheduler")
futures = scheduler.compute(delayed_funcs)  # submit functions for computation, return futures immediately

_distributed.progress(futures) # wait for computation to complete
print("\n")

x = scheduler.gather(futures) # get the results
print(x)



