import time


class PeriodicPrint(object):
    def __init__(self, limitSecs = 60):
        self._limit = limitSecs
        self._last = time.time() - self._limit

    def __call__(self):
        now = time.time()
        if( (now - self._last) > self._limit ):
            self._last = now
            return True
        else:
            return False

i = 0
a = PeriodicPrint(5)
print(a)


while(i<1000):
    if( a() ):
        print(i)
    i += 1
    time.sleep(0.1)
