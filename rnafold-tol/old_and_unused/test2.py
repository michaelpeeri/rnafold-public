from __future__ import print_function
import sys
import time
import threading


class Per(threading.Thread):

    def __init__(self):
        super(Per, self).__init__()
        self._stop = threading.Event()

    def run(self):
        i = 0

        while(True):
            if(self._stop.isSet()):
                return

            i += 1

            if( i==1 ):
                print("hello")
            elif( i==5 ):
                i = 0


            self._stop.wait(timeout=5.0)

    def stop(self):
        self._stop.set()


t = Per()
t.start()


try:

    i = 0
    while(i<15):
        print(i)
        i += 1
        time.sleep(1.1)

        if( i==11 ):
            raise Exception("done")
finally:
    t.stop()
    t.join()


