from __future__ import print_function
from random import shuffle

def mergedIntervals(intervals):
    completedTo = 0
    #prevInterval = None

    sortedIntervals = list(intervals)
    sortedIntervals.sort( key=lambda x:x[0] )  # sort by start point

    mergedInterval = [-1,-1, ""]
    out = []

    for interval in sortedIntervals:
        #print("Examining interval %s" % str(interval))

        if interval[0] > mergedInterval[1]:
            # emit prev interval
            if( mergedInterval[0] > -1 ):
                out.append(mergedInterval)
            #print("Completed: ", mergedInterval)

            # start a new merged interval
            mergedInterval = list(interval)
        else:
            mergedInterval[1] = interval[1]
            mergedInterval[2] = "%s; %s" % (mergedInterval[2], interval[2])
            #print("Merged: ", mergedInterval)

        #print (prevInterval, interval)
        #assert(prevInterval is None or interval[0] >= prevInterval[0])
        

        #prevInterval = interval

    out.append(mergedInterval)
    #print("Merged: ", mergedInterval)
    return out


def test1(intervals):
    print("-------------------------------------")
    t = list(intervals)
    shuffle(t)
    print("In: ", t)
    ret = mergedIntervals(t)
    print("Out: ", ret)
    return 0


def tests_main():
    test1(((10,20, "a"), (30,40, "b"), (50,60, "c"), (70,80, "d"), (90, 100, "e")))

    test1(((10,20, "a"), (18,30, "b"), (28,40, "c"), (52,60, "d"), (58, 100, "e"), (54,55, "d2"), (20,22,"b2")))

    test1(((1,100, "a"), (18,30, "b"), (28,40, "c"), (52,60, "d"), (58, 100, "e"), (54,55, "d2"), (20,22,"b2")))


if __name__=="__main__":
    import sys
    sys.exit(tests_main())

