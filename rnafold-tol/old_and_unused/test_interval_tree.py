from intervaltree import IntervalTree, Interval


def buildTestTree(N=1000000):
    t = IntervalTree()
    step = N/10

    while(step >= 1):
        print("Step=%d" % step)
        print(".1")
        t2 = IntervalTree.from_tuples( [(i,i+step,None) for i in range(0, N, step)] )
        print(".2")
        t = t | t2
        print(".3")
        step = step/10
    return t

t = buildTestTree()
print("Tree has %d intervals"%len(t))
for i in range(10000):
    t.search(10000,10001)

print(t.search(10010,10020))

