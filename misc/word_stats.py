# Counts words (overlapping and non-overlapping) in strings. Multiple strings can be accumulated.
# Note: Multiple strings (a.k.a segments) are not interpreted as concatenated, e.g.:
# a= OverlappingWordStats(3)
# b= NonoverlappingWordStats(3)
# a.addSegment("aaaaa")
# a.addSegment("bbbbb")
# b.addSegment("aaaaa")
# b.addSegment("bbbbb")
# assert("aab" not in a.counts())
# assert("abb" not in a.counts())
# assert("aab" not in b.counts())
# assert("abb" not in b.counts())

class OverlappingWordStats(object):
    def __init__(self, N):
        assert(N>0)
        self._N = N
        self._counts = {}

    def addSegment(self, segment):
        for word in [segment[i:i+self._N] for i in range(0,len(segment)-self._N+1)]:
            if word in self._counts:
                self._counts[word] += 1
            else:
                self._counts[word] = 1

    def counts(self):
        return self._counts


class NonoverlappingWordStats(object):
    def __init__(self, N):
        self._N = N
        self._counts = {}

    def addSegment(self, segment):
        lengthRoundedToWordBoundary = len(segment)-divmod(len(segment),self._N)[1]
        for word in [segment[i:i+self._N] for i in range(0, lengthRoundedToWordBoundary, self._N)]:
            if word in self._counts:
                self._counts[word] += 1
            else:
                self._counts[word] = 1

    def counts(self):
        return self._counts


def tests_nonoverlapping():
    from random import randint, shuffle
    words = ("aaa", "bbb", "ccc", "ddd", "eee")
    test = []
    actualCounts = {}

    for w in words:
        n = randint( 3, 50 )
        for i in range(n):
            test.append(w)
        actualCounts[w] = n

    shuffle(test)

    nws3 = NonoverlappingWordStats(3)
    nws1 = NonoverlappingWordStats(1)

    nws3.addSegment( "".join(test))
    nws1.addSegment( "".join(test))

    assert(nws3.counts() == actualCounts )

    for w,n in nws1.counts().items():
        assert( n == actualCounts[w*3]*3 )

    for w,n in actualCounts.items():
        assert( n*3 == nws1.counts()[w[0]] )

    return 0


def tests_overlapping():
    test = "0123456789abcd123456"
    expectedresult = {"012":1, "123":2, "234":2, "345":2, "456":2, "567":1, "678":1, "789":1, "89a":1, "9ab":1, "abc":1, "bcd":1, "cd1":1, "d12":1}
    ows3 = OverlappingWordStats(3)

    for i in range(1,100):
        ows3.addSegment(test)

        for w,n in ows3.counts().items():
            assert( n == expectedresult[w]*i )

        for w,n in expectedresult.items():
            assert( n*i == ows3.counts()[w] )

    return 0


def tests_overlapping02():
    test = "aaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbcccccccccccccccccccccccccccc"
    expectedresult = {"aaa":19, "aab":1, "abb":1, "bbb":29, "bbc":1, "bcc":1, "ccc":26 }
    ows3 = OverlappingWordStats(3)

    for i in range(1,100):
        ows3.addSegment(test)

        for w,n in ows3.counts().items():
            assert( n == expectedresult[w]*i )

        for w,n in expectedresult.items():
            assert( n*i == ows3.counts()[w] )

    return 0


def tests_main():

    print("Test 1...")
    for i in range(10000):
        ret = tests_nonoverlapping()
        if( ret != 0 ):   return ret

    print("Test 2...")
    for i in range(1000):
        ret = tests_overlapping()
        if( ret != 0 ):   return ret

    print("Test 3...")
    for i in range(1000):
        ret = tests_overlapping02()
        if( ret != 0 ):   return ret

    return 0


if __name__=="__main__":
    import sys
    sys.exit(tests_main())
