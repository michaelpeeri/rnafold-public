# Functions for on-line (streaming, O(1) memory) calculation of mean, variance and stdev
# Refs:
# https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
from math import sqrt
from random import uniform


class RunningStats_alt(object):
    def __init__(self):
        self.clear()

    def clear(self):
        self._n = 0;

    def push(self, x):
        self._n += 1;
        x = float(x)

        if( self._n == 1):
            self._oldM = x
            self._newM = x
            self._oldS = 0.0
        else:
            self._newM = self._oldM + (x - self._oldM) / self._n
            self._newS = self._oldS + (x - self._oldM)*(x - self._newM)
            self._oldM = self._newM
            self._oldS = self._newS
            
    def count(self):
        return self._n

    def mean(self):
        if( self._n > 0):
            return self._newM
        else:
            return 0.0
    
    def variance(self):
        if( self._n > 1):
            return self._newS / (self._n - 1)
        else:
            return 0.0

    def stdev(self):
        return sqrt(self.variance())


class RunningStats(object):
    def __init__(self):
        self.clear()

    def clear(self):
        self._n = 0;
        self._mean = 0.0
        self._m2 = 0.0

    def push(self, x):
        self._n += 1;
        x = float(x)
        
        delta = x - self._mean
        self._mean += delta / self._n
        self._m2 += delta * (x - self._mean)

    def count(self):
        return self._n

    def mean(self):
        if( self._n > 0 ):
            return self._mean
        else:
            raise Exception("Can't return the mean of 0 numbers")
    
    def variance(self):
        if( self._n > 1):
            if( self._n >= 1000 ):
                return self._m2 / (self._n - 1)
            else:
                raise Exception("Refusing to return estimated variance for N<1000. Use exact non-streaming calculation.")
        else:
            return float('nan')

    def stdev(self):
        return sqrt(self.variance())


def ref_mean_std(vals):
    N = len(vals)
    sumX = float(sum(vals))
    meanX = sumX / N

    sumResiduals = 0.0
    for x in vals:
        sumResiduals += (x - meanX)**2
    return( meanX, sqrt(sumResiduals / N))


def test():
    N = 50000
    test = [uniform(-100,100) for x in range(N)]
    a = RunningStats()
    b = RunningStats_alt()
    for i in test:
        a.push(i)
        b.push(i)
    assert(a.count()==N)
    assert(b.count()==N)
        
    test_mean, test_std = a.mean(), a.stdev()
    test2_mean, test2_std = b.mean(), b.stdev()
    ref_mean, ref_std = ref_mean_std(test)

    testOk = True

    if( abs(test_mean - ref_mean) > 1e-8 ):
        testOk = False
    if( abs((test_std - ref_std) / ref_std) > 0.01 ):
        testOk = False
    
    if( not testOk ):
        print("Test values:\t\tmean=%f, std=%f (error: %g%%)" % (test_mean, test_std, (test_std - ref_std) / ref_std * 100))
        print("Test2 values:\t\tmean=%f, std=%f" % (test2_mean, test2_std))
        print("Reference values:\tmean=%f, std=%f" % (ref_mean, ref_std))
        print("N=%d" % N)
        raise Exception("Test failed!")
    else:
        pass

if __name__=="__main__":
    test()
        
