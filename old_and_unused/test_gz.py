import gzip
from cStringIO import StringIO
import random
from runningstats import RunningStats


def build_x(length = 2210-40):
    out = []
    out.append("{")
    out.append("windows:[")

    for x in range(length):
        out.append("%.3g," % (random.random()*-100,))

    out.append("]}")

    return "".join(out)

#with gzip.open('test.txt.gz', 'wb') as f:
#     for i in range(10):
#         f.write('111222333444555666777888999000'*100)

bstat = RunningStats()
sstat = RunningStats()

def test():
    big = build_x()
    bstat.push(len(big))

    sio = StringIO()
    f = gzip.GzipFile("", "wb", 9, sio)
    f.write(big)
    f.close()
    
    sstat.push(len(sio.getvalue()))
    sio.close()

for i in range(1005):
    test()

print("%.4g %.3g +-%.3g %.4g" % (bstat.min(), bstat.mean(), 2*bstat.stdev(), bstat.max()))
print("%.4g %.3g +-%.3g %.4g" % (sstat.min(), sstat.mean(), 2*sstat.stdev(), sstat.max()))

