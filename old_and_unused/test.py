import numpy as np
import pandas as pd
import config
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt



hist = [ 5400, 13650, 18600,  1200, 29700, 29400, 23550, 12750,   450,  5700,  1500,   750,
         150,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
         0, 0]
histBins = [ 0.24,  0.26,  0.28,  0.3,   0.32,  0.34,  0.36,  0.38,  0.4,   0.42,  0.44,  0.46,
             0.48,  0.5,   0.52,  0.54,  0.56,  0.58,  0.6,   0.62,  0.64,  0.66,  0.68,  0.7,
             0.72,  0.74]

#histBins = np.arange(0.24, 0.76, 0.02)*100
#    hist, _ = np.histogram( medianGCContent, bins = histBins )
#print(len(medianGCContent))
print(hist)
print(histBins)

fig = plt.figure()
#plt.hist( np.asarray(hist), histBins, normed=0, histtype='stepfilled', range=(20, 80))
#plt.hist( np.asarray(hist), bins=histBins*100.0, histtype='bar', range=(20, 80))
plt.bar( histBins, hist, width=0.02)
plt.xticks(np.arange(0.24, 0.76, 0.04))
plt.savefig("mediangc_%s.pdf" % "test")
plt.savefig("mediangc_%s.svg" % "test")
plt.close(fig)
