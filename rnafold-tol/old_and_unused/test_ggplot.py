import numpy as np
from ggplot import *
import pandas as pd




def plotProfile(taxid, data):
    #plt = ggplot(data, aes(y="native", x="position")) + geom_line()
    plotdata = pd.melt(data, id_vars=["position"], value_vars=["native", "shuffled"])
    print(plotdata)
    plt = ggplot(plotdata, aes(y="value", x="position", color="variable")) \
        + geom_line() \
        + xlim(1,150)
    ggsave(filename="test.pdf", plot=plt)
    ggsave(filename="test.svg", plot=plt)


def plotXY(xvals, yvals, labels):
    pass




#xvals = pd.DataFrame((1,5, 2,3,0,1,8), dtype="int32")
#yvals = np.asarray((2,3,-3,3,2,6,2), dtype="int32")
#labels = ("A", "B", "C", "D", "E", "F", "G")

#plotXY(xvals, yvals, labels)


native = np.random.rand(150)*-5 -5
shuffled = np.random.rand(150)*-3 -2
gc = np.random.rand(150)*0.1+0.45
xrange = np.arange(1,151,1)
df = pd.DataFrame( { "native": native, "shuffled": shuffled, "gc": gc, "position": xrange}, index=xrange )
#df.head()
print(pd.melt(df, id_vars=("position",)))

plotProfile(1234, df)
