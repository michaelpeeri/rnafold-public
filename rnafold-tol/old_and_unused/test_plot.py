import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


xaxis = np.arange(1,151,1)

def addMargins(x0, x1, margin=0.1):
    if( x0 > x1 ):
        (x0,x1) = (x1,x0)
    assert(x1 >= x0)

    _margin = abs(x1-x0) * 0.1
    y0, y1 = 0, 0
    y0 = x0 - _margin

    y1 = x1 + _margin

    return (y0,y1)


# def plotProfile_(identifier, native, shuffle):
#     #fig = plt.figure()

#     f, (ax1,ax2) = plt.subplots(2, sharex=True)
#     #f.subplots_adjust(height_ratios=1.5)

#     plt.xlabel('Position (nt, window start, from cds)')
#     plt.title(identifier)

#     plt.subplots_adjust(top=0.1,bottom=0.7)

#     ax1.plot(xaxis, native, label="native")
#     ax1.plot(xaxis, shuffle, label="shuffle")
#     ax1.set_title('Mean LFE')
#     ax1.legend()
#     ax1.grid(True)

#     ax2.plot(xaxis, np.random.rand(150,1)*0.1+0.45)

#     plt.savefig("mfe_40nt_cds_%s.pdf" % identifier )
#     plt.close(fig)


def plotProfile(taxid, data):
    fig, (ax1,ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [2, 1]})

    data[['native', 'shuffled']].plot(ax=ax1)

    plt.title(taxid)

    plt.xlabel('Position (nt, window start, from cds)')

    ax1.set_title('Mean LFE')
    ax1.legend()
    ax1.grid(True)

    ax1.set_ylabel('Mean LFE')

    data['gc'].plot(ax=ax2)
    #plt.subplots_adjust(top=0.4, bottom=0.1)
    ax2.set_title("GC%")
    ax2.grid(True)


    plt.savefig("mfe_40nt_cds_%s.pdf" % taxid )
    plt.savefig("mfe_40nt_cds_%s.svg" % taxid )
    plt.close(fig)



def plotXY(xvals, yvals, _labels):
    fig = plt.figure()
    plt.scatter(xvals, yvals )

    for x,y,l in zip(xvals, yvals, _labels):
        plt.annotate(l, (x-0.1, y+0.2))

    plt.xlabel('window start (nt, from cds)')
    plt.ylabel('MFE')
    plt.grid(True)
    plt.savefig("scatter_xy.pdf")
    plt.close(fig)


xvals = np.asarray((1,5, 2,3,0,1,8), dtype="int32")
yvals = np.asarray((2,3,-3,3,2,6,2), dtype="int32")
labels = ("A", "B", "C", "D", "E", "F", "G")

plotXY(xvals, yvals, labels)


native = np.random.rand(150)*-5 -5
shuffled = np.random.rand(150)*-3 -2
gc = np.random.rand(150)*0.1+0.45
xrange = np.arange(1,151,1)
df = pd.DataFrame( { "native": native, "shuffled": shuffled, "gc": gc, "position": xrange}, index=xrange )


plotProfile(1234, df )
