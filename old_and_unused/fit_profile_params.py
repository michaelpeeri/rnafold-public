import numpy as np
import scipy.stats


def getGaussianEstimate(xs, ys, center, scale):
    weights = scipy.stats.norm( loc=center, scale=scale )

    sigma = 0.0
    sigmaW = 0.0
    for n in range(1,len(xs)):
        w = weights.cdf(xs[n])-weights.cdf(xs[n-1])
        sigma += ys[n] * w
        sigmaW += w
        if( sigmaW > 0.4999 and w < 1e-6 ):
            break
        #print("n=%d w=%f sigma=%f ys=%f" % (n, w, sigma, ys[n]))

    #print("SigmaW: %f" % sigmaW)
    #print(sigmaW)
    return sigma / sigmaW


def getEstimatedParams(xs, ys, locs, scales):
    estimates = []
    for loc, scale in zip(locs, scales):
        estimates.append( getGaussianEstimate(xs, ys, loc, scale))
    return estimates
