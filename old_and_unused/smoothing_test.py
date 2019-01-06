import numpy as np
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
from scipy.interpolate import PchipInterpolator
from fit_profile_params import getEstimatedParams


#x = np.linspace(0.0, 1.0, 150)
x = np.array([0.0, 20, 60, 120])
y = np.random.uniform(size=(10,4)) * 2 - 1
y[:,3] = y[:,2]

splines = []

for i in range(y.shape[0]):
    problem = y[i,:]
    print(problem)
    #s = InterpolatedUnivariateSpline(x, problem, k=3 )
    s = PchipInterpolator(x, problem, extrapolate=True )
    splines.append(s)

#print(s.get_coeffs())
#derivs = s.derivatives(0.0)
#print(derivs)



xs = np.linspace(0.0, 120.0, 300)

for i in range(y.shape[0]):
    ys = splines[i](xs)
    plt.plot(xs, ys, zorder=10,  label='smoothed_%d' % i)

#ys = s(xs)

#ys2 = derivs[0] + xs*derivs[1] + np.power(xs, 2)*derivs[2]/2 + np.power(xs,3)*derivs[3]/6


samplesX = np.arange(0, 120.0, 1.0)
samplesY = splines[0](samplesX) + np.random.normal(size=len(samplesX))*0.02


def estimateProfileParams(xs, ys):
    return getEstimatedParams(xs, ys, [0.0, 20.0, 80.0], [1.0, 1.0, 8.0])

(p0, p1, p2) = estimateProfileParams(samplesX, samplesY)
print(p0,p1,p2)

plt.plot(samplesX, samplesY, label='simulated_0')

plt.scatter([0.0, 20.0, 60.0], [p0,p1,p2], label='estimated')

#plt.plot(x,y, '.-', label='data')
#plt.plot(xs, ys, zorder=10,  label='smoothed')
#plt.plot(xs, ys2, alpha=0.5, zorder=0, label='taylor')
#plt.ylim([-2,2])
plt.xlim([0,120])
plt.legend()
plt.grid()


plt.savefig("smoothing_test.pdf")
plt.savefig("smoothing_test.svg")
plt.close()

