import numpy as np
#from sklearn import decomposition
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
from scipy.interpolate import UnivariateSpline


x = np.linspace(0.0, 1.0, 150)
y = 2*np.sin((x+0.5)*5)/(x+0.1)
s = UnivariateSpline(x, y, k=3 )

print(s.get_coeffs())
derivs = s.derivatives(0.0)
print(derivs)



xs = np.linspace(0.0, 1.0, 1000)
ys = s(xs)

ys2 = derivs[0] + xs*derivs[1] + np.power(xs, 2)*derivs[2]/2 + np.power(xs,3)*derivs[3]/6


plt.plot(x,y, '.-', label='data')
plt.plot(xs, ys, zorder=10,  label='smoothed')
plt.plot(xs, ys2, alpha=0.5, zorder=0, label='taylor')
plt.ylim([-6,12])
plt.legend()


plt.savefig("spline_test.pdf")
plt.savefig("spline_test.svg")
plt.close()
