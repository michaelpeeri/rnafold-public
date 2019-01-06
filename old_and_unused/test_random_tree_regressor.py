import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use("ggplot") # Use the ggplot style
import seaborn as sns



Nfeatures = 10

X, y = make_regression(n_features=Nfeatures, n_informative=1, n_samples=5000,
                       random_state=None, shuffle=False)

regr = RandomForestRegressor(n_estimators=100, max_depth=10, random_state=0)
regr.fit(X, y)



print(regr.feature_importances_)

print(regr.predict([[0]*Nfeatures]) )
#print(regr.decision_path([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]))

data = np.vstack((X[:,0], y))

Npoints = 500
xs = np.hstack( (np.expand_dims(np.linspace(min(X[:,0]), max(X[:,0]), Npoints ), 1), np.zeros( (Npoints,Nfeatures-1)) ) )
ys = regr.predict(xs)

plt.plot( X[:,0], y )
plt.plot( xs[:,0], ys )

plt.grid(True)
plt.savefig("test_random_tree_regressor.out.pdf")
plt.savefig("test_random_tree_regressor.out.svg")
plt.close()
