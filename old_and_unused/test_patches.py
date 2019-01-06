import matplotlib #; matplotlib.use('Agg')#; matplotlib.rc('text', usetex=True)
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt


fig, ax = plt.subplots()

patches = []
polygon = Polygon( [[0, 0], [0, 1], [0.5, 0.5], [1, 1], [1, 0]], True )
patches.append(polygon)

p = PatchCollection( patches, alpha=0.5 )
ax.add_collection(p)

plt.show()
