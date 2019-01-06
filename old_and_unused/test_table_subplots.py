import matplotlib
#from matplotlib.collections import PatchCollection
#from matplotlib.patches import Rectangle
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style


fig, axes = plt.subplots(7, 6) #, gridspec_kw={'height_ratios':[1,1,1,1,1,1,1], 'width_ratios':[1,1,1,1,1,1]})

headings = ["Element", "Refs", "All", "Bacteria", "Archaea", "Eukaryota"]
elements = ["", "Start", "\"Gu\" peak", "Mid-CDS", "End", "End", "\"Model\""]

for i in range(0,7):
    myAxis = axes[i][0]
    myAxis.annotate(elements[i], xy=(0.0, 0.0))
    myAxis.set_xlim([0.0, 1.0])
    myAxis.set_ylim([0.0, 1.0])
    myAxis.set_xticks([])
    myAxis.set_yticks([])
    

for i in range(0,6):
    myAxis = axes[0][i]
    myAxis.annotate(headings[i], xy=(0.0, 0.0))
    myAxis.set_xlim([0.0, 1.0])
    myAxis.set_ylim([0.0, 1.0])
    myAxis.set_xticks([])
    myAxis.set_yticks([])
    
for i in range(1,7):
    for j in range(1,6):
        myAxis = axes[i][j]
        myAxis.barh([0.0], [0.14*i])
        myAxis.set_xlim([0.0, 1.0])
        myAxis.set_xticks([0.0, 0.5, 1.0])
        myAxis.set_yticks([])


plt.savefig("test.png", orientation='portrait', bbox_inches='tight')
plt.close(fig)
