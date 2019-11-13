import numpy as np
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style


N = 20
fig, ax = plt.subplots()

X_reduced = np.random.random((N,2))*200 - 100


class CenterPreservingNormlizer(matplotlib.colors.Normalize):
    def __init__(self, negativeRange, positiveRange):
        self._negativeRange = negativeRange
        assert(positiveRange>0.0)
        self._positiveRange = positiveRange
        assert(negativeRange<0.0)


        # Declare to the rest of the API what is our output range
        self.vmin=0
        self.vmax=1

    #def __call__(self, values, clip=None):
    def __call__(self, values):
        outHalfScale = 0.5*(float(self.vmax)-float(self.vmin))
        outCenter = self.vmin + outHalfScale

        out = values.copy()

        
        factor = self._positiveRange/outHalfScale
        #print("+factor: %g" % (factor))
        values[values > 0.0] /= factor
        factor = self._negativeRange/outHalfScale*-1
        #print("-factor: %g" % (factor))
        values[values <= 0.0] /= factor

        values += outCenter

        values = 1-values

        assert(np.all(values >= self.vmin))
        assert(np.all(values <= self.vmax))
        
        # Use the logistic function (https://en.wikipedia.org/wiki/Logistic_function) to increase 'contrast'
        # To test using gnuplot:
        # plot [0:1] 1/(1+exp(-15*(x-0.5)))
        #
        steepness = 15.0
        return 1/(1+np.exp(-steepness*(values-0.5)))

cmapNormalizer = CenterPreservingNormlizer(-2.9, 2.9) # TODO - USE REAL RANGE !!!!!!
    
plt.scatter( X_reduced[:,1], X_reduced[:,0], label=[str(i) for i in range(N)] )
scaleX = max(X_reduced[:,1]) - min(X_reduced[:,1])
scaleY = max(X_reduced[:,0]) - min(X_reduced[:,0])
for i in range(N):
    x = X_reduced[i,1]
    y = X_reduced[i,0]
    label = str(i)
    plt.annotate(label, (x+0.05, y-0.05), fontsize=7)
    
    #plt.imshow( np.array( biasProfiles[taxId] ).reshape(1,-1), cmap='coolwarm', norm=cmapNormalizer, extent=(x-scaleX*0.05, x+scaleX*0.05, y-scaleY*0.008, y+scaleY*0.008) )
    #axicon = fig.add_axes( (x-scaleX*0.05, x+scaleX*0.05, y-scaleY*0.008, y+scaleY*0.008) )
    #axicon.imshow( np.array( biasProfiles[taxId] ).reshape(1,-1), cmap='coolwarm', norm=cmapNormalizer )
    #axicon.set_xticks(())
    #axicon.set_yticks(())
    #plt.imshow( np.random.random((10,)).reshape(1,-1)*2.9*2-2.9, extent=(y-scaleY*0.05, y+scaleY*0.05, x-scaleY*0.008, x+scaleY*0.008), alpha=0.3, cmap='coolwarm', norm=cmapNormalizer )
    imw = scaleX*0.08
    imh = scaleY*0.01
    imofy = 4
    plt.imshow( np.random.random((10,)).reshape(1,-1)*1.1*2-1.1, extent=(x-imw, x+imw, -y+imh+imofy, -y-imh+imofy), cmap='coolwarm', norm=cmapNormalizer, interpolation='gaussian', filternorm=1.0, filterrad=100 )
    

ax.set_aspect("equal")
ax.set_xlim((-100,100))
ax.set_ylim((-100,100))
    
plt.grid(False)
plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)

plt.savefig("pca_test.pdf", dpi=800)
plt.savefig("pca_test.png", dpi=100)
plt.close(fig)
