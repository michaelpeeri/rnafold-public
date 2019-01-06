import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns; sns.set()


plt.title(  u"\u0394LFE" ) 
plt.ylabel( u"\u0394LFE" )
plt.xlabel( r"$\delta LFE" )

plt.savefig("test_matplotlib_unicode.out.pdf")
plt.savefig("test_matplotlib_unicode.out.png")
plt.close()
