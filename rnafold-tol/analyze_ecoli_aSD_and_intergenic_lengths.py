# For RNAfold SI
from csv import reader
import pandas as pd
import matplotlib
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
import seaborn as sns



aSDfn                           = 'three_prime_regions.fna.aSD.out'
threePrimeIntergenicDistancesfn = '3prime_intergenic_distances.csv'

def readaSD(fn):
    #AAC74669,-3.700000
    ret = {}
    with open(fn) as csvfile:
        for row in reader( csvfile, delimiter=',' ):
            ret[ row[0] ] = float(row[1])
    return ret


def readIntergenicDistances(fn):
    #AAC74669,b1597,+,423,AACAAATTGAGGGTATGACAATG
    ret = {}
    with open(fn) as csvfile:
        for row in reader( csvfile, delimiter=',' ):
            ret[ row[0] ] = int(row[3])
    return ret

intergenicDistances = readIntergenicDistances(threePrimeIntergenicDistancesfn)
aSD                 = readaSD(aSDfn)


df = pd.DataFrame({'threePrimeDist':pd.Series(dtype='int'), 'aSD':pd.Series(dtype='float') }, index=aSD.keys())

for k,v in aSD.items():
    df.loc[k,'aSD'] = v
    
for k,v in intergenicDistances.items():
    df.loc[k,'threePrimeDist'] = v

print(df.head())


fig, ax1 = plt.subplots()
g = sns.distplot(df.aSD, kde=False, bins=(-9,-8,-7,-6,-5,-4,-3,-2,-1,0), ax=ax1)
plt.savefig("additional_controls_aSD_histogram.pdf")
plt.close(fig)

fig, ax1 = plt.subplots()
g = sns.distplot(df.threePrimeDist, kde=False, bins=(-30,0,50,100,150,200,300,400,500,600,700,800,900,1000,1500,2000,2500,3000,4000,5000,6000,7000), ax=ax1)
plt.savefig("additional_controls_3prime_intergenic_histogram.pdf")
plt.close(fig)

fig, ax1 = plt.subplots()
g = sns.jointplot("threePrimeDist", "aSD", data=df, kind="reg")
plt.savefig("additional_controls_3prime_joint.pdf")
plt.close(fig)

