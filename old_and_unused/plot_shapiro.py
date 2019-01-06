import csv
import matplotlib
matplotlib.use("AGG")  # use raster back-end
import matplotlib.pyplot as plt
import pandas as pd
from math import log10

f = 'shapiropval_taxid_223926.csv'
#f = 'skew_taxid_223926.csv'

df = pd.DataFrame()



with open(f, 'r') as csvfile:
    for row in csv.reader(csvfile, delimiter=','):
        data = list(map(lambda x: None if (x=='None') else float(x), row[1:]))

        df.loc[:, df.shape[1]] = data

fig, ax = plt.subplots()

cax = plt.imshow(df, interpolation='none', cmap=matplotlib.cm.coolwarm)
cbar = fig.colorbar(cax, ticks=[0.0, 0.05, 1.0], orientation='horizontal')
plt.savefig('shapiro.png', dpi=900)
plt.close()



