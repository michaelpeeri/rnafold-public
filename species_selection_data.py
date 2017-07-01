import pandas as pd
import numpy as np
import sys


filename = 'Species selection.csv'
df = None

with open(filename, 'r') as csvfile:

    #df = pd.read_csv(fn1, header=None, names=('ProtId','Statistic1', 'Pval1', 'Statistic2', 'Pval2'), na_values='None')
    df = pd.read_csv(csvfile, sep=',', na_values='None', dtype={'Tax-id':np.str})

#print(df[df['Tax-id'].eq('3055')])

def findByTaxid(taxid):
    return df[df['Tax-id'].eq(str(taxid))]

def tests_main():
    # TODO Impl. this
    pass

if __name__=="__main__":
    sys.exit(tests_main())
