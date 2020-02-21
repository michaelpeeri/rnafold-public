# RNAFold - Analyze mRNA folding bias (Local Folding Energy) based on randomizations.
# Copyright (C) 2016-2020 Michael Peeri
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
