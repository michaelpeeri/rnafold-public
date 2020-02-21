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
import numpy as np
import pandas as pd
from data_helpers import getSpeciesName, allSpeciesSource # CDSHelper, SpeciesCDSSource, getSpeciesName, getSpeciesTranslationTable, getSpeciesGenomeAnnotationsFile, getSpeciesGenomeAnnotationsVariant
from process_series_data import readSeriesResultsForSpecies, convertResultsToMFEProfiles, sampleProfilesFixedIntervals, profileElements #, profileLe
from mfe_plots import plotProfilesBoxplots
from mysql_rnafold import Sources
from ncbi_taxa import getDomainForSpecies
from pcache import pcache


# configuration
OutputPositions   = (6,31)
DisplayPositions  = (-250,0)
speciesToExclude = frozenset((405948,999415,946362,470, 158189, 1307761, 456481, 505682, 1280, 4932, 508771, 2850, 753081 ))

@pcache("series-results")
def readSeriesResultsForSpecies_cached( seriesSourceNumber, species, minShuffledGroups, maxShuffledGroups, shuffleType, cdsFilter=None, returnCDS=True ):
    return list(readSeriesResultsForSpecies(seriesSourceNumber, species, minShuffledGroups, maxShuffledGroups, shuffleType, cdsFilter=cdsFilter, returnCDS=returnCDS ))



def getGeneNativeLFEProfiles(taxId, args):
    for result in sampleProfilesFixedIntervals(
            convertResultsToMFEProfiles(
                readSeriesResultsForSpecies_cached( (args.computation_tag,), taxId, 0, 0, shuffleType=args.shuffle_type )
                , 0)
            , args.profile[3], args.profile[0], args.profile[1], args.profile[2]):
        nativeLFE   = result['profile-data'][0]
        yield nativeLFE
    

def processGenome(taxId, args, outputPositions=OutputPositions, displayPositions=DisplayPositions):
    buffer = []
    for result in getGeneNativeLFEProfiles(taxId, args):
        buffer.append(result)

    allLFEs = np.stack(buffer)

    #mean = allLFEs.mean(axis=0)
    #std  = allLFEs.std( axis=0)
    mean = np.apply_along_axis( lambda x: x[~np.isnan(x)].mean(), axis=0, arr=allLFEs )
    std  = np.apply_along_axis( lambda x: x[~np.isnan(x)].std(),  axis=0, arr=allLFEs )

    vals = []
    
    vals.append( ('Species', 
                  getSpeciesName(taxId),
                  'str') )
    vals.append( ('Domain', 
                  getDomainForSpecies(taxId),
                  'str') )
    
    for pos, displayPos in zip( outputPositions, displayPositions ):
        vals.append( ('Mean at {}'.format(displayPos),
                      mean[pos],
                      'float') )
        vals.append( ('Std at {}'.format(displayPos),
                      std[pos],
                      'float') )

    ssetup = [(label, pd.Series( [value], index=[taxId], dtype=dtype)) for label, value, dtype in vals ]
                      
    return pd.DataFrame( dict(ssetup), index=[taxId] )



def parseProfileSpec():
    def convert(value):
        o = value.split(':')
        assert(len(o) >= 3 and len(o) <= 4)
        
        o[0] = int(o[0])
        assert(o[0]>0)
        
        o[1] = int(o[1])
        assert(o[1]>0)
        
        assert(o[2]=="begin" or o[2]=="end")

        if( len(o) == 4 ):
            o[3] = int(o[3])
        else:
            o.append(0)
        
        return (o[0], o[1], o[2], o[3])
    return convert


def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert
    
if __name__=="__main__":
    import sys
    import argparse
    import os.path
    
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument( "--taxid",           type=parseList(int) )
    argsParser.add_argument( "--all-taxa",        default=False, action="store_true")
    argsParser.add_argument( "--profile",         type=parseProfileSpec(), default=parseProfileSpec()('310:10:end:0') )
    argsParser.add_argument( "--computation-tag", type=int,                default=Sources.RNAfoldEnergy_SlidingWindow40_v2 )
    argsParser.add_argument( "--shuffle-type",    type=int,                default=Sources.ShuffleCDSv2_python )
    #argsParser.add_argument( "--Ecoli-workaround",   default=False, action="store_true" )
    args = argsParser.parse_args()

    taxonsToProcess = []
    if not args.all_taxa:
        taxonsToProcess = args.taxid
    else:
        taxonsToProcess= frozenset( allSpeciesSource() ) - speciesToExclude

    assert(taxonsToProcess)
    print("Processing {} taxons".format(len(taxonsToProcess)))

    out = pd.DataFrame()
    for taxId in taxonsToProcess:
        df = processGenome(taxId, args)
        print(df)
        out = out.append( df )

    # write data to file
    out.to_csv( "output_native_LFE_statistics_table.csv" )
    # print domain stats
    out = out.append( out.assign(Domain=lambda x:"All") )
    print( out.assign(Count=lambda x:1).groupby('Domain').sum().Count )
    print( out.groupby('Domain').mean() )

    plotProfilesBoxplots( out.assign(domain=lambda x:x.Domain), positions=tuple(['Mean at {}'.format(x) for x in DisplayPositions ]) )

    sys.exit()
        
