# Create correlation plots for Supp. Fig. S2 and S12 - correlations between dLFE profiles calculated using different parameters (randomization type or temperature)
#
# Example (compare randomization types):
# python2 correlation_between_series.py --use-profile-data "gcdata_v2_taxid_*_profile_310_10_begin_0_t11.h5,gcdata_v2_taxid_*_profile_310_10_begin_0_t12.h5"
import os
from glob import glob
import argparse
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
import pandas as pd
from mfe_plots import loadProfileData
from data_helpers import getSpeciesProperty


def getProfilesForAllGroups(args):
    files = [[x for x in glob(profilesGlob) if os.path.exists(x)] for profilesGlob in args.use_profile_data]

    profiles = []

    for group, files in enumerate(files):
        (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics, _1, _2) = loadProfileData(files)
        #self.xdata.append( xdata )
        #self.ydata.append( ydata )
        #self.ydata_nativeonly.append( ydata_nativeonly )
        #self.ydata_shuffleonly.append( ydata_shuffledonly )
        #self.label.append( labels )
        #self.groups.append( groups )
        #self.filesUsed.append( filesUsed )
        print("{} -> {} profiles".format(group, len(biasProfiles)))

        profiles.append(biasProfiles)

    return profiles


def seriesCorrelationAnalysisByPosition(args, profiles, positionsRange=(0,300,10)):
    assert(len(profiles) >= 2)  # we must have two groups to compare
    
    # limit the comparison to species included in both groups
    t1 = frozenset(profiles[0].keys())
    t2 = frozenset(profiles[1].keys())
    taxIdsToCompare = t1.intersection(t2)
    assert( len(taxIdsToCompare) > 2 ) # there must be sufficient overlap to calculate correlation
    print("Comparing {} species".format(len(taxIdsToCompare)))

    print(profiles[0][list(taxIdsToCompare)[0]])

    corr = pd.DataFrame({"Pearson":pd.Series(dtype='float'), "Spearman":pd.Series(dtype='float')}, index=range(positionsRange[0], positionsRange[1], positionsRange[2]))

    plotRange = [positionsRange[0], positionsRange[1], positionsRange[2]]
    if not args.plot_range is None:
        plotRange[0] = args.plot_range[0]
        plotRange[1] = args.plot_range[1]
        
    
    for pos in range(plotRange[0], plotRange[1], plotRange[2]):
        data = []
        
        for profileId in (0,1):
            vals = dict([(taxId, profiles[profileId][taxId][pos]) for taxId in taxIdsToCompare])
            data.append(vals)

        df = pd.DataFrame({"a":pd.Series( data[0] ), "b":pd.Series( data[1] ) } )
        print(df)
        corrForPosition_Pearson  = df.corr(method="pearson")
        assert( corrForPosition_Pearson.iloc[0,1]==corrForPosition_Pearson.iloc[1,0] )
        corr.loc[pos,"Pearson"] = corrForPosition_Pearson.iloc[0,1]
        
        corrForPosition_Spearman = df.corr(method="spearman")
        assert( corrForPosition_Spearman.iloc[0,1]==corrForPosition_Spearman.iloc[1,0] )
        corr.loc[pos,"Spearman"] = corrForPosition_Spearman.iloc[0,1]


    print(corr)
    corr.plot(linewidth=2, alpha=0.7)
    plt.annotate( s="N={}".format(len(taxIdsToCompare)), xy=(255, 0.58), fontsize=14 )
    #plt.annotate( s="N={}".format(len(taxIdsToCompare)), xy=(120, 0.77), fontsize=14 )
    plt.xlim( [ plotRange[0], plotRange[1] ])
    
    plt.title(u"Correlation between {}LFE calculated using \"CDS-wide\" and\n\"position-specific\" randomizations".format(u"\u0394"))
    plt.xlabel('Position (nt, window start, from CDS start)')
    plt.ylabel('Correlation')
    
    plt.savefig("correlation_between_series.out.svg")
    plt.savefig("correlation_between_series.out.pdf")
    plt.savefig("correlation_between_series.out.png")



def seriesCorrelationAnalysisByPosition_ByTrait(args, profiles, positionsRange=(0,300,10)):
    assert(len(profiles) >= 2)  # we must have two groups to compare
    
    # limit the comparison to species included in both groups
    t1 = frozenset(profiles[0].keys())
    t2 = frozenset(profiles[1].keys())
    taxIdsToCompare = t1.intersection(t2)
    assert( len(taxIdsToCompare) > 2 ) # there must be sufficient overlap to calculate correlation
    print("Comparing {} species".format(len(taxIdsToCompare)))

    traitData = getTraitData( args.plot_vs_trait, taxIdsToCompare )
    taxIdsByGroup = []
    tempGroupNames = []
    for r0, r1 in ((20,35),(35,45),(45,75),(75,110)):
        newgroup = [taxId for taxId in taxIdsToCompare if traitData.get(taxId, None)>=r0 and traitData.get(taxId, None)<r1  ]
        taxIdsByGroup.append(frozenset(newgroup))
        tempGroupNames.append(u"Optimum temp {}-{}{}c (N={})".format(r0, r1, u"\xb0", len(newgroup)))
    #print(taxIdsByGroup)
    

    #print(profiles[0][list(taxIdsToCompare)[0]])

    plotRange = [positionsRange[0], positionsRange[1], positionsRange[2]]
    if not args.plot_range is None:
        plotRange[0] = args.plot_range[0]
        plotRange[1] = args.plot_range[1]
        

    fig, ax = plt.subplots()
    c = ["blue", "green", "orange", "red"]
    
    for group in range(4):
        groupName = tempGroupNames[group]

        corr = pd.DataFrame({groupName:pd.Series(dtype='float')}, index=range(positionsRange[0], positionsRange[1], positionsRange[2]))

        for pos in range(plotRange[0], plotRange[1], plotRange[2]):
            data = []

            if len(taxIdsByGroup[group])<4: continue

            for profileId in (0,1):
                vals = dict([(taxId, profiles[profileId][taxId][pos]) for taxId in taxIdsByGroup[group]])
                data.append(vals)

            #print(data)
            if(len(data[0])<4): continue

            print("============= group {} ===== len {} ==========".format(group, len(data[0])))


            df = pd.DataFrame({"a":pd.Series( data[0] ), "b":pd.Series( data[1] ) } )
            #print(df)
            corrForPosition_Pearson  = df.corr(method="pearson")
            assert( corrForPosition_Pearson.iloc[0,1]==corrForPosition_Pearson.iloc[1,0] )
            corr.loc[pos,groupName] = corrForPosition_Pearson.iloc[0,1]
            #corr.loc[pos,"group"] = group
            
            if pos==40:
                print(df)
                print(corrForPosition_Pearson)

            #corrForPosition_Spearman = df.corr(method="spearman")
            #assert( corrForPosition_Spearman.iloc[0,1]==corrForPosition_Spearman.iloc[1,0] )
            #corr.loc[pos,"Spearman"] = corrForPosition_Spearman.iloc[0,1]

        corr.plot(ax=ax, color=c[group], linewidth=2.0, alpha=0.7 )

    #print(corr)
    #corr.plot(ax=ax, by='group')
        
    #plt.annotate( s="N={}".format(len(taxIdsToCompare)), xy=(255, 0.58), fontsize=14 )
    plt.annotate( s="N={}".format(len(frozenset.union(*taxIdsByGroup))), xy=(120, 1.02), fontsize=14 )
    plt.annotate( s=u"Standard temp. = 37{}c".format(u"\xb0"), xy=(10, 1.02), fontsize=14 )
    plt.xlim( [ plotRange[0], plotRange[1] ])
    plt.ylim( [0.55, 1.05] )

    plt.title(u"Correlation between standard-temperature and native-temperature {}LFE profiles".format(u"\u0394"))
    plt.xlabel('Position (nt, window start, from CDS start)')
    plt.ylabel('Pearson\'s correlation')
    
    plt.savefig("correlation_between_series.out.svg")
    plt.savefig("correlation_between_series.out.pdf")
    plt.savefig("correlation_between_series.out.png")
    

def getTraitData(trait, species):
    ret = {}

    for taxId in species:
        optimumTemp = getSpeciesProperty(taxId, 'optimum-temperature')[0]
        if not optimumTemp is None:
            ret[taxId] = float(optimumTemp)

    return ret
    
    

# def seriesCorrelationAnalysisVsTrait(args, profiles):
#     assert(len(profiles) >= 2)  # we must have two groups to compare
    
#     # limit the comparison to species included in both groups
#     t1 = frozenset(profiles[0].keys())
#     t2 = frozenset(profiles[1].keys())
#     taxIdsToCompare = t1.intersection(t2)
#     assert( len(taxIdsToCompare) > 2 ) # there must be sufficient overlap to calculate correlation
#     print("Comparing {} species".format(len(taxIdsToCompare)))

#     #print(profiles[0])
#     #print(profiles[0][list(taxIdsToCompare)[0]])

#     traitData = getTraitData( args.plot_vs_trait, taxIdsToCompare )

#     df = pd.DataFrame({"trait":pd.Series( dtype="float" ), "pos":pd.Series( dtype="int"), "corr":pd.Series( dtype="float" ) } )

#     for pos in (0,60,120):
#         vals_a = dict([(taxId, profiles[0][taxId][pos]) for taxId in taxIdsToCompare])
#         vals_b = dict([(taxId, profiles[1][taxId][pos]) for taxId in taxIdsToCompare])
#         print(vals_a)
#         print(vals_b)


#     #fig, ax1 = plt.subplots()
        

#     #corr = pd.DataFrame({"Pearson":pd.Series(dtype='float'), "Spearman":pd.Series(dtype='float')}, index=range(positionsRange[0], positionsRange[1], positionsRange[2]))

    
    
    


def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert

def standalone():
    import argparse
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument("--use-profile-data", type=parseList(str), required=True)
    argsParser.add_argument("--plot-range", type=parseList(int), required=False)
    argsParser.add_argument("--plot-vs-trait", type=str, required=False)
    args = argsParser.parse_args()

    profiles = getProfilesForAllGroups(args)
    if( args.plot_vs_trait is None ):
        seriesCorrelationAnalysisByPosition( args, profiles )
    else:
        seriesCorrelationAnalysisByPosition_ByTrait( args, profiles )
        
    return 0

    
if __name__=="__main__":
    import sys
    sys.exit(standalone())
        
