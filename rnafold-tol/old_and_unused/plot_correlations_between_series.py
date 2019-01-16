# unfinished/deprecated
from mfe_plots import loadProfileData



def standalone():
    import os
    import sys
    from glob import glob
    files = sys.argv[1:]
    
    #----------------------------------------------------------
    # Load profiles
    profileSpecs = sys.argv[1].split(",")
    
    allDataFiles = [[x for x in glob(profilesGlob) if os.path.exists(x)] for profilesGlob in profileSpecs]
    

    groups = []
    for group, files in enumerate(allDataFiles):
        (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, _groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics, wilcoxon, peakWilcoxon) = loadProfileData(files, loadWilcoxon=True, loadTransitionPeakWilcoxon=True if group==0 else False)
        groups.append( (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, _groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics, wilcoxon, peakWilcoxon) )
        assert(len(groups[-1])==12)

    print(len(groups))

    
if __name__=="__main__":
    from sys import exit
    exit( standalone() )

