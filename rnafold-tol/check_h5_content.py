# Check h5 files for missing content (to allow make-up runs)
import os
import argparse
from mfe_plots import loadProfileData
from glob import glob

def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert

argsParser = argparse.ArgumentParser()
argsParser.add_argument("--use-profile-data", type=parseList(str), default=())
args = argsParser.parse_args()

allFiles = [[x for x in glob(profilesGlob) if os.path.exists(x)] for profilesGlob in args.use_profile_data]

for group, files in enumerate(allFiles):
    (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics, dLFEwilcoxon, transitionPeakWilcoxon, edgeWilcoxon ) = loadProfileData(files, loadWilcoxon=True, loadTransitionPeakWilcoxon=True, loadEdgeWilcoxon=True)

    print("=========== Group {} ===========".format(group))
    print("Total files: {}".format(len(files)))
    taxIdsWithEdgeWilx = [k for k,v in edgeWilcoxon.items() if not v is None and len(v.shape)==2]
    print("Files with edge-wilcoxon data: {}   {}".format( len(taxIdsWithEdgeWilx), ",".join( map(str, taxIdsWithEdgeWilx) )) )


