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


