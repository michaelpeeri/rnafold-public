from data_helpers import db
from process_series_data import readSeriesResultsForSpecies, convertResultsToMFEProfiles, sampleProfilesFixedIntervals
from collections import Counter

# Configuration - what to extract
taxId=511145 # E. coli 
computationTag = db.Sources.RNAfoldEnergy_SlidingWindow40_v2  # 102/103
shuffleType =  db.Sources.ShuffleCDSv2_python # 11
numShuffledGroups = 20
#profile = (1000, 10, "begin", 0)

# for result in sampleProfilesFixedIntervals(
#         convertResultsToMFEProfiles(
#             readSeriesResultsForSpecies((computationTag,), taxId, numShuffledGroups, numShuffledGroups, shuffleType=shuffleType )
#             , numShuffledGroups)
#         , profile[3], profile[0], profile[1], profile[2]):

c = Counter()
for result in readSeriesResultsForSpecies((computationTag,), taxId, numShuffledGroups, numShuffledGroups, shuffleType=shuffleType ):

    # result -> {"content":...., "taxid":..., "cds":}
    # result["content"] -> [{"MFE-profile":..., "Mean-MFE":..., "MeanMFE":..., "seq-crc":..., "id":"taxid:acc_id:seq_id:shuffle_id"}, ...]
    #print(result.keys())
    #profileData = result["profile-data"]
    #print(profileData)
    for result in result["content"]:
        seqId = result["id"]
        seqId = seqId.replace(":", "/").split("/")
        if(len(seqId)) != 4:
            print(seqId)
        print("{},{},{}".format(seqId[1], "native" if int(seqId[3]) == -1 else "randomized", ",".join(map(str, map(lambda x: x if not x is None else "", result["MFE-profile"])))))
        c[seqId[3]] += 1
print(c)



