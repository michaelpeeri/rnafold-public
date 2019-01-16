from mysql_rnafold import Sources
from process_series_data import readSeriesResultsForSpeciesWithSequence, convertResultsToMFEProfiles, sampleProfilesFixedIntervals, profileLength, profileElements, MeanProfile, calcSampledGCcontent


for result in sampleProfilesFixedIntervals(
        convertResultsToMFEProfiles(
            readSeriesResultsForSpeciesWithSequence(
                (Sources.RNAfoldEnergy_SlidingWindow40_v2,),
                85962,
                20,
                20 )
            , 20)
        , 150, 600, 10):
    
    fullCDS = result["cds-seq"]
    seq = fullCDS[150:600]
    
    print(seq)

