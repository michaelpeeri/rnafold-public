# Read all computed series data from DB; write summary to hdf5 
# Replaces convert_data_for_plotting.py
#
from builtins import map
from builtins import zip
from builtins import range
from builtins import object
from collections import Counter  # Testing only
import argparse
import csv
import numpy as np
import pandas as pd
from math import log10
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from scipy.stats import wilcoxon, spearmanr
import statsmodels.api as sm
from statsmodels.formula.api import ols
from data_helpers import decompressNucleicSequence, checkSpeciesExist, getSpeciesFileName, allSpeciesSource, getCDSProperty, CDSHelper
from mysql_rnafold import Sources, getWindowWidthForComputationTag
from process_series_data import readSeriesResultsForSpeciesWithSequence, convertResultsToMFEProfiles, sampleProfilesFixedIntervals, profileLength, profileElements, MeanProfile, calcSampledGCcontent, profileEdgeIndex
from mfe_plots import plotMFEProfileWithGC, plotMFEProfileV3, plotXY, scatterPlot, scatterPlotWithColor, plotMFEProfileByPA, plotMFEProfileMultiple, scatterPlotWithKernel, plotMFEProfileForMultipleRandomizations, plot2WayMFEComparison
from codonw import readCodonw, meanCodonwProfile


# ------------------------------------------------------------------------------------
# Command-line args

def parseList(conversion=str):
    def convert(values):
        return list(map(conversion, values.split(",")))
    return convert
    
def parseProfileSpec():
    def convert(value):
        o = value.split(':')
        assert(len(o) >= 3 and len(o) <= 4)
        
        o[0] = int(o[0])
        assert(o[0]>0)
        
        o[1] = int(o[1])
        assert(o[1]>0)
        
        assert(o[2]=="begin" or o[2]=="end" or o[2]=="stop3utr")

        if( len(o) == 4 ):
            o[3] = int(o[3])
        else:
            o.append(0)
        
        return (o[0], o[1], o[2], o[3])
    return convert


class PlottingDestination(object):
    def __init__(self, args, shuffleType, taxid, subclass="all"):
        self.args = args
        self.taxid = taxid
        self.shuffleType = shuffleType
        self.subclass = subclass

        self.reset()

    def accumulate(self, profileData, seq, fullCDS):
        args = self.args

        self.nativeMeanProfile.add( profileData[0,None] )
        self.shuffledMeanProfile.add( profileData[1:] )
        #self.deltaLFEMeanProfile.add( profileData[0,None] - profileData[1:] )


        # Prepare GC profile
        gc = calcSampledGCcontent( seq, args.profile[1] )
        if( gc.size > profileLength(args.profile) ):  # truncate the profile if necessary
            gc = np.resize( gc, (profileLength(args.profile),))
        self.GCProfile.add( np.expand_dims(gc,0) )

        # Prepare percentile mean profiles
        self.shuffled25Profile.add( np.expand_dims( np.percentile(profileData[1:], 25, axis=0), 0) )
        self.shuffled75Profile.add( np.expand_dims( np.percentile(profileData[1:], 75, axis=0), 0) )

        cds_length_nt = len(fullCDS)
        self.cdsLengths.append(cds_length_nt)
        


    def reset(self):
        args = self.args
        profileLen = profileLength(args.profile)
        
        self.nativeMeanProfile   = MeanProfile( profileLen )
        self.shuffledMeanProfile = MeanProfile( profileLen )
        #self.deltaLFEMeanProfile = MeanProfile( profileLen )
        self.shuffled25Profile   = MeanProfile( profileLen )
        self.shuffled75Profile   = MeanProfile( profileLen )
        self.GCProfile           = MeanProfile( profileLen )

        self.cdsLengths = []


    def printDebugInfo(self):
        print("subclass = {}".format(self.subclass))
        print(self.nativeMeanProfile.counts())
        print(self.shuffledMeanProfile.counts())

    def finalize(self):
        # -----------------------------------------------------------------------------
        # Save mean profiles as H5

        # Format (for compatible with plot_xy.py and old convert_data_for_plotting.py:
        #         gc  native  position  shuffled
        # 1    0.451  -4.944         1    -5.886
        # 2    0.459  -5.137         2    -6.069
        # 3    0.473  -5.349         3    -6.262

        args = self.args
        taxid = self.taxid
        shuffleType = self.shuffleType

        xRange = profileElements(args.profile)
        
        df = pd.DataFrame( {
            "native": self.nativeMeanProfile.value(),
            "shuffled": self.shuffledMeanProfile.value(),
            "gc":self.GCProfile.value(),
            "position": xRange,
            "shuffled25":self.shuffled25Profile.value(),
            "shuffled75":self.shuffled75Profile.value()},
            index=xRange )

        statisticsDF = pd.DataFrame({
            'mean_mean_gc': pd.Series([np.mean(self.GCProfile.value())]),
            'taxid': pd.Series([taxid], dtype='int'),
            'cds_count': pd.Series([len(self.cdsLengths)], dtype='int'),
            'media_cds_length_nt': pd.Series([np.median(self.cdsLengths)])
        })
        
        
        if( args.computation_tag == Sources.RNAfoldEnergy_SlidingWindow40_v2 ):
            h5fn = "gcdata_v3_taxid_{}_profile_{}_{}_{}_{}_t{}_{}.h5".format(taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3], shuffleType, self.subclass)
        else:
            h5fn = "gcdata_v3_taxid_{}_profile_{}_{}_{}_{}_t{}_{}_series{}.h5".format(taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3], shuffleType, self.subclass, args.computation_tag)

        # Compression parameters are described here:  http://www.pytables.org/usersguide/libref/helper_classes.html#filtersclassdescr
        # ...and discussed thoroughly in the performance FAQs
        with pd.io.pytables.HDFStore(h5fn, complib="zlib", complevel=1) as store:
            store["df_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = df
            #store["deltas_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = deltasForWilcoxon
            #store["spearman_rho_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = spearman_rho
            store["statistics_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = statisticsDF
            #if( args.codonw ):
            #    store["profiles_spearman_rho_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = dfProfileCorrs
            #store["wilcoxon_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = wilcoxonDf
            #store["transition_peak_wilcoxon_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = guPeakDf
            #store["edge_wilcoxon_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = edgeWilcoxonDf

            store.flush()

        return df

    def getNativeProfile(self):
        return self.nativeMeanProfile.value()
    
    def getShuffledProfile(self):
        return self.shuffledMeanProfile.value()

    #def getDeltaLFEProfile(self):
    #    return self.deltaLFEMeanProfile.value()
    
class ProfilePlot(object):
    def __init__(self, taxId, args):

        self._taxId = taxId
        self._args = args
        
        pa = {}
        if( args.pax_db ):
            with open(args.pax_db, "rb") as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                for row in reader:
                    rank = float(row[2])
                    assert(rank >= 0.0 and rank <= 1.0 )
                    pa[row[0]] = rank

        # Determine the window width
        self.windowWidth = getWindowWidthForComputationTag(args.computation_tag)
                    

    def performPlots(self):
        
        args = self._args
        taxid = self._taxId

        # ------------------------------------------------------------------------------------

        numShuffledGroups = args.num_shuffles
        shuffleTypes = args.shuffle_types
        print("*********** {} ***********".format(shuffleTypes))

        combinedData = {}

        stopCodonFreq = Counter()

        frout = open("out.csv", "wt")

        for shuffleType in shuffleTypes:
            n = 0

            x1 = Counter()
            x2 = Counter()
            x3 = Counter()
            
            print("Processing species %d (%s), shuffleType=%d" % (taxid, getSpeciesFileName(taxid), shuffleType))

            nativeMeanProfile = MeanProfile( profileLength(args.profile) )
            shuffledMeanProfile = MeanProfile( profileLength(args.profile) )

            shuffled25Profile = MeanProfile( profileLength(args.profile) )
            shuffled75Profile = MeanProfile( profileLength(args.profile) )

            h5destination     = PlottingDestination(args, shuffleType, taxid)
            h5destination_tag = PlottingDestination(args, shuffleType, taxid, subclass="endcodon_tag")
            h5destination_taa = PlottingDestination(args, shuffleType, taxid, subclass="endcodon_taa")
            h5destination_tga = PlottingDestination(args, shuffleType, taxid, subclass="endcodon_tga")
            
            h5destination_readthrough_pos = PlottingDestination(args, shuffleType, taxid, subclass="readthrough_pos")
            h5destination_readthrough_neg = PlottingDestination(args, shuffleType, taxid, subclass="readthrough_neg")

            h5destination_next_cds_same     = PlottingDestination(args, shuffleType, taxid, subclass="next_cds_same")
            h5destination_next_cds_opposite = PlottingDestination(args, shuffleType, taxid, subclass="next_cds_opposite")

            h5destination_integenic_short    = PlottingDestination(args, shuffleType, taxid, subclass="intergenic_positive_short")
            h5destination_integenic_long     = PlottingDestination(args, shuffleType, taxid, subclass="intergenic_positive_long")
            h5destination_integenic_overlap  = PlottingDestination(args, shuffleType, taxid, subclass="intergenic_positive_overlap")

            dLFEvsFlankingLength = pd.DataFrame({'flanking_length':pd.Series(dtype='int'), 'dLFE':pd.Series(dtype='float'), 'nextCDSOppositeStrand':pd.Series(dtype='int')})

           
            xRange = profileElements(args.profile)

            nativeMeanProfile_HighPAOnly = None
            nativeMeanProfile_MediumPAOnly = None
            nativeMeanProfile_LowPAOnly = None
            if( args.pax_db ):
                nativeMeanProfile_HighPAOnly = MeanProfile( profileLength(args.profile) )
                nativeMeanProfile_MediumPAOnly = MeanProfile( profileLength(args.profile) )
                nativeMeanProfile_LowPAOnly = MeanProfile( profileLength(args.profile) )

            GCProfile = MeanProfile( profileLength(args.profile) )

            #deltasForWilcoxon = np.zeros((0,2), dtype=float)
            deltasForWilcoxon = pd.DataFrame({'pos':pd.Series(dtype='int'), 'delta':pd.Series(dtype='float')})

            fullDeltas = []

            geneLevelScatter = pd.DataFrame({'gc':pd.Series(dtype='float'), 'logpval':pd.Series(dtype='float'), 'abslogpval':pd.Series(dtype='float'), 'protid':pd.Series(dtype='str')})

            cdsLengths = []

            fullSeqs = []
            dfCodonw = None

            # ------------------------------------
            # Process all CDS for this species
            # ------------------------------------
            for result in sampleProfilesFixedIntervals(
                    convertResultsToMFEProfiles(
                        readSeriesResultsForSpeciesWithSequence((args.computation_tag,), taxid, numShuffledGroups, numShuffledGroups, shuffleType=shuffleType )
                        , numShuffledGroups)
                    , args.profile[3], args.profile[0], args.profile[1], args.profile[2]):

                fullCDS = result["cds-seq"]

                if args.profile[2]=="stop3utr":
                    # TODO FIX THIS
                    seq = fullCDS[args.profile[3]:args.profile[0]]
                else:
                    seq = fullCDS[args.profile[3]:args.profile[0]]

                if not seq:
                    continue

                                    

                stopCodonPos = result['content'][0]['stop-codon-pos']
                assert(stopCodonPos%3==0)
                stopCodon = fullCDS[stopCodonPos:stopCodonPos+3]
                stopCodonFreq.update( (stopCodon, ) )
                
                protId = result["cds"].getProtId()
                #print("Length: {}nt".format(result["cds"].length()))

                fullSeqs.append( SeqRecord( Seq(fullCDS, NucleotideAlphabet), id=protId) )


                profileData = result["profile-data"]
                #print("ppl: {}".format( profileData.shape))
                assert(profileData.shape[0] >= numShuffledGroups)
                #print(profileData.shape)
                #print(profileData)

                #print(profileData[:,0].T)

                #print( profileData[:, [0,99,-1]] )
                #print(profileData.shape)

                # Prepare mean MFE profiles

                h5destination.accumulate( profileData, seq, fullCDS )
                
                if stopCodon=="tag":
                    h5destination_tag.accumulate( profileData, seq, fullCDS )
                    
                elif stopCodon=="taa":
                    h5destination_taa.accumulate( profileData, seq, fullCDS )
                    
                elif stopCodon=="tga":
                    h5destination_tga.accumulate( profileData, seq, fullCDS )


                #
                nextCDSOnOppositeStrand = CDSHelper( taxid, protId ).nextCDSOnOppositeStrand()
                assert(isinstance(nextCDSOnOppositeStrand, bool))
                if nextCDSOnOppositeStrand:
                    h5destination_next_cds_opposite.accumulate( profileData, seq, fullCDS )
                else:
                    h5destination_next_cds_same.accumulate( profileData, seq, fullCDS )

                flanking3UTRRegionLengthNt = result['cds'].flankingRegion3UtrLength()

                if profileData.shape[1] > 105:
                    dLFE_105 = profileData[0,105] - np.mean(profileData[1:,105]) 
                    dLFEvsFlankingLength = dLFEvsFlankingLength.append(
                        pd.DataFrame({'flanking_length':pd.Series([flanking3UTRRegionLengthNt]),
                                      'dLFE':pd.Series([dLFE_105]),
                                      'nextCDSOppositeStrand':pd.Series([int(nextCDSOnOppositeStrand)])}) )
                else:
                    print("Warning: skipping {}...".format( result['cds'].getProtId() ))

                if flanking3UTRRegionLengthNt > 0     and flanking3UTRRegionLengthNt <= 100:
                    h5destination_integenic_short.accumulate( profileData, seq, fullCDS )
                elif flanking3UTRRegionLengthNt > 100 and flanking3UTRRegionLengthNt <= 400:
                    h5destination_integenic_long.accumulate( profileData, seq, fullCDS )
                elif flanking3UTRRegionLengthNt <= 0:
                    h5destination_integenic_overlap.accumulate( profileData, seq, fullCDS )
                    
                
                exid = 0
                experiment_name = ("ribo_MG1655_MOPS_rep1", "ribo_MG1655_MOPS_rep2", "ribo_rich", "WT_rep1", "WT_rep2", "WT_rep3")[exid]
                readthroughVal = getCDSProperty( taxid, protId, "readthrough-v2.ex{}".format(exid) )
                #print("readthrough-v2.ex0:{}".format(readthroughVal))
                if readthroughVal=="1":
                    h5destination_readthrough_pos.accumulate( profileData, seq, fullCDS )
                    
                elif readthroughVal=="0":
                    h5destination_readthrough_neg.accumulate( profileData, seq, fullCDS )

                # Prepare data for genome-wide wilcoxon test
                #newDeltas = profileData[0,0::4] - np.mean(profileData[1:,0::4], axis=0)
                newDeltas = profileData[0,0::1] - np.mean(profileData[1:,0::1], axis=0)
                #print("newDeltas: {}".format(newDeltas.shape))
                #newPositions = range(args.profile[3], profileLength(args.profile), 40)
                newPositions = list(range(args.profile[3], args.profile[0], args.profile[1]))
                deltaspd = pd.DataFrame({'pos':pd.Series(newPositions, dtype='int'), 'delta':pd.Series(newDeltas, dtype='float')})
                #print("deltaspd: {}".format(deltaspd.shape))
                deltasForWilcoxon = deltasForWilcoxon.append(deltaspd)

                fullDeltas.append( profileData[0,0::1] - profileData[1:,0::1] )  # store the 20x31 matrix of deltas for full wilcoxon test

                # Prepare data for GC vs. selection test
                meanGC = calcSampledGCcontent( seq, 10000)[0]
                if( not (meanGC >= 0.05 and meanGC <= 0.95)):
                    meanGC = None
                #deltas = profileData[0,0::4] - np.mean(profileData[1:,0::4], axis=0)
                deltas = profileData[0,0::1] - np.mean(profileData[1:,0::1], axis=0)
                #print("deltas: {}".format(deltas.shape))
                pvalue = wilcoxon(deltas).pvalue
                direction = np.mean(deltas)
                directedLogPval = None

                if( pvalue > 0.0 ):
                    directedLogPval = log10(pvalue) * direction * -1.0
                else:
                    directedLogPval = -250.0      * direction * -1.0

                paval = None
                if( args.pax_db ):
                    paval = pa.get(protId)
                    if( paval >= 0.8 ):
                        nativeMeanProfile_HighPAOnly.add( profileData[0,None] )
                    elif( paval <= 0.2 ):
                        nativeMeanProfile_LowPAOnly.add( profileData[0,None] )
                    elif( not paval is None ):
                        nativeMeanProfile_MediumPAOnly.add( profileData[0,None] )

                cds_length_nt = len(fullCDS)
                cdsLengths.append(cds_length_nt)



                frline = "{}\t{}\t{}".format(protId, stopCodon, deltas)
                frout.write(frline)
                print(frline)
                
                geneLevelScatter = geneLevelScatter.append(pd.DataFrame({'gc':pd.Series([meanGC]), 'logpval': pd.Series([directedLogPval]), 'abslogpval': pd.Series([pvalue]), 'protid':pd.Series([protId]), 'pa':pd.Series([paval]), 'cds_length_nt':pd.Series([cds_length_nt])}))


                x1.update((fullCDS[0],))
                x2.update((fullCDS[1],))
                x3.update((fullCDS[2],))
                del fullCDS

                del result
                n += 1
            del(pvalue); del(direction); del(seq); del(deltas)

            # Refuse to proceed if the data found is unreasonably small
            if( n < 100 ):
                raise Exception("Found insufficient data to process taxid=%d (n=%d)" % (taxid, n))


            CUBmetricsProfile = None
            if( args.codonw ):
                fFullSeqs = NamedTemporaryFile(mode="w")         # create a temporary file
                SeqIO.write( fullSeqs, fFullSeqs.name, "fasta")  # write the full sequences into the file
                dfCodonw = readCodonw( fFullSeqs.name )          # run codonw and get the gene-level results

                print('****************************************************')
                print(dfCodonw.columns)
                print(dfCodonw.head())

                print(geneLevelScatter.columns)
                print(geneLevelScatter.head())

                geneLevelScatter = pd.merge(dfCodonw, geneLevelScatter, left_index=True, right_index=False, right_on='protid')
                print(geneLevelScatter.corr())

                #args.profile[3], args.profile[0], args.profile[1]
                CUBmetricsProfile = meanCodonwProfile(fullSeqs, self.windowWidth, 'begin', args.profile[3], args.profile[0], args.profile[1]) # TODO - use real values!
                print(CUBmetricsProfile)

            #else:
            #    geneLevelScatter['CAI'] = pd.Series( np.zeros(len(geneLevelScatter)), index=geneLevelScatter.index)



            # ------------------------------------
            # Display summary for this species
            # ------------------------------------
            #print("Native:")
            #print(nativeMeanProfile.value())
            #print(nativeMeanProfile.counts())

            #print("Shuffled:")
            #print(shuffledMeanProfile.value())
            #print(shuffledMeanProfile.counts())

            #print(deltasForWilcoxon.shape)

            #------------------------------------------------------------------------------------------------------------------
            # Test for significance of the mean dLFE (postive or negative) at each position along the genome
            # (This will answer questions such as "how many genomes have (significantly) negative dLFE at position X?")
            #------------------------------------------------------------------------------------------------------------------
            # wilcoxonDf = pd.DataFrame({'pos':pd.Series(dtype='int'), 'logpval':pd.Series(dtype='float'), 'N':pd.Series(dtype='int') })
            # if( True ):
            #     print("Processing full deltas...")

            #     # Perform statistical tests based on the deltas for each position (from all proteins)
            #     for pos in range(profileLength(args.profile)):

            #         # Collect all deltas for this position (data will be an list of arrays of length 20 - one for each protein long enough to have deltas at this position)
            #         print("pos: {} len: {}".format( pos, len(fullDeltas)))
            #         data = [x[:,pos] for x in fullDeltas if x.shape[1]>pos]
            #         print(data.shape)
            #         dataar = np.concatenate(data)  # flatten all deltas
            #         assert(dataar.ndim == 1)

            #         # Perform 1-sample Wilcoxon signed-rank test on the deltas (testing whether the deltas are symmetrical)
            #         wilcoxonPval = wilcoxon(dataar).pvalue  # 2-sided test
            #         if wilcoxonPval>0.0:
            #             logWilcoxonPval = log10(wilcoxonPval)
            #         else:
            #             logWilcoxonPval = -324.0 # ~minimum value for log10(0.000.....)

            #         N = dataar.shape[0]
                    
            #         #wilcoxonDf = wilcoxonDf.append(pd.DataFrame({'pos':pd.Series(xRange[pos]), 'N':pd.Series([N]), 'logpval': pd.Series([logWilcoxonPval]) } ) )
                    
            #         #alldeltas = np.concatenate(fullDeltas)
            #     #print(wilcoxonDf)
            #     del(data); del(dataar)

            #------------------------------------------------------------------------------------------------------------------

            #------------------------------------------------------------------------------------------------------------------
            # Find "transition peak"
            #------------------------------------------------------------------------------------------------------------------
            # Calculate the dLFE
            # print("-TransitionPeak-TransitionPeak-TransitionPeak-TransitionPeak-")
            # meanDeltaLFE = nativeMeanProfile.value() - shuffledMeanProfile.value()
            # peakPos = np.argmin( meanDeltaLFE )
            # peakVal = meanDeltaLFE[peakPos]
            
            # guPeakDf = pd.DataFrame({'pos':pd.Series(dtype='int'), 'logpval':pd.Series(dtype='float') })
            
            # if peakVal <= 0.0:
            #     print("{} {}".format(peakPos, peakVal))
            #     if not wilcoxonDf[wilcoxonDf['pos']==peakPos*10].empty:
            #         logpval = wilcoxonDf[wilcoxonDf['pos']==peakPos*10].logpval.loc[0]
            #         print(type(logpval))
            #         #print(logpval.shape)
            #         print(logpval)

            #         if logpval < -2.0:
            #             #

            #             # Collect all deltas for this position (data will be an list of arrays of length 20 - one for each protein long enough to have deltas at this position)

            #             for otherPos in range(profileLength(args.profile)):

            #                 data1 = [x[:,peakPos] for x in fullDeltas if x.shape[1]>max(peakPos,otherPos)]
            #                 peakData = np.concatenate(data1)  # flatten all deltas
            #                 assert(peakData.ndim == 1)

            #                 data2 = [x[:,otherPos] for x in fullDeltas if x.shape[1]>max(peakPos,otherPos)]
            #                 otherData = np.concatenate(data2)  # flatten all deltas

            #                 assert(len(peakData)==len(otherData))
            #                 datax = otherData-peakData

            #                 print("/-: {} {} {}".format(peakPos, otherPos, np.mean(datax)))

            #                 #if( peakPos==otherPos ):
            #                 #    print(datax)

            #                 wilcoxonPval = None
            #                 if np.allclose(otherData, peakData):
            #                     logWilcoxonPval = 0.0
            #                 else:
            #                     # Perform 1-sample Wilcoxon signed-rank test on the deltas (testing whether the deltas are symmetrical)
            #                     wilcoxonPval = wilcoxon(peakData, otherData).pvalue  # 2-sided test (not ideal in this situation...)
            #                     if wilcoxonPval>0.0:
            #                         logWilcoxonPval = log10(wilcoxonPval)
            #                     elif wilcoxonPval==0.0:
            #                         logWilcoxonPval = -324.0 # ~minimum value for log10(0.000.....)
            #                     else:
            #                         logWilcoxonPval = None

            #                 if not logWilcoxonPval is None:
            #                     #guPeakDf = guPeakDf.append(pd.DataFrame({'pos':pd.Series(xRange[otherPos]), 'logpval': pd.Series([logWilcoxonPval]) } ) )


            #             print(guPeakDf)


            #------------------------------------------------------------------------------------------------------------------
            # Calculate edge-referenced wilcoxon
            #------------------------------------------------------------------------------------------------------------------

            # edgePos = profileEdgeIndex(args.profile)
            # data0 = [x[:,edgePos] if x.shape[1]>pos else None for x in fullDeltas]
            # edgeWilcoxonDf = pd.DataFrame({'pos':pd.Series(dtype='int'), 'logpval':pd.Series(dtype='float') })
            
            # for pos in range(profileLength(args.profile)):
            #     data1 = [x[:,pos] if x.shape[1]>pos else None for x in fullDeltas]
            #     assert(len(data0)==len(data1))
            #     if not data1[0] is None:
            #         print("]]]]]]]]]]]]] {}".format(data1[0].shape))

            #     diffs = []
            #     for d0, d1 in zip(data0, data1):
            #         if (not d0 is None) and (not d1 is None):
            #             #print("{} {}".format(d0.shape, d1.shape))
            #             d = d1-d0 
            #             diffs.append( d[~np.isnan(d)] )

            #     alldiffs = np.concatenate( diffs )
            #     #print(alldiffs.shape)
            #     print(pos)
            #     #print(alldiffs[:100])
            #     print(alldiffs.shape)
                
            #     wilcoxonPval = None
            #     if np.allclose(alldiffs, 0):
            #         logWilcoxonPval = 0.0
            #     else:
            #         # Perform 1-sample Wilcoxon signed-rank test on the deltas (testing whether the deltas are symmetrical)
            #         wilcoxonPval = wilcoxon(alldiffs).pvalue  # 2-sided test (not ideal in this situation...)
            #         if wilcoxonPval>0.0:
            #             logWilcoxonPval = log10(wilcoxonPval)
            #         elif wilcoxonPval==0.0:
            #             logWilcoxonPval = -324.0 # ~minimum value for log10(0.000.....)
            #         else:
            #             logWilcoxonPval = None

            #     #if not logWilcoxonPval is None:
            #         #edgeWilcoxonDf = edgeWilcoxonDf.append(pd.DataFrame({'pos':pd.Series(xRange[pos]), 'logpval': pd.Series([logWilcoxonPval]) } ))
            # #print(edgeWilcoxonDf)

            # Count the mininum number of sequences
            minCount = 1000000
            for pos in range(profileLength(args.profile)):
                countAtPos = sum([1 if x.shape[1]>pos else 0 for x in fullDeltas])
                if countAtPos < minCount: minCount = countAtPos

            #------------------------------------------------------------------------------------------------------------------
            # Store the results
            #------------------------------------------------------------------------------------------------------------------
            
            #native = np.asarray(nativeMean[1:], dtype="float")
            #shuffled = np.asarray(shuffledMean[1:], dtype="float")
            #gc = np.asarray(gcMean[1:], dtype="float")
            #xrange = [x for x in args.profile.Elements() if x<profileInfo.cdsLength()]
            profileId = "%d_%d_%s_t%d" % (args.profile[0], args.profile[1], args.profile[2], shuffleType)


            #print(nativeMeanProfile.value())
            #print(xRange)
            
            #df = pd.DataFrame( { "native": nativeMeanProfile.value(), "shuffled": shuffledMeanProfile.value(), "gc":GCProfile.value(), "position": xRange, "shuffled25":shuffled25Profile.value(), "shuffled75":shuffled75Profile.value()}, index=xRange )
            #print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
            #if( not CUBmetricsProfile is None ):
            #    df = pd.merge( df, CUBmetricsProfile, how='left', left_on='position', right_on='windowStart')
            #    df = df.set_index('position')


            df = None
            for dest in (h5destination, h5destination_tag, h5destination_taa, h5destination_tga, h5destination_readthrough_pos, h5destination_readthrough_neg):
                dest.printDebugInfo()
                newDf = dest.finalize()
                if df is None:
                    df = newDf # only save the first df...
                
            
            combinedData[shuffleType] = df
            #print(df)

            dfProfileCorrs = None
            if( args.codonw ):
                plotMFEProfileMultiple(taxid, profileId, df, ('GC', 'Nc', 'CAI', 'CBI', 'Fop'), scaleBar=self.windowWidth)

                smfe = df['native'] - df['shuffled']
                spearman_gc  = spearmanr( df['GC'],  smfe )
                spearman_Nc  = spearmanr( df['Nc'],  smfe )
                spearman_CAI = spearmanr( df['CAI'], smfe )
                spearman_Fop = spearmanr( df['Fop'], smfe )
                dfProfileCorrs = pd.DataFrame( { "spearman_smfe_gc_rho":   spearman_gc.correlation,
                                                 "spearman_smfe_gc_pval":  spearman_gc.pvalue,
                                                 "spearman_smfe_Nc_rho":   spearman_Nc.correlation,
                                                 "spearman_smfe_Nc_pval":  spearman_Nc.pvalue,
                                                 "spearman_smfe_CAI_rho":  spearman_CAI.correlation,
                                                 "spearman_smfe_CAI_pval": spearman_CAI.pvalue,
                                                 "spearman_smfe_Fop_rho":  spearman_Fop.correlation,
                                                 "spearman_smfe_Fop_pval": spearman_Fop.pvalue },
                                                 index=(taxid,) )


            lengthsDist = np.array(cdsLengths)
            statisticsDF = pd.DataFrame({
                'mean_mean_gc': pd.Series([np.mean(GCProfile.value())]),
                'taxid': pd.Series([taxid], dtype='int'),
                'cds_count': pd.Series([len(cdsLengths)], dtype='int'),
                'media_cds_length_nt': pd.Series([np.median(cdsLengths)])
                })


            plotMFEProfileWithGC(taxid, profileId, df, computationTag=args.computation_tag)

            # plot2WayMFEComparison( taxid,
            #                        profileId,
            #                        pd.DataFrame({'native':   h5destination_readthrough_neg.getNativeProfile(),
            #                                      'shuffled': h5destination_readthrough_neg.getShuffledProfile()},
            #                                     index=xRange),
            #                        pd.DataFrame({'native':   h5destination_readthrough_pos.getNativeProfile(),
            #                                      'shuffled': h5destination_readthrough_pos.getShuffledProfile()},
            #                                     index=xRange),
            #                        computationTag=args.computation_tag,
            #                        comment=experiment_name)
            
            # plot2WayMFEComparison( taxid,
            #                        profileId,
            #                        pd.DataFrame({'native':   h5destination_next_cds_same.getNativeProfile(),
            #                                      'shuffled': h5destination_next_cds_same.getShuffledProfile()},
            #                                     index=xRange),
            #                        pd.DataFrame({'native':   h5destination_next_cds_opposite.getNativeProfile(),
            #                                      'shuffled': h5destination_next_cds_opposite.getShuffledProfile()},
            #                                     index=xRange),
            #                        computationTag=args.computation_tag,
            #                        yRange=(-13,1),
            #                        comment="Next CDS on same strand")

            # numGenesWithShortIntergenicRegions = np.min(h5destination_integenic_short.nativeMeanProfile.counts())
            # numGenesWithLongIntergenicRegions  = np.min(h5destination_integenic_long.nativeMeanProfile.counts())
            
            # plot2WayMFEComparison( taxid,
            #                        profileId,
            #                        pd.DataFrame({'native':   h5destination_integenic_short.getNativeProfile(),
            #                                      'shuffled': h5destination_integenic_short.getShuffledProfile()},
            #                                     index=xRange),
            #                        pd.DataFrame({'native':   h5destination_integenic_long.getNativeProfile(),
            #                                      'shuffled': h5destination_integenic_long.getShuffledProfile()},
            #                                     index=xRange),
            #                        computationTag=args.computation_tag,
            #                        yRange=(-13,1),
            #                        comment="Intergenic region lengths\nShort (N>={}) Long (N>={}) Total (N>={}) ".format( numGenesWithShortIntergenicRegions, numGenesWithLongIntergenicRegions, numGenesWithShortIntergenicRegions + numGenesWithLongIntergenicRegions ) )

            numGenesWithShortIntergenicRegions = np.min(h5destination_integenic_short.nativeMeanProfile.counts())
            numGenesWithOverlaps               = np.min(h5destination_integenic_overlap.nativeMeanProfile.counts())
            
            plot2WayMFEComparison( taxid,
                                   profileId,
                                   pd.DataFrame({'native':   h5destination_integenic_short.getNativeProfile(),
                                                 'shuffled': h5destination_integenic_short.getShuffledProfile()},
                                                index=xRange),
                                   pd.DataFrame({'native':   h5destination_integenic_overlap.getNativeProfile(),
                                                 'shuffled': h5destination_integenic_overlap.getShuffledProfile()},
                                                index=xRange),
                                   computationTag=args.computation_tag,
                                   yRange=(-13,1),
                                   comment="Intergenic region lengths\nShort (N>={}) Overlapping (N>={}) Total (N>={}) ".format( numGenesWithShortIntergenicRegions, numGenesWithOverlaps, numGenesWithShortIntergenicRegions + numGenesWithOverlaps ) )



            scatterPlot( taxid, profileId, dLFEvsFlankingLength, xvar='flanking_length', yvar='dLFE', colorvar='nextCDSOppositeStrand', title="flanking_length - %s")

            #plotMFEProfileV3(taxid, profileId, df, dLFEData=meanDeltaLFE, ProfilesCount=minCount)
            #plotMFEProfileV3(taxid, profileId, df, dLFEData=meanDeltaLFE, wilcoxon=wilcoxonDf, transitionPeak=guPeakDf, transitionPeakPos=peakPos*10, edgeWilcoxon=edgeWilcoxonDf, ProfilesCount=minCount)

            # Plot the number of genes included in each profile position
            plotXY(
                taxid,
                profileId,
                pd.DataFrame( { "num_genes": h5destination.nativeMeanProfile.counts() }, index=xRange ),
                "position",
                "num_genes",
                "Number of genes included, per starting position (all)",
                computationTag=args.computation_tag
                )
            plotXY(
                taxid,
                profileId,
                pd.DataFrame( { "num_genes_opposite": h5destination_next_cds_opposite.nativeMeanProfile.counts() }, index=xRange ),
                "position",
                "num_genes_opposite",
                "Number of genes included, per starting position (opposite)",
                computationTag=args.computation_tag
                )
            
            # scatterPlotWithKernel(
            #     taxid,
            #     profileId,
            #     geneLevelScatter,
            #     "gc",
            #     "logpval",
            #     "GC vs. MFE selection - %s"
            #     )

            # if( args.codonw ):
            #     scatterPlot(
            #         taxid,
            #         profileId,
            #         geneLevelScatter,
            #         "GC3s",
            #         "logpval",
            #         "GC3s vs. MFE selection - %s"
            #     )

            # if( args.codonw ):
            #     scatterPlot(
            #         taxid,
            #         profileId,
            #         geneLevelScatter,
            #         "gc",
            #         "Nc",
            #         "GC vs. ENc - %s"
            #     )

            # if( args.codonw ):
            #     scatterPlot(
            #         taxid,
            #         profileId,
            #         geneLevelScatter,
            #         "GC3s",
            #         "Nc",
            #         "GC3s vs. ENc - %s"
            #     )

            # if( args.codonw ):
            #     scatterPlot(
            #         taxid,
            #         profileId,
            #         geneLevelScatter,
            #         "Nc",
            #         "logpval",
            #         "ENc vs. MFE selection - %s"
            #     )

            # if( args.codonw ):
            #     scatterPlot(
            #         taxid,
            #         profileId,
            #         geneLevelScatter,
            #         "CBI",
            #         "logpval",
            #         "CBI vs. MFE selection - %s"
            #     )


            # if( args.pax_db ):
            #     #print(geneLevelScatter.head())
            #     scatterPlotWithColor(
            #         taxid,
            #         profileId,
            #         shuffleType,
            #         geneLevelScatter,
            #         "gc",
            #         "logpval",
            #         "pa",
            #         "GC vs. PA - %s"
            #     )

            #     if( args.codonw ):
            #         scatterPlot(
            #             taxid,
            #             profileId,
            #             geneLevelScatter,
            #             "Nc",
            #             "pa",
            #             "ENc vs. PA - %s"
            #         )


                # dfProfileByPA = pd.DataFrame( { "native": nativeMeanProfile.value(), "shuffled": shuffledMeanProfile.value(), "position": xRange, "shuffled25":shuffled25Profile.value(), "shuffled75":shuffled75Profile.value(), "native_pa_high":nativeMeanProfile_HighPAOnly.value(), "native_pa_med":nativeMeanProfile_MediumPAOnly.value(), "native_pa_low":nativeMeanProfile_LowPAOnly.value() }, index=xRange )

                # plotMFEProfileByPA(taxid, profileId, dfProfileByPA)

            # # Try to fit a linear model to describe the gene-level data
            # if( args.codonw ):
            #     if( args.pax_db ):
            #         model = ols("logpval ~ gc + cds_length_nt + Nc + GC3s + CAI + pa", data=geneLevelScatter).fit()
            #     else:
            #         model = ols("logpval ~ gc + cds_length_nt + Nc + GC3s + CAI", data=geneLevelScatter).fit()
            # else:
            #     model = ols("logpval ~ gc + cds_length_nt", data=geneLevelScatter).fit()

            # print(model.params)
            # print(model.summary())
            # print("r     = %f" % model.rsquared**.5)
            # print("r_adj = %f" % model.rsquared_adj**.5)



            spearman_rho = geneLevelScatter.corr(method='spearman')
            print(spearman_rho)
            spearman_rho.to_csv('mfe_v2_spearman_%d_%s_t%d.csv' % (taxid, profileId, shuffleType))



            # vars = ['gc', 'logpval', 'pa', 'cds_length_nt']
            # spearman_rho  = np.zeros((len(vars),len(vars)), dtype=float)
            # spearman_pval = np.zeros((len(vars),len(vars)), dtype=float)
            # for n1, var1 in enumerate(vars):
            #     for n2, var2 in enumerate(vars):
            #         rho, pval = spearmanr(geneLevelScatter[var1], geneLevelScatter[var2], nan_policy='omit')
            #         spearman_rho[n1,n2] = rho
            #         spearman_pval[n1,n2] = pval
            # print(spearman_rho)
            # print(spearman_pval)



            print(statisticsDF)


            # ------------------------------------------------------------------------------------
            # Print final report

            print("Got %d results" % n)

            print(x1)
            print(x2)
            print(x3)

        print("//"*20)
        print(list(combinedData.keys()))

        print("Found stop codons: {}".format(stopCodonFreq))

        if len(combinedData)>1:
            profileId = "%d_%d_%s" % (args.profile[0], args.profile[1], args.profile[2])
            plotMFEProfileForMultipleRandomizations(taxid, profileId, combinedData)

        return (taxid, (x1,x2,x3))


def calcProfilesForSpeciesX(taxid, args):

    # TEST ONLY ### TEST ONLY ### TEST ONLY ### TEST ONLY ### TEST ONLY ### TEST ONLY ### TEST ONLY ### TEST ONLY #
    #return (taxid,)
    # TEST ONLY ### TEST ONLY ### TEST ONLY ### TEST ONLY ### TEST ONLY ### TEST ONLY ### TEST ONLY ### TEST ONLY #
    
    p = ProfilePlot(taxid, args)
    return p.performPlots()

def runDistributed(args):
    import _distributed
    import dask

    scheduler = _distributed.open()

    results = {}

    taxids = []
    delayedCalls = []

    for taxid in args.taxid:
        call = dask.delayed( calcProfilesForSpeciesX )(taxid, args)
        delayedCalls.append( call )
        taxids.append(taxid)

    futures = scheduler.compute(delayedCalls) # submit all delayed calculations; obtain futures immediately

    try:
        _distributed.progress(futures) # wait for all calculations to complete
    except Exception as e:
        print(E)
    print("\n")

    print("Waiting for all tasks to complete...")
    _distributed.wait(futures)

    results = {}
    errorsCount = 0
    for taxid, f in zip(taxids, futures):
        try:
            r = scheduler.gather(f)
            returnedTaxId = r[0]
            assert(taxid==returnedTaxId)
            results[taxid] = r
            
        except Exception as e:
            print(e)
            results[taxid] = None
            errorsCount += 1

    print("Finished with %d errors" % errorsCount)
    return results
        

if __name__=="__main__":
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument("--taxid", type=parseList(int))
    argsParser.add_argument("--all-taxa", type=bool, default=False)
    argsParser.add_argument("--profile", type=parseProfileSpec())
    argsParser.add_argument("--computation-tag", type=int, default=Sources.RNAfoldEnergy_SlidingWindow40_v2)
    argsParser.add_argument("--shuffle-types", type=parseList(int) )
    argsParser.add_argument("--num-shuffles", type=int, default=20)
    argsParser.add_argument("--pax-db", type=str, required=False)
    argsParser.add_argument("--codonw", action="store_true", default=False)
    argsParser.add_argument("--distributed", action="store_true", default=False)
    args = argsParser.parse_args()

    if( args.all_taxa ):
        args.taxid = list(allSpeciesSource())

    if( args.taxid is None ):
        raise Exception("No species requested (use '--taxid tax1,tax2,tax3' or '--all-taxa')")
    
    # ------------------------------------------------------------------------------------
    # Argument validity checks
    if( len(args.taxid) > len(frozenset(args.taxid)) ):
        raise Exception("Duplicate taxid encountered in list %s" % args.taxid)  # Make sure no taxid was specified twice (will skew calculations...)

    checkSpeciesExist(args.taxid)  # Check for non-existant taxids to avoid doomed runs

    if( args.profile[2] != "begin" and args.profile[2] != "end" and args.profile[2] != "stop3utr" ):
        raise Exception("Unsupported profile reference '%s'" % args.profile[2]) # Currently only profile with reference to CDS 'begin' are implemented...

    results = None
    if( not args.distributed ):
        # run locally
        results = {}
        for taxid in args.taxid:
            ret = calcProfilesForSpeciesX(taxid, args)
            results[taxid] = ret
    else:
        results = runDistributed(args)

    print(results)
    print("Total results: %d" % len(results))
    print("Succeeded: %d" % len([1 for x in list(results.values()) if not x is None]))
    failed = [k for k,v in list(results.items()) if v is None]
    print("Failed: %d (%s)" % (len(failed), failed))
    

            
    
