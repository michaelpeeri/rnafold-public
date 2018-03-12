# Replaces convert_data_for_plotting.py
#
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
from data_helpers import decompressNucleicSequence, checkSpeciesExist, getSpeciesFileName, allSpeciesSource
from mysql_rnafold import Sources
from process_series_data import readSeriesResultsForSpeciesWithSequence, convertResultsToMFEProfiles, sampleProfilesFixedIntervals, profileLength, profileElements, MeanProfile, calcSampledGCcontent
from mfe_plots import plotMFEProfileWithGC, plotXY, scatterPlot, scatterPlotWithColor, plotMFEProfileByPA, plotMFEProfileMultiple, scatterPlotWithKernel
from codonw import readCodonw, meanCodonwProfile


# ------------------------------------------------------------------------------------
# Configuration
confWindowWidth = 40


# ------------------------------------------------------------------------------------
# Command-line args

def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert
    
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
                    

    def performPlots(self):
        
        args = self._args
        taxid = self._taxId

        # ------------------------------------------------------------------------------------

        x1 = Counter()
        x2 = Counter()
        x3 = Counter()

        numShuffledGroups = args.num_shuffles

        n = 0

        print("Processing species %d (%s)" % (taxid, getSpeciesFileName(taxid)))

        nativeMeanProfile = MeanProfile( profileLength(args.profile) )
        shuffledMeanProfile = MeanProfile( profileLength(args.profile) )

        shuffled25Profile = MeanProfile( profileLength(args.profile) )
        shuffled75Profile = MeanProfile( profileLength(args.profile) )

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

        geneLevelScatter = pd.DataFrame({'gc':pd.Series(dtype='float'), 'logpval':pd.Series(dtype='float'), 'abslogpval':pd.Series(dtype='float'), 'protid':pd.Series(dtype='string')})

        cdsLengths = []

        fullSeqs = []
        dfCodonw = None

        # ------------------------------------
        # Process all CDS for this species
        # ------------------------------------
        for result in sampleProfilesFixedIntervals(
                convertResultsToMFEProfiles(
                    readSeriesResultsForSpeciesWithSequence(args.computation_tags, taxid, numShuffledGroups, numShuffledGroups )
                    , numShuffledGroups)
                , args.profile[3], args.profile[0], args.profile[1]):

            fullCDS = result["cds-seq"]
            seq = fullCDS[args.profile[3]:args.profile[0]]

            if not seq:
                continue

            protId = result["cds"].getProtId()

            fullSeqs.append( SeqRecord( Seq(fullCDS, NucleotideAlphabet), id=protId) )


            profileData = result["profile-data"]
            assert(profileData.shape[0] >= numShuffledGroups)

            # Prepare mean MFE profiles
            nativeMeanProfile.add( profileData[0,None] )
            shuffledMeanProfile.add( profileData[1:] )

            # Prepare GC profile
            gc = calcSampledGCcontent( seq, args.profile[1] )
            if( gc.size > profileLength(args.profile) ):  # truncate the profile if necessary
                gc = np.resize( gc, (profileLength(args.profile),))
            GCProfile.add( np.expand_dims(gc,0) )

            # Prepare percentile mean profiles
            shuffled25Profile.add( np.expand_dims( np.percentile(profileData[1:], 25, axis=0), 0) )
            shuffled75Profile.add( np.expand_dims( np.percentile(profileData[1:], 75, axis=0), 0) )

            # Prepare data for genome-wide wilcoxon test
            #newDeltas = profileData[0,0::4] - np.mean(profileData[1:,0::4], axis=0)
            newDeltas = profileData[0,0::1] - np.mean(profileData[1:,0::1], axis=0)
            #newPositions = range(args.profile[3], profileLength(args.profile), 40)
            newPositions = range(args.profile[3], args.profile[0], args.profile[1])
            deltaspd = pd.DataFrame({'pos':pd.Series(newPositions, dtype='int'), 'delta':pd.Series(newDeltas, dtype='float')})
            deltasForWilcoxon = deltasForWilcoxon.append(deltaspd)

            # Prepare data for GC vs. selection test
            meanGC = calcSampledGCcontent( seq, 10000)[0]
            if( not (meanGC >= 0.05 and meanGC <= 0.95)):
                meanGC = None
            #deltas = profileData[0,0::4] - np.mean(profileData[1:,0::4], axis=0)
            deltas = profileData[0,0::1] - np.mean(profileData[1:,0::1], axis=0)
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


            geneLevelScatter = geneLevelScatter.append(pd.DataFrame({'gc':pd.Series([meanGC]), 'logpval': pd.Series([directedLogPval]), 'abslogpval': pd.Series([pvalue]), 'protid':pd.Series([protId]), 'pa':pd.Series([paval]), 'cds_length_nt':pd.Series([cds_length_nt])}))


            x1.update((fullCDS[0],))
            x2.update((fullCDS[1],))
            x3.update((fullCDS[2],))
            del fullCDS

            del result
            n += 1

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
            CUBmetricsProfile = meanCodonwProfile(fullSeqs, confWindowWidth, 'begin', args.profile[3], args.profile[0], args.profile[1]) # TODO - use real values!
            print(CUBmetricsProfile)

        #else:
        #    geneLevelScatter['CAI'] = pd.Series( np.zeros(len(geneLevelScatter)), index=geneLevelScatter.index)



        # ------------------------------------
        # Display summary for this species
        # ------------------------------------
        print("Native:")
        print(nativeMeanProfile.value())
        print(nativeMeanProfile.counts())

        print("Shuffled:")
        print(shuffledMeanProfile.value())
        print(shuffledMeanProfile.counts())

        print(deltasForWilcoxon.shape)


        #native = np.asarray(nativeMean[1:], dtype="float")
        #shuffled = np.asarray(shuffledMean[1:], dtype="float")
        #gc = np.asarray(gcMean[1:], dtype="float")
        #xrange = [x for x in args.profile.Elements() if x<profileInfo.cdsLength()]
        xRange = profileElements(args.profile)
        profileId = "%d_%d_%s" % (args.profile[0], args.profile[1], args.profile[2])


        df = pd.DataFrame( { "native": nativeMeanProfile.value(), "shuffled": shuffledMeanProfile.value(), "gc":GCProfile.value(), "position": xRange, "shuffled25":shuffled25Profile.value(), "shuffled75":shuffled75Profile.value()}, index=xRange )
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        if( not CUBmetricsProfile is None ):
            df = pd.merge( df, CUBmetricsProfile, how='left', left_on='position', right_on='windowStart')
            df = df.set_index('position')

        #print(df)

        dfProfileCorrs = None
        if( args.codonw ):
            plotMFEProfileMultiple(taxid, profileId, df, ('GC', 'Nc', 'CAI', 'CBI', 'Fop'), scaleBar=confWindowWidth)

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


        plotMFEProfileWithGC(taxid, profileId, df)

        # Plot the number of genes included in each profile position
        plotXY(
            taxid,
            profileId,
            pd.DataFrame( { "num_genes": nativeMeanProfile.counts() }, index=xRange ),
            "position",
            "num_genes",
            "Number of genes included, per starting position"
            )

        scatterPlotWithKernel(
            taxid,
            profileId,
            geneLevelScatter,
            "gc",
            "logpval",
            "GC vs. MFE selection - %s"
            )

        if( args.codonw ):
            scatterPlot(
                taxid,
                profileId,
                geneLevelScatter,
                "GC3s",
                "logpval",
                "GC3s vs. MFE selection - %s"
            )

        if( args.codonw ):
            scatterPlot(
                taxid,
                profileId,
                geneLevelScatter,
                "gc",
                "Nc",
                "GC vs. ENc - %s"
            )

        if( args.codonw ):
            scatterPlot(
                taxid,
                profileId,
                geneLevelScatter,
                "GC3s",
                "Nc",
                "GC3s vs. ENc - %s"
            )

        if( args.codonw ):
            scatterPlot(
                taxid,
                profileId,
                geneLevelScatter,
                "Nc",
                "logpval",
                "ENc vs. MFE selection - %s"
            )

        if( args.codonw ):
            scatterPlot(
                taxid,
                profileId,
                geneLevelScatter,
                "CBI",
                "logpval",
                "CBI vs. MFE selection - %s"
            )


        if( args.pax_db ):
            #print(geneLevelScatter.head())
            scatterPlotWithColor(
                taxid,
                profileId,
                geneLevelScatter,
                "gc",
                "logpval",
                "pa",
                "GC vs. PA - %s"
            )

            if( args.codonw ):
                scatterPlot(
                    taxid,
                    profileId,
                    geneLevelScatter,
                    "Nc",
                    "pa",
                    "ENc vs. PA - %s"
                )


            dfProfileByPA = pd.DataFrame( { "native": nativeMeanProfile.value(), "shuffled": shuffledMeanProfile.value(), "position": xRange, "shuffled25":shuffled25Profile.value(), "shuffled75":shuffled75Profile.value(), "native_pa_high":nativeMeanProfile_HighPAOnly.value(), "native_pa_med":nativeMeanProfile_MediumPAOnly.value(), "native_pa_low":nativeMeanProfile_LowPAOnly.value() }, index=xRange )

            plotMFEProfileByPA(taxid, profileId, dfProfileByPA)

        # Try to fit a linear model to describe the gene-level data
        if( args.codonw ):
            if( args.pax_db ):
                model = ols("logpval ~ gc + cds_length_nt + Nc + GC3s + CAI + pa", data=geneLevelScatter).fit()
            else:
                model = ols("logpval ~ gc + cds_length_nt + Nc + GC3s + CAI", data=geneLevelScatter).fit()
        else:
            model = ols("logpval ~ gc + cds_length_nt", data=geneLevelScatter).fit()

        print(model.params)
        print(model.summary())
        print("r     = %f" % model.rsquared**.5)
        print("r_adj = %f" % model.rsquared_adj**.5)



        spearman_rho = geneLevelScatter.corr(method='spearman')
        print(spearman_rho)
        spearman_rho.to_csv('mfe_v2_spearman_%d_%s.csv' % (taxid, profileId))



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


        # -----------------------------------------------------------------------------
        # Save mean profiles as H5

        # Format (for compatible with plot_xy.py and old convert_data_for_plotting.py:
        #         gc  native  position  shuffled
        # 1    0.451  -4.944         1    -5.886
        # 2    0.459  -5.137         2    -6.069
        # 3    0.473  -5.349         3    -6.262


        # Compression parameters are described here:  http://www.pytables.org/usersguide/libref/helper_classes.html#filtersclassdescr
        # ...and discussed thoroughly in the performance FAQs
        with pd.io.pytables.HDFStore("gcdata_v2_taxid_%d_profile_%d_%d_%s_%d.h5" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3]), complib="zlib", complevel=1) as store:
            store["df_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = df
            store["deltas_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = deltasForWilcoxon
            store["spearman_rho_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = spearman_rho
            store["statistics_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = statisticsDF
            if( args.codonw ):
                store["profiles_spearman_rho_%d_%d_%d_%s_%d" % (taxid, args.profile[0], args.profile[1], args.profile[2], args.profile[3])] = dfProfileCorrs

            store.flush()


        # ------------------------------------------------------------------------------------
        # Print final report

        print("Got %d results" % n)

        print(x1)
        print(x2)
        print(x3)

        return (taxid, (x1,x2,x3))


def calcProfilesForSpecies(taxid, args):
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
        call = dask.delayed( calcProfilesForSpecies )(taxid, args)
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
    argsParser.add_argument("--computation-tags", type=int, default=(Sources.RNAfoldEnergy_SlidingWindow40_v2,))
    argsParser.add_argument("--num-shuffles", type=int, default=20)
    argsParser.add_argument("--pax-db", type=str, required=False)
    argsParser.add_argument("--codonw", type=bool, default=False)
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

    if( args.profile[2] != "begin"):
        raise Exception("Unsupported profile reference '%s'" % args.profile[2]) # Currently only profile with reference to CDS 'begin' are implemented...

    results = None
    if( not args.distributed ):
        # run locally
        results = {}
        for taxid in args.taxid:
            ret = calcProfilesForSpecies(taxid, args)
            results[taxid] = ret
    else:
        results = runDistributed(args)

    print(results)
    print("Total results: %d" % len(results))
    print("Succeeded: %d" % len([1 for x in results.values() if not x is None]))
    failed = [k for k,v in results.items() if v is None]
    print("Failed: %d (%s)" % (len(failed), failed))
    

            
    
