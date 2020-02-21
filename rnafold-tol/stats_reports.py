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
# Create report tables summarizing the data (and adding stats useful for evaluation data validity and coverage)
import pandas as pd
import dask
from random import randint
from Bio.Seq import Seq
from tabulate import tabulate
from ncbi_phyla_report import parseReport
from data_helpers import allSpeciesSource, getSpeciesName, countSpeciesCDS, getSpeciesProperty, nativeSequencesSource, getSpeciesTranslationTable
from mfe_plots import getSpeciesShortestUniqueNamesMapping_memoized
from process_series_data import readSeriesResultsForSpecies, readSeriesResultsForSpeciesWithSequence, convertResultsToMFEProfiles, sampleProfilesFixedIntervals, profileLength, profileElements, MeanProfile, calcSampledGCcontent
from collections import Counter
from ncbi_taxa import ncbiTaxa

# Configuration
numShuffledGroups = 20

#speciesToExclude = frozenset((405948,999415,946362,470, 1280, 4932, 508771, 2850, 753081, 195065, 641309 ))
speciesToExclude = frozenset(( 195065, 641309 ))
# 158189,456481,272632,1307761,505682

shortNames = getSpeciesShortestUniqueNamesMapping_memoized()

 
def formatLineage(lineage, names):
    out = []
    for l in lineage:
        name = names[l]
        out.append(name)
    return ", ".join(out)
    

"""
Generate phylum-level species count report
This may be needed (as a supplementary table) for the paper, and is also useful for planning dataset coverage

To display the reference list of taxons of a given rank, use the following type of Entrez query:
Firmicutes[Subtree] AND Order[Rank] 

"""
def speciesByPhylaTable():
    allPhyla = parseReport() # get all existing phyla

    domainCounts = Counter()
    phylaCounts = Counter()
    skippedCounts = Counter()
    #classesByPhyla = {}   # Disable tallying by class, since these are not used for many taxons
    ordersByPhyla = {}
    familiesByPhyla = {}
    genusesByPhyla = {}


    phylaDf = pd.DataFrame({'Domain': pd.Categorical([]),                # Bacteria, Eukaryota, Archaea
                            'Phylum': pd.Categorical([]),                # Phylum name (string)
                            'TaxId': pd.Series([], dtype='int'),         # Phylum TaxId
                            'ParentTaxId': pd.Series([], dtype='int'),   # Parent TaxId
                            'NumSpecies': pd.Series([], dtype='int'),    # Species count for this phyla
     #                       'NumClasses': pd.Series([], dtype='int'),    # Species count for this phyla
                            'NumOrders': pd.Series([], dtype='int'),     # Orders count for this phyla
                            'NumFamilies': pd.Series([], dtype='int'),   # Families count for this phyla
                            'NumGenuses': pd.Series([], dtype='int'),    # Genuses count for this phyla
                            'RowType': pd.Categorical([])  })            # Phylum count or total

    for group, phyla in allPhyla.items():
        for phylum, record in phyla.items():
            # Add item for each phylum
            taxId = record['taxId']
            
            phylaDf = phylaDf.append(pd.DataFrame({'Domain': pd.Categorical([group]),
                                                   'Phylum': pd.Categorical([phylum]),
                                                   'TaxId': pd.Series([taxId], dtype='int'),
                                                   'ParentTaxId': pd.Series([record['parentTaxId']], dtype='int'),
                                                   'NumSpecies': pd.Series([0], dtype='int'),
            #                                       'NumClasses': pd.Series([0], dtype='int'),
                                                   'NumOrders': pd.Series([0], dtype='int'),
                                                   'NumFamilies': pd.Series([0], dtype='int'),
                                                   'NumGenuses': pd.Series([0], dtype='int'),
                                                   'RowType': pd.Categorical(['Phylum']) } ))
            #classesByPhyla[record['taxId']]  = set()
            ordersByPhyla[record['taxId']]   = set()
            familiesByPhyla[record['taxId']] = set()
            genusesByPhyla[record['taxId']]  = set()


    # Create "special" items
    pid = 1
    for group in allPhyla.keys():
        # Add "Unknown phylum" tally for each domain
        phylaDf = phylaDf.append(pd.DataFrame({'Domain': pd.Categorical([group]),
                                               'Phylum': pd.Categorical(['[Unknown]']),
                                               'TaxId': pd.Series([pid], dtype='int'),
                                               'ParentTaxId': pd.Series([0], dtype='int'),
                                               'NumSpecies': pd.Series([0], dtype='int'),
        #                                       'NumClasses': pd.Series([0], dtype='int'),
                                               'NumOrders': pd.Series([0], dtype='int'),
                                               'NumFamilies': pd.Series([0], dtype='int'),
                                               'NumGenuses': pd.Series([0], dtype='int'),
                                               'RowType': pd.Categorical(['Total']) }))
        pid += 1
        # Add totals tally for each domain
        phylaDf = phylaDf.append(pd.DataFrame({'Domain': pd.Categorical([group]),
                                               'Phylum': pd.Categorical(['[Total]']),
                                               'TaxId': pd.Series([pid], dtype='int'),
                                               'ParentTaxId': pd.Series([0], dtype='int'),
                                               'NumSpecies': pd.Series([0], dtype='int'),
        #                                       'NumClasses': pd.Series([0], dtype='int'),
                                               'NumOrders': pd.Series([0], dtype='int'),
                                               'NumFamilies': pd.Series([0], dtype='int'),
                                               'NumGenuses': pd.Series([0], dtype='int'),
                                               'RowType': pd.Categorical(['Total']) }))
        pid += 1
    # Add overally totals items
    phylaDf = phylaDf.append(pd.DataFrame({'Domain': pd.Categorical(['[All]']),
                                           'Phylum': pd.Categorical(['[Total]']),
                                           'TaxId': pd.Series([pid], dtype='int'),
                                           'ParentTaxId': pd.Series([0], dtype='int'),
                                           'NumSpecies': pd.Series([0], dtype='int'),
    #                                       'NumClasses': pd.Series([0], dtype='int'),
                                           'NumOrders': pd.Series([0], dtype='int'),
                                           'NumFamilies': pd.Series([0], dtype='int'),
                                           'NumGenuses': pd.Series([0], dtype='int'),
                                           'RowType': pd.Categorical(['Total']) }))
            
    phylaDf.set_index('TaxId', inplace=True)
    skippedSpecies = []

    # Count the number of species under each phylum
    for taxId in allSpeciesSource():
        if taxId in speciesToExclude: continue
        lineage = ncbiTaxa.get_lineage(taxId)
        names = ncbiTaxa.get_taxid_translator(lineage)

        ranks = ncbiTaxa.get_rank(lineage)

        # Determine kingdom/domain
        kingdomTaxId = [t for t,rank in ranks.items() if rank=='superkingdom']
        if not kingdomTaxId:
            kingdomTaxId = [t for t,rank in ranks.items() if rank=='kingdom']
        domain = names[kingdomTaxId[0]]
        domainCounts.update([domain])
        
        # Determine phylum
        phylumTaxId = [t for t,rank in ranks.items() if rank=='phylum']
        if not phylumTaxId:
            skippedSpecies.append(taxId)
            skippedCounts.update([domain])
            print("Skipping %d: (%s) missing phylum" % (taxId, names[taxId]))
            #print(formatLineage(lineage, names))
            continue  # This table is structured by phylum; information will be missing for any species missing a phylum; it will be included in the "species missing phylum" ([Unknown]) row.
        else:
            phylumTaxId = phylumTaxId[0]

        if phylumTaxId:
            phylaCounts.update([phylumTaxId]) # tally this species under the specified phylum
            
        #classTaxId = [t for t,rank in ranks.items() if rank=='class']
        #if classTaxId:
        #    classesByPhyla[phylumTaxId].add( classTaxId[0] )

        orderTaxId = [t for t,rank in ranks.items() if rank=='order']
        if orderTaxId:
            ordersByPhyla[phylumTaxId].add( orderTaxId[0] )

        familyTaxId = [t for t,rank in ranks.items() if rank=='family']
        if familyTaxId:
            familiesByPhyla[phylumTaxId].add( familyTaxId[0] )

        genusTaxId = [t for t,rank in ranks.items() if rank=='genus']
        if genusTaxId:
            genusesByPhyla[phylumTaxId].add( genusTaxId[0] )
            


    assert(sum(skippedCounts.values()) == len(skippedSpecies))


    # Update the phyla counts
    for phylaTaxId, counts in phylaCounts.items():
        #phylaDf.loc[phylaTaxId, 'NumClasses']  = len(classesByPhyla[phylaTaxId])
        phylaDf.loc[phylaTaxId, 'NumOrders']   = len(ordersByPhyla[phylaTaxId])
        phylaDf.loc[phylaTaxId, 'NumFamilies'] = len(familiesByPhyla[phylaTaxId])
        phylaDf.loc[phylaTaxId, 'NumGenuses']  = len(genusesByPhyla[phylaTaxId])
        phylaDf.loc[phylaTaxId, 'NumSpecies']  = counts


    # Update the "Unknown phyla" count for each domain
    for group, countMissing in skippedCounts.items():
        #print('-'*20)
        #print("%s - %d missing" % (group, countMissing))
        dummyTaxIdForBasketGroup = phylaDf[ (phylaDf.Domain==group) & (phylaDf.Phylum=='[Unknown]') ].index[0]
        phylaDf.loc[dummyTaxIdForBasketGroup, 'NumSpecies'] = countMissing

    # Update the total for each domain
    for group, totalCount in domainCounts.items():
        dummyTaxIdForBasketGroup = phylaDf[ (phylaDf.Domain==group) & (phylaDf.Phylum=='[Total]') ].index[0]
        phylaDf.loc[dummyTaxIdForBasketGroup, 'NumSpecies'] = totalCount

    # Update the overall total count
    dummyTaxIdForBasketGroup = phylaDf[ (phylaDf.Domain=="[All]") & (phylaDf.Phylum=='[Total]') ].index[0]
    phylaDf.loc[dummyTaxIdForBasketGroup, 'NumSpecies']  = sum(domainCounts.values())
    phylaDf.loc[dummyTaxIdForBasketGroup, 'NumOrders']   = sum([len(x) for x in ordersByPhyla.values()])
    phylaDf.loc[dummyTaxIdForBasketGroup, 'NumFamilies'] = sum([len(x) for x in familiesByPhyla.values()])
    phylaDf.loc[dummyTaxIdForBasketGroup, 'NumGenuses']  = sum([len(x) for x in genusesByPhyla.values()])


    # Prepare and save the final table
    phylaReportDf = phylaDf[phylaDf['NumSpecies'] > 0]              # remove "empty" items
    phylaReportDf = phylaReportDf.sort_values(by=['Domain', 'RowType', 'Phylum' ])    # sort rows
    print(phylaReportDf)
    phylaReportDf.to_html( 'phyla_report.html', columns=['Phylum', 'NumOrders', 'NumFamilies', 'NumGenuses', 'NumSpecies', 'Domain'])
    phylaReportDf.to_excel('phyla_report.xlsx', sheet_name='Phyla Summary')

    with open("phyla_report.rst", "w") as f:
        f.write( phylaReportDf.drop(['RowType', 'NumFamilies', 'NumGenuses', 'NumOrders', 'ParentTaxId'], axis=1).pipe( tabulate, headers='keys', tablefmt='rst' ) )
    

    # Prepare the "Missing phyla" report
    missingPhylaReportDf = phylaDf[phylaDf['NumSpecies'] == 0]
    missingPhylaReportDf = missingPhylaReportDf.sort_values(by=['Domain', 'RowType', 'Phylum' ])    # sort rows
    missingPhylaReportDf.to_html( 'phyla_report_missing.html', columns=['Phylum', 'NumSpecies', 'Domain'])
    missingPhylaReportDf.to_excel('phyla_report_missing.xlsx', sheet_name='Missing Phyla Summary')

    # print counts
    print(domainCounts)
    #print(phylaCounts)


    # Display "skipped items" warning
    if( skippedSpecies):
        print("="*50)
        print("Warning: Skipped %d species" % len(skippedSpecies))
        print(skippedCounts)
        print("="*50)


        

def summarizeCounter(counter):
    ret = []
    for item, count in counter.items():
        ret.append("{}:{}".format( item, count))
    return ','.join(ret)

def countShuffledProfiles(taxId, profile, computationTag, shuffleType):
    
    shuffledMeanProfile = MeanProfile( profileLength(profile) )

    for result in sampleProfilesFixedIntervals(
            convertResultsToMFEProfiles(
                readSeriesResultsForSpeciesWithSequence((computationTag,), taxId, numShuffledGroups, numShuffledGroups, shuffleType=shuffleType )
                , numShuffledGroups)
            , profile[3], profile[0], profile[1], profile[2]):
        
        profileData = result["profile-data"]
        
        shuffledMeanProfile.add( profileData[1:] )

    print(shuffledMeanProfile.counts())

    numShuffledSeqs =  shuffledMeanProfile.counts()[0] / numShuffledGroups
    return (taxId, numShuffledSeqs)

@dask.delayed
def calcNativeSequencesStatistics(taxId, fraction, numFractions):

    #countPairedNucleotides = 0
    #countTotalNucleotides  = 0
    cdsCount    = 0
    gcCount     = 0
    totalCount  = 0
    cdsWarnings = 0
    warnings = Counter()
    firstAA  = Counter()
    lastAA   = Counter()
    

    geneticCode = getSpeciesTranslationTable(taxId)
    
    
    for seqId, seq in nativeSequencesSource(taxId, fraction, numFractions):
        seq = seq.lower()
        seqHasWarnings = False

        gcCount    += sum([1 for x in seq if (x=='c' or x=='g')])
        totalCount += sum([1 for x in seq if (x=='c' or x=='g' or x=='a' or x=='t')])  # don't count 'N's

        if len(seq)%3 != 0:
            seqHasWarnings = True
            warnings['cds-length'] += 1

        xlation = Seq(seq).translate(table=geneticCode).lower()
        if xlation[0] != 'm':
            seqHasWarnings = True
            warnings['translation-methionine'] += 1

        if xlation[-1] != '*':
            seqHasWarnings = True
            warnings['translation-stop-codon'] += 1

        if seqHasWarnings:
            cdsWarnings += 1
            
        firstAA.update(xlation[0])
        lastAA.update(xlation[-1])
            
        cdsCount += 1

        
    #print("Total:  %d" % countTotalNucleotides)
    #print("Paired: %d (%.3g%%)" % (countPairedNucleotides, float(countPairedNucleotides)/countTotalNucleotides*100))
    

    return (taxId, fraction, cdsCount, gcCount, totalCount, cdsWarnings, warnings, firstAA, lastAA)

#@dask.delayed
#def countComputedResultsForSeries( ):
#    
#    for a in readSeriesResultsForSpecies( seriesSourceNumber, species, minShuffledGroups=20, maxShuffledGroups=20, shuffleType=db.Sources.ShuffleCDSv2_python, cdsFilter=None, returnCDS=True )
#
        

def speciesStatisticsAndValidityReport(args):
    import _distributed

    speciesDf = pd.DataFrame({
        'TaxId': pd.Series([], dtype='int'),           # Species TaxId
        'Species': pd.Series([], dtype='str'),         # Species binomial name
        'Nickname': pd.Series([], dtype='str'),
        'Domain': pd.Categorical([]),                  # Bacteria, Eukaryota, Archaea
        'Phylum': pd.Categorical([]),                  # Phylum name (string)
        'NumCDSs': pd.Series([], dtype='int'),         # CDS count for this species
        'NumCDSsInProfile': pd.Series([], dtype='int'),         # Num seqs with 20 shuffled profiles for this species
        'AnnotatedNumCDSs': pd.Series([], dtype='int'),   # 
        'CDSDifference': pd.Series([], dtype='float'), # 
        'NumNativeSeqs': pd.Series([], dtype='int'),      # 
        'GCContentInCDS': pd.Series([], dtype='float'), # 
        'AnnotatedGCContent': pd.Series([], dtype='float'), # 
        'RowType': pd.Categorical([]),                 # Species count or total
        'Warnings': pd.Series([], dtype='str'),        # 
        'CDSWarnings': pd.Series([],  dtype='int'),    # 
        'CDSWarnings_': pd.Series([], dtype='str'),    # 
        'FirstAA': pd.Series([], dtype='str'),         # 
        'LastAA' : pd.Series([], dtype='str')          # 
    })
    
    scheduler = _distributed.open()

    results = {}
    delayedCalls_native = []
    
    shuffledCounts = {}
    delayedCalls_shuffledProfiles  = []

    
    for taxId in allSpeciesSource():
        if taxId in speciesToExclude: continue # always exclude species from the blacklist
        if args.taxid and taxId not in args.taxid: continue  # if a whitelist is specified, skip other species

        warnings = []

        ## DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ##
        #if randint(0, 20) > 0:
        #    continue
        ## DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ##
   
        cdsCountInRedis = countSpeciesCDS(taxId)

        #cdsCountProfiles = countx(taxId, (310, 10, "begin", 0), 102, 11)
        

        annotatedProteinCount = getSpeciesProperty(taxId, 'protein-count')[0]

        annotatedGCContent = getSpeciesProperty(taxId, 'gc-content')[0]

        proteinDifference = None
        if not annotatedProteinCount is None:
            proteinDifference = (1.0 - float(cdsCountInRedis) / float(annotatedProteinCount)) * 100.0

            if abs(proteinDifference) > 9.9:
                warnings.append("CDS_count")
        else:
            warnings.append("No_CDS_count")


        # Determine phylum
        lineage = ncbiTaxa.get_lineage(taxId)
        names = ncbiTaxa.get_taxid_translator(lineage)

        ranks = ncbiTaxa.get_rank(lineage)

        # Determine kingdom/domain
        domain = ""
        kingdomTaxId = [t for t,rank in ranks.items() if rank=='superkingdom']
        if not kingdomTaxId:
            kingdomTaxId = [t for t,rank in ranks.items() if rank=='kingdom']
        domain = names[kingdomTaxId[0]]

        phylumName = ""
        # Determine phylum
        phylumTaxId = [t for t,rank in ranks.items() if rank=='phylum']
        if phylumTaxId:
            phylumName = names[phylumTaxId[0]]

        
                

        speciesDf = speciesDf.append(pd.DataFrame({
            'TaxId': pd.Series([taxId], dtype='int'),                    # Species TaxId
            'Species': pd.Series([getSpeciesName(taxId)], dtype='str'),
            'Nickname': pd.Series([shortNames[taxId]], dtype='str'),
            'Domain': pd.Categorical([domain]),                             # Bacteria, Eukaryota, Archaea
            'Phylum': pd.Categorical([phylumName]),                             # Phylum name (string)
            'NumCDSs': pd.Series([cdsCountInRedis], dtype='int'),        # CDS count for this species
            'NumCDSsInProfile': pd.Series([0], dtype='int'),          # Num seqs with 20 shuffled profiles
            'AnnotatedNumCDSs': pd.Series([0 if annotatedProteinCount is None else annotatedProteinCount], dtype='int'),   # 
            'CDSDifference': pd.Series([proteinDifference], dtype='float'), # 
            'NumNativeSeqs': pd.Series([0], dtype='int'),                   # 
            'GCContentInCDS': pd.Series([0.0], dtype='float'), # 
            'AnnotatedGCContent': pd.Series([annotatedGCContent], dtype='float'), # 
            'RowType': pd.Categorical(["species"]),                      # Species count or total
            'Warnings': pd.Series([", ".join(warnings)], dtype='str'),                       #
            'CDSWarnings': pd.Series([0], dtype='int'),
            'CDSWarnings_': pd.Series([""], dtype='str'),
            'FirstAA': pd.Series([""], dtype='str'),
            'LastAA' : pd.Series([""], dtype='str'),
            'Source':  pd.Series([""], dtype='str')
            }))

        fractionSize = 1000   # How many sequences (roughly) to process in each task
        numFractions = cdsCountInRedis/fractionSize
        if numFractions == 0: numFractions = 1
                
        for i in range(numFractions):
            # DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #
            #if i%100!=5: continue
            # DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #
            
            call = dask.delayed( calcNativeSequencesStatistics )(taxId, i, numFractions)
            delayedCalls_native.append( call )

        call = dask.delayed( countShuffledProfiles )(taxId, (310, 10, "begin", 0), 102, 11)
        delayedCalls_shuffledProfiles.append( call )
        

    speciesDf.set_index('TaxId', inplace=True)

    print("Starting {} calls...".format(len(delayedCalls_native)+len(delayedCalls_shuffledProfiles)) )

    futures = scheduler.compute(delayedCalls_native + delayedCalls_shuffledProfiles) # submit all delayed calculations; obtain futures immediately

    try:
        _distributed.progress(futures) # wait for all calculations to complete
    except Exception as e:
        print(E)
    print("\n")

    print("Waiting for all tasks to complete...")
    _distributed.wait(futures)

    results = {}

    errorsCount = 0
    for f in futures:
        try:
            ret = scheduler.gather(f)
            if( len(ret)==9 ):
                (taxId, fraction, cdsCount, gcCounts, totalCounts, cdsWarnings, warnings, firstAA, lastAA) = ret

                current = None
                if taxId in results:
                    current = results[taxId]
                else:
                    current = (0, 0, 0, 0, Counter(), Counter(), Counter())

                current = (current[0] + cdsCount,
                           current[1] + gcCounts,
                           current[2] + totalCounts,
                           current[3] + cdsWarnings,
                           current[4] + warnings,
                           current[5] + firstAA,
                           current[6] + lastAA)

                results[taxId] = current
                
            elif( len(ret)==2 ):
                (taxId, numShuffledSeqs) = ret
                shuffledCounts[taxId] = numShuffledSeqs
                
            else:
                assert(False)
                
            
        except Exception as e:
            print(e)
            errorsCount += 1

    for taxId, result in results.items():
        (numNativeSeqs, gcCounts, totalCounts, cdsWarnings, warnings, firstAA, lastAA) = result
        speciesDf.at[taxId, 'NumNativeSeqs'] = numNativeSeqs

        speciesDf.at[taxId, 'GCContentInCDS'] = round( float(gcCounts)/float(totalCounts)*100.0, 1)
        
        speciesDf.at[taxId, 'CDSWarnings'] = cdsWarnings
        
        speciesDf.at[taxId, 'CDSWarnings_'] = summarizeCounter(warnings)
        speciesDf.at[taxId, 'FirstAA']      = summarizeCounter(firstAA)
        speciesDf.at[taxId, 'LastAA']       = summarizeCounter(lastAA)
        
        #if numNativeSeqs < species.at[taxId, 'NumCDSs']:
        #    pass

    for taxId, result in shuffledCounts.items():
        speciesDf.at[taxId, 'NumCDSsInProfile']       = result

    

    speciesDf = speciesDf.sort_values(by=['Domain', 'Species' ])    # sort rows
    speciesDf.to_html( 'species_report.html', float_format='{0:.1f}'.format, columns=['Species', 'Nickname', 'NumCDSs', 'NumCDSsInProfile', 'AnnotatedNumCDSs', 'CDSDifference', 'NumNativeSeqs', 'GCContentInCDS', 'AnnotatedGCContent', 'Phylum', 'Domain', 'Warnings', 'CDSWarnings', 'CDSWarnings_', 'FirstAA', 'LastAA'])

    with open("species_report_simple.rst", "w") as f:
        f.write( speciesDf.drop(['RowType', 'Warnings', 'CDSWarnings', 'CDSWarnings_', 'FirstAA', 'LastAA', 'CDSDifference'], axis=1).pipe( tabulate, headers='keys', tablefmt='rst' ) )
        
    speciesDf.to_html( 'species_report_simple.html', float_format='{0:.1f}'.format, columns=['Species', 'Nickname', 'NumCDSs', 'NumCDSsInProfile', 'AnnotatedNumCDSs', 'CDSDifference', 'NumNativeSeqs', 'GCContentInCDS', 'AnnotatedGCContent', 'Phylum', 'Domain'])

    
    
    speciesDf.to_excel('species_report.xlsx', sheet_name='Species summary')
    
    
def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert


def standalone():
    import argparse
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument("--taxid", type=parseList(int))
    args = argsParser.parse_args()
    
    speciesByPhylaTable()
    speciesStatisticsAndValidityReport(args)
    return 0
    

if __name__=="__main__":
    import sys
    sys.exit(standalone())

