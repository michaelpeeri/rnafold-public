# Create report tables summarizing the data
from ncbi_phyla_report import parseReport
from data_helpers import allSpeciesSource
from collections import Counter
import pandas as pd
from ncbi_taxa import ncbiTaxa


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
            phylaDf = phylaDf.append(pd.DataFrame({'Domain': pd.Categorical([group]),
                                                   'Phylum': pd.Categorical([phylum]),
                                                   'TaxId': pd.Series([record['taxId']], dtype='int'),
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

    
    


def standalone():
    speciesByPhylaTable()
    return 0
    

if __name__=="__main__":
    import sys
    sys.exit(standalone())

