from collections import Counter
import pandas as pd
from tabulate import tabulate
from data_helpers import allSpeciesSource, getSpeciesName, getSpeciesTranslationTable, getSpeciesProperty
from endosymbionts import isEndosymbiont
from reference_trees import pruneReferenceTree_Nmicrobiol201648
from mfe_plots import getSpeciesShortestUniqueNamesMapping_memoized

#speciesToExclude = frozenset((470, 1280, 4932, 2850, 753081, 195065, 641309 ))
speciesToExclude = frozenset(( 195065, 641309 ))


speciesDf = pd.DataFrame({
    'TaxId': pd.Series([], dtype='int'),           # Species TaxId
    'Species': pd.Series([], dtype='str'),         # Species binomial name
    'Nickname': pd.Series([], dtype='str'),        #
    'Source' : pd.Series([], dtype='str'),         #
    'TranslationTbl': pd.Series([], dtype='int'),  #
    'InPhyloTree': pd.Series([], dtype='bool'),
    'GenomicGC%': pd.Series([], dtype='float'),
    'GenomicENc\'': pd.Series([], dtype='float'),
    'GrowthTempC': pd.Series([], dtype='float'),
    'GenomeSizeMb': pd.Series([], dtype='float'),
    'GrowthTimeHours': pd.Series([], dtype='float'),
    'IsEndosymbiont': pd.Series([], dtype='bool'),
    'EndosymbiontRef': pd.Series([], dtype='str')
})

def getSpeciesFromTree(tree):
    taxids = []
    for node in prunedTree.traverse():
        if not node.name:
            if( node.is_leaf() ):
                print("Warning: node.name missing for node {}".format(node.taxId))
            continue
        
        taxids.append( int(node.name) )
    return frozenset(taxids)

def getSpeciesToInclude():
    ret = []
    for taxId in allSpeciesSource():
        if taxId in speciesToExclude:
            continue
        ret.append(taxId)
    return ret

# get tree
(_, prunedTree) = pruneReferenceTree_Nmicrobiol201648(getSpeciesToInclude()) # prune complete reference phylogenetic tree to include only dataset species
speciesInTree = getSpeciesFromTree(prunedTree)


shortNames = getSpeciesShortestUniqueNamesMapping_memoized()


stats = Counter()
for taxId in getSpeciesToInclude():
    
    genomicGC = getSpeciesProperty(taxId, 'gc-content')[0]
    if not genomicGC is None:
        genomicGC = float(genomicGC)

    genomicENcprime = getSpeciesProperty(taxId, 'ENc-prime')[0]
    if not genomicENcprime is None:
        genomicENcprime = float(genomicENcprime)

    optimumTemp = getSpeciesProperty(taxId, 'optimum-temperature')[0]
    if not optimumTemp is None:
        optimumTemp = float(optimumTemp)

    genomeSizeMb = getSpeciesProperty(taxId, 'genome-size-mb')[0]
    if not genomeSizeMb is None:
        genomeSizeMb = float(genomeSizeMb)

    growthTimeHours = getSpeciesProperty(taxId, 'growth-time-hours-v2')[0]
    if not growthTimeHours is None:
        growthTimeHours = float(growthTimeHours)

    inPhyloTree = taxId in speciesInTree
    if inPhyloTree:
        stats.update(['tree'])

    speciesDf = speciesDf.append( pd.DataFrame({
        'TaxId': pd.Series([taxId], dtype='int'),
        'Species': pd.Series([getSpeciesName(taxId)], dtype='str'),
        'Nickname': pd.Series([shortNames[taxId]], dtype='str'),
        'Source': pd.Series([''], dtype='str'),
        'TranslationTbl': pd.Series([getSpeciesTranslationTable(taxId)], dtype='int'),
        'InPhyloTree': pd.Series([inPhyloTree], dtype='bool'),
        'GenomicGC%': pd.Series([genomicGC], dtype='float'),
        'GenomicENc\'': pd.Series([genomicENcprime], dtype='float'),
        'GrowthTempC': pd.Series([optimumTemp], dtype='float'),
        'GenomeSizeMb': pd.Series([genomeSizeMb], dtype='float'),
        'GrowthTimeHours': pd.Series([growthTimeHours], dtype='float'),
        'IsEndosymbiont': pd.Series([isEndosymbiont(taxId)], dtype='bool'),
        'EndosymbiontRef': pd.Series([''], dtype='str')
         }) )
    stats.update(["included"])

speciesDf = speciesDf.sort_values(by=['Species' ])    # sort rows
speciesDf.to_html( 'environment_report.html', float_format='{0:.1f}'.format )

with open("environment_report.rst", "w") as f:
    #f.write( speciesDf.drop(['RowType', 'Warnings', 'CDSWarnings', 'CDSWarnings_', 'FirstAA', 'LastAA', 'CDSDifference'], axis=1).pipe( tabulate, headers='keys', tablefmt='rst' ) )
    f.write( speciesDf.pipe( tabulate, headers='keys', tablefmt='rst' ) )

speciesDf.to_excel('environment_report.xlsx', sheet_name='Genomic and env. attributes')
    
print(speciesDf)
print(stats)

