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
# Create a tax-id to taxonomic-group memberships table for all species (for use by R code...)
import pandas as pd
import csv  # for QUOTE_NONNUMERIC
from data_helpers import allSpeciesSource, getSpeciesShortestUniqueNamesMapping, getSpeciesName
from ncbi_entrez import getKingdomForSpecies
from ncbi_taxa import ncbiTaxa


# Configuration
minimalTaxonSize = 9   # Note - this should match the value in tree_traits_analysis_with_taxgroups.r
csvFilename = "TaxidToKingdom.csv"


print("Processing all species...")
taxidToKingdom = {}
taxidToLineage = {}
for taxId in allSpeciesSource():    # read kingdoms and taxonomic lineages for all species
    kingdom = getKingdomForSpecies(taxId)
    taxidToKingdom[taxId] = kingdom

    lineage = ncbiTaxa.get_lineage( taxId )
    taxidToLineage[taxId] = lineage


"""
Return a list of "major taxonomic groups" (i.e., groups having at least a minimum number of species)
Input - map taxid -> lineage
Output - collection of pairs (taxId, num_species)
"""
def getMajorTaxonomicGroups(taxidToLineage):
    bigGroups = {}
    for taxId, lineage in taxidToLineage.items():
        for item in lineage:
            if item in bigGroups:
                bigGroups[item] += 1
            else:
                bigGroups[item] = 1

    totalSpecies = len(taxidToLineage)

    out = [x for x in bigGroups.items() if x[1] < totalSpecies and x[1] >= minimalTaxonSize ]  # keep all groups with at least 20 species
    return sorted( out, key=lambda x: -(x[1]) )  # Sort by decreasing size


shortNames = getSpeciesShortestUniqueNamesMapping()

# Create an empty data-frame for csv output
df = pd.DataFrame({
    'kingdom': pd.Series(dtype="string"),
    'full.name': pd.Series(dtype="string"),
    'short.name': pd.Series(dtype="string")
    }, index=pd.Index([], name='tax_id', dtype='int') )

# Add kingdom data to the data-frame
for k,v in taxidToKingdom.items():
    df.loc[k, 'kingdom'   ] = v
    df.loc[k, 'full.name' ] = getSpeciesName(k)
    df.loc[k, 'short.name'] = shortNames[k]
    assert(df.loc[k, 'kingdom'] == v)

# Get list of large taxonomic groups (based on the lineages of all species)
majorGroups = getMajorTaxonomicGroups(taxidToLineage)

# Add a binary membership column for each major group
for groupTaxId, _ in majorGroups:
    groupName = ncbiTaxa.get_taxid_translator( [groupTaxId] )[groupTaxId]
    groupName = "Member_%s_%d" % (groupName.replace(" ", "_").replace("/", "_").replace("-", "_"), groupTaxId)

    groupDf = pd.DataFrame(
        {groupName: pd.Series(dtype='bool')},
        index = pd.Index(df.index.values, name='tax_id', dtype='int')
    )

    for taxId, lineage in taxidToLineage.items():
        isMember = int(groupTaxId in lineage)
        groupDf.loc[taxId, groupName] = isMember
    
    df = pd.merge( df, groupDf, how='inner', left_index=True, right_index=True )  # Add the new column (is there an easier way to do this?)

print(df.shape)

print(df.loc[511145,])
    
df.to_csv(csvFilename, quoting=csv.QUOTE_NONNUMERIC)
print("Wrote %s" % csvFilename)

#majorGroupsTree = ncbiTaxa.get_topology( [x[0] for x in majorGroups] )
#majorGroupsTree.write(format=1, outfile=nwFilename)
#print("Wrote %s" % nwFilename)
# Moved to plot_tree_effects_anaylsis_results_on_tree.py
