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
# Create species mapping table for vieira (internal species vs. table species) for manual confirmation.
#
# Full import process:
# > Supp. table A1 (format 2 - through DOI, orig. filename is journal.pgen.1000808.s005.DOC) is in word format - OCR not needed this time (just copy-paste and save as .csv)...
#
# > First, use edit-distanct (fuzzy_text_matching.py) to automatically map to internal species
# python2 create_auto_species_mapping_for_vieira_silva_data_table_format_2.py > ./data/journal.pgen.1000808.s005.DOC.csv.with.mapping.csv
#
# > Second, manually edit automatic species mapping to accept/reject inclusion of species. Data is matched at the species level (i.e., all strains belonging to the same genus and species are deemed equivalent). In some cases, NCBI synonyms are used to make an acceptable mapping.
#
# > Third, store the growth rates (for records which have manually accepted species matches)
# python store_vieiera_silva_generation_times_format_2.py
#
import csv
from data_helpers import allSpeciesSource, getSpeciesName
from fuzzy_text_matching import sw

slowFastGrowthCutoffHours = 2.5

viera_silva_supp_tableA1_filename = "./data/journal.pgen.1000808.s005.DOC.csv"


speciesMapping = {}
def loadSpeciesMapping():
    for taxId in allSpeciesSource():
        speciesName = getSpeciesName(taxId)
        speciesMapping[speciesName] = taxId

loadSpeciesMapping()

def findMatchingName(name):
    bestScore = 0.0
    bestMatch = ""
    bestMatchTaxId = 0
    isBestMatchATie = False
    for target, targetTaxId in speciesMapping.items():
        match = sw.align(name, target)
        #if match.score>0:
        #    print("{}\t{}".format( match.score, target) )

        if match.score==bestScore: # tie
            isBestMatchATie = True
            
        elif match.score > bestScore and match.score > 0:
            bestMatch = target
            bestMatchTaxId = targetTaxId
            bestScore = match.score
            isBestMatchATie = False
    
    return (bestMatch, bestMatchTaxId, isBestMatchATie)
        

#print(findMatchingName("Escherichia coli K12"))
#import sys
#sys.exit(0)

totalRows = 0
numRecordsWith1overMuValue = 0
bestScore = 0
bestMatch = ""
with open( viera_silva_supp_tableA1_filename, "r") as csvfile:
    for row in csv.DictReader(csvfile, delimiter="\t"):
        totalRows += 1
        try:
            # Perform some validity checks (since the csv file was converted from PDF using OCR)
            assert(len(row)==5)
            
            growthTime = float(row["d (h)"])

            numRecordsWith1overMuValue += 1

            speciesName = row['Species name']
            matchingName, matchingTaxId, isBestMatchATie = findMatchingName(speciesName)
            print("{}\t{}\t{}\t{}\t{}".format(speciesName, matchingName, matchingTaxId, isBestMatchATie, growthTime))
                
        except Exception as e:
            print(e)
            print(row)
            raise e

print("#Found {} records with 1/mu value (from a total of {} records).".format(numRecordsWith1overMuValue, totalRows))
print(bestMatch)
                
