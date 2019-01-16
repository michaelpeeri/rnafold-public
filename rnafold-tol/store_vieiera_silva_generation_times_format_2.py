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
from data_helpers import setSpeciesProperty

# Configuration
viera_silva_supp_tableA1_mapping_filename = "./data/journal.pgen.1000808.s005.DOC.csv.with.mapping.csv"
dryRun = False
overwriteValues = False


#{'Accept2 (repeat)': '0', '#Vieira_Silva_Species': 'Yersinia pestis', 'Internal_species': 'Aspergillus niger', 'TieWarning': 'True', 'CheckConcensus': '0', 'GrowthTime': '1.25', 'Accept (manual decision \xe2\x80\x93 store into redis yes/no)': '0', 'Internal_taxid': '5061'}
rowAcceptYesNo = 'Decision (manual decision yes/no)'
rowAcceptYesNo2 = 'Decision2'
rowDisagreement = 'Disagreement'

rowTaxId = 'Interal_tax_id'
rowGrowthTime = 'GrowthTimeHours'

numAcceptedRows = 0
with open( viera_silva_supp_tableA1_mapping_filename, "r") as csvfile:
    for row in csv.DictReader(csvfile):

        assert(len(row)==9)
        assert(row[rowAcceptYesNo]==row[rowAcceptYesNo2])
        assert(row[rowDisagreement]=="0")

        if int(row[rowAcceptYesNo]) != 1: continue

        growthTimeHours = float(row[rowGrowthTime]) # make sure value is valid floating-point number (however original string repr. will be stored)

        taxId = int(row[rowTaxId])

        if not dryRun:
            setSpeciesProperty( taxId, "growth-time-hours-v2", row[rowGrowthTime], "Vieira-Silva table A1 (format2)", overwrite=overwriteValues )
        else:
            print( taxId, "growth-time-hours-v2", row[rowGrowthTime], "Vieira-Silva table A1 (format2)", overwriteValues )
        
        numAcceptedRows += 1


print("numAcceptedRows: {}".format(numAcceptedRows))
            
                
