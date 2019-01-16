import csv
from data_helpers import setSpeciesProperty

# Configuration
viera_silva_supp_tableA1_mapping_filename = "./data/Viera_Silva_Supp_tableA1.csv.mapping.csv"
dryRun = False
overwriteValues = False


#{'Accept2 (repeat)': '0', '#Vieira_Silva_Species': 'Yersinia pestis', 'Internal_species': 'Aspergillus niger', 'TieWarning': 'True', 'CheckConcensus': '0', 'GrowthTime': '1.25', 'Accept (manual decision \xe2\x80\x93 store into redis yes/no)': '0', 'Internal_taxid': '5061'}
rowAcceptYesNo = 'Accept (manual decision \xe2\x80\x93 store into redis yes/no)'
rowTaxId = 'Internal_taxid'
rowGrowthTime = 'GrowthTime'

numAcceptedRows = 0
with open( viera_silva_supp_tableA1_mapping_filename, "r") as csvfile:
    for row in csv.DictReader(csvfile, delimiter="\t"):

        assert(len(row)==8)

        if int(row[rowAcceptYesNo]) != 1: continue

        growthTimeHours = float(row[rowGrowthTime]) # make sure value is valid floating-point number (however original string repr. will be stored)

        taxId = int(row[rowTaxId])

        if not dryRun:
            setSpeciesProperty( taxId, "growth-time-hours", row[rowGrowthTime], "Vieira-Silva table A1", overwrite=overwriteValues )
        else:
            print( taxId, "growth-time-hours", row[rowGrowthTime], "Vieira-Silva table A1", overwriteValues )
        
        numAcceptedRows += 1


print("numAcceptedRows: {}".format(numAcceptedRows))
            
                
