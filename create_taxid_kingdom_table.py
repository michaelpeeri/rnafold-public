# Create a tax-id to kingdom table for all species (for use by R code...)
from data_helpers import allSpeciesSource
from ncbi_entrez import getKingdomForSpecies
import pandas as pd


# Configuration
csvFilename = "TaxidToKingdom.csv"


taxidToKingdom = {}
for taxId in allSpeciesSource():
    kingdom = getKingdomForSpecies(taxId)
    taxidToKingdom[taxId] = kingdom


df = pd.DataFrame({
#    'kingdom': pd.Categorical([], categories=['Bacteria', 'Eukaryota', 'Archaea'], ordered=False)  # Do pandas categorical vars ever work?
    'kingdom': pd.Series(dtype="string")
    }, index=pd.Index([], name='tax_id', dtype='int') )

for k,v in taxidToKingdom.items():
    df.loc[k, 'kingdom'] = v
    assert(df.loc[k, 'kingdom'] == v)

df.to_csv(csvFilename)
print("Wrote %s" % csvFilename)

