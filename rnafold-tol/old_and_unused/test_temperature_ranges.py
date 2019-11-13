from collections import Counter
from data_helpers import getSpeciesProperty, allSpeciesSource, getSpeciesTemperatureInfo
import pandas as pd
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style


categories = Counter()
salinityCategories = Counter()
oxygenReqCategories = Counter()

habitatCategories = Counter()


#Counter({'Mesophilic': 79, 'Hyperthermophilic': 20, 'Thermophilic': 13, 'Psychrophilic': 4, 'Unknown': 1})


# Temperatures and categories for all species *that have temperatures*
temperatureVsCategoryStatistics = pd.DataFrame({
    'tax_id':pd.Series(dtype='int'),
    'temperature':pd.Series(dtype='float'),
    'category':pd.Categorical([])
    })


# Plot raw data
for taxId in allSpeciesSource():

    category = None
    temperatureRange = getSpeciesProperty( taxId, 'temperature-range')
    if not temperatureRange[0] is None:
        category = temperatureRange[0]
        categories.update((category,))
    else:
        category = "Unknown"
    assert(not category is None)

    optimalTemperatureData = getSpeciesProperty( taxId, 'optimum-temperature')
    optimalTemperature = None
    if not optimalTemperatureData[0] is None:
        optimalTemperature = float(optimalTemperatureData[0])

        temperatureVsCategoryStatistics = temperatureVsCategoryStatistics.append(pd.DataFrame({
            'tax_id':pd.Series([taxId], dtype='int'),
            'temperature':pd.Series([optimalTemperature], dtype='float'),
            'category':pd.Categorical([category]) }))

        

print(categories)

print(temperatureVsCategoryStatistics)

fig, ax = plt.subplots()

temperatureVsCategoryStatistics.boxplot(column='temperature', by='category', ax=ax)

plt.savefig("temperatures.pdf")
plt.savefig("temperatures.svg")
plt.close(fig)



# Temperatures and categories for all species 
temperatureVsCategoryStatistics2 = pd.DataFrame({
    'tax_id':pd.Series(dtype='int'),
    'temperature':pd.Series(dtype='float'),
    'category':pd.Categorical([]),
    'category_source':pd.Categorical([]),
    'salinity':pd.Categorical([]),
    'oxygen_req':pd.Categorical([]),
    'genome_size':pd.Series(dtype='float'),
    'protein_count':pd.Series(dtype='float'),
    'habitat':pd.Categorical([])
    })

# Plot raw data
for taxId in allSpeciesSource():

    (numericalProp, categoricalProp) = getSpeciesTemperatureInfo(taxId)
    print("%s %s" % (numericalProp, categoricalProp))

    assert(not categoricalProp[0] is None)

    optimalTemperature = None
    if not numericalProp[0] is None:
        optimalTemperature = float(numericalProp[0])

    salinity = getSpeciesProperty(taxId, 'salinity')
    if not salinity[0] is None:
        salinityCategories.update([salinity[0]])
        salinty = salinity[0]

        if salinity=="Mesophilic":
            salinity=="NonHalophilic"
    else:
        salinity = "Unknown"

    oxygenReq = getSpeciesProperty(taxId, 'oxygen-req')
    if not oxygenReq[0] is None:
        oxygenReqCategories.update([oxygenReq[0]])
        oxygenReq = oxygenReq[0]

    proteinCount = getSpeciesProperty(taxId, 'protein-count')
    if not proteinCount[0] is None:
        proteinCount = int(proteinCount[0])

    genomeSizeMB = getSpeciesProperty(taxId, 'genome-size-mb')
    if not genomeSizeMB[0] is None:
        genomeSizeMB = float(genomeSizeMB[0])

    genomicGC = getSpeciesProperty(taxId, 'gc-content')
    if not genomicGC[0] is None:
        genomicGC = float(genomicGC[0])

    habitat = getSpeciesProperty(taxId, 'habitat')
    if not habitat[0] is None:
        habitatCategories.update([habitat[0]])
        habitat = habitat[0]
        
    
    #{'NonHalophilic': 16, 'Mesophilic': 6, 'ModerateHalophilic': 3, 'ExtremeHalophilic': 2}

    temperatureVsCategoryStatistics2 = temperatureVsCategoryStatistics2.append(pd.DataFrame({
        'tax_id':pd.Series([taxId], dtype='int'),
        'temperature':pd.Series([optimalTemperature], dtype='float'),
        'category':pd.Categorical([categoricalProp[0]]),
        'category_source':pd.Categorical([categoricalProp[1]]),
        'salinity':pd.Categorical([salinity]),
        'oxygen_req':pd.Categorical([oxygenReq]),
        'genome_size':pd.Series([genomeSizeMB], dtype='float'),
        'protein_count':pd.Series([proteinCount], dtype='float'),
        'habitat':pd.Categorical([habitat])
    }))
    
print(temperatureVsCategoryStatistics2)

fig, ax = plt.subplots()

temperatureVsCategoryStatistics2.boxplot(column='temperature', by='category', ax=ax)

plt.savefig("temperatures2.pdf")
plt.savefig("temperatures2.svg")
plt.close(fig)

print(salinityCategories)
print(oxygenReqCategories)
print(habitatCategories)

