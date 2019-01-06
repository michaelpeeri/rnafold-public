from data_helpers import setSpeciesProperty, allSpeciesSource, getSpeciesProperty



def runDistributed():

    for taxId in allSpeciesSource():

        currentProp = getSpeciesProperty(taxId, 'paired-mRNA-fraction')

        if currentProp[0] is None:
            continue

        if currentProp[1] == "computed":
            origVal = float(currentProp[0])
            fixedVal = origVal*2
            setSpeciesProperty( taxId, "paired-mRNA-fraction", "%.4g"%fixedVal, "computed (v2)", overwrite=True )
            print("Fixed %d: %.4g -> %.4g" % (taxId, origVal, fixedVal))
        
runDistributed()
