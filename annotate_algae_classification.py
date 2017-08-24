from random import randint
from data_helpers import getAllNativeCDSsForSpecies, decompressNucleicSequence, countSpeciesCDS, setSpeciesProperty, allSpeciesSource, getSpeciesProperty
from ete3 import NCBITaxa


ncbiTaxa = NCBITaxa()


# configuration


# From google doc:
positiveGroup = ((3067,'None'),(3068,'None'),(3055,'None'),(3046,'None'),(436017,'None'),(70448,'None'),(564608,'None'),(574566,'None'),(2769,'None'),(280699,'None'),(130081,'None'),(35688,'None'),(2788,'None'),(753081,'None'),(641309,'None'),(280463,'None'),(195065,'None'),(905079,'None'),(2898,'None'),(464988,'None'),(556484,'None'),(296543,'None'),(159749,'None'),(635003,'None'),(1519565,'None'),(44056,'None'))

negativeGroup = ((559292,'None'),(284811,'None'),(436907,'None'),(237561,'None'),(336722,'None'),(1047168,'None'),(1397361,'None'),(4927,'None'),(1041607,'None'),(1069680,'None'),(420778,'None'),(1287680,'None'),(284812,'None'),(5061,'None'),(330879,'None'),(242507,'None'),(418459,'None'),(214684,'None'),(511145,'None'),(585056,'None'),(386585,'None'),(685038,'None'),(198214,'None'),(99287,'None'),(420890,'None'),(505682,'None'),(347257,'None'),(184922,'None'),(353152,'None'), (164328, 'None'), (4781, 'None'), (67593, 'None'), (1223560, 'None'), (695850, 'None'), (65357, 'None'), (352472, 'None'))


# Species belonging to the following groups are not algae, by definition (source: http://enwp.org/Algae)
algaeDefinition_ExcludedGroups = frozenset((4751, # Fungi
                                           3193, # Embryophyta
                                           33208 # Metazoa
))

algaeDefinition_IncludedGroups = frozenset((3041,  # Chlorophyta
                                            2763,  # Rhodophyta
                                            38254, # Glaucocystophyceae
                                            29197, # Chlorarachniophyceae
                                            3035,  # Euglenida
                                            33634, # Stramenopiles
                                            2830,  # Haptophyceae
                                            3027,  # Cryptophyta
                                            2864   # Dinophyceae
))



def run():

    positiveDict = dict(positiveGroup)
    negativeDict = dict(negativeGroup)

    totalCount = 0
    positiveCount = 0
    negativeCount = 0

    for taxId in allSpeciesSource():

        totalCount += 1

        #if not getSpeciesProperty(taxId, 'algae')[0] is None:
        #    continue

        lineage = frozenset(ncbiTaxa.get_lineage(taxId))

        algeaClassification = None

        if lineage.intersection( algaeDefinition_ExcludedGroups ):
            algeaClassification = ('No', 'Excluded taxonomic group')

        elif taxId in positiveDict:
            algeaClassification = ('Yes', positiveDict[taxId])
            
        elif taxId in negativeDict:
            algeaClassification = ('No', negativeDict[taxId])

        if not algeaClassification is None:
            setSpeciesProperty( taxId, "algae", algeaClassification[0], algeaClassification[1], overwrite=True )

            # Done; update counts
            if algeaClassification[0]=='Yes':
                positiveCount += 1
                
                if lineage.intersection( algaeDefinition_ExcludedGroups ):
                    print("Warning: possible false annotation: %d" % taxId)
                if not lineage.intersection( algaeDefinition_IncludedGroups ):
                    print("Warning: possible false annotation: %d" % taxId)
            elif  algeaClassification[0]=='No':
                negativeCount += 1
            else:
                assert(False)
        else:
            if lineage.intersection( algaeDefinition_IncludedGroups ):
                print("Warning: check unannotated possible algae: %d" % taxId)
    

    print("Finished %d species (%d annotated; %d positive, %d negative)" % (totalCount, positiveCount+negativeCount, positiveCount, negativeCount))
    
#test()
run()
