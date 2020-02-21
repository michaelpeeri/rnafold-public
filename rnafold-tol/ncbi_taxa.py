from ete3 import NCBITaxa


ncbiTaxa = NCBITaxa()

def getDomainForSpecies(taxId):
    lineage = frozenset(ncbiTaxa.get_lineage(taxId))
    if 2 in lineage:
        return "Bacteria"
    elif 2157 in lineage:
        return "Archaea"
    elif 2759 in lineage:
        return "Eukaryotes"
    else:
        raise Exception("Couldn't determine domain of species {}".format(taxId))

