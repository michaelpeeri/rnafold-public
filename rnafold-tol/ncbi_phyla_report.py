# Parse NCBI Entrez taxonomy report to determine all phyla-rank taxonomic groups
# Note: This could probably also be done by processing NCBI taxa data, but this way seemed easier...
#
# The report to be parsed was created using NCBI Entrez taxonomy query "Phylum[Rank]" (saved as XML format)
#
import xml.etree.ElementTree as ET

# configuartion
reportPath = "./data/taxonomy_result.xml"

# Define group names to group the data under; each group will collect items identified by a prefix of the <Lineage> element
# Note: speciesByPhylaTable() (in stats_reports.py) assumes the keys match NCBI taxonomic group names.
groupsOfInterest = {'Bacteria': 'cellular organisms; Bacteria', 'Eukaryota': 'cellular organisms; Eukaryota', 'Archaea': 'cellular organisms; Archaea'}

"""
Parse the taxonomy report; return a json-like dictionary mapping group names (e.g., 'Bacteria') to dictionaries further mapping phyla names to details
"""
def parseReport(reportFilename=reportPath):
    tree = None
    try:
        tree = ET.parse(reportFilename)
    except ET.ParseError as e:
        print(e)
        raise e
        
    count = 0

    phyla = {}
    for x in groupsOfInterest.keys():
        phyla[x] = {}
    

    for phylum in tree.findall("./Taxon"):
        taxId = int(phylum.find("./TaxId").text)
        assert(taxId)

        parentTaxId = int(phylum.find("./ParentTaxId").text)
        assert(parentTaxId)

        name = phylum.find("./ScientificName").text
        assert(name)

        assert(phylum.find("./Rank").text =="phylum")

        # Determine which domain this phylum belongs to. Note that many phyla are found inside unranked 'groups' (i.e., not at the top of the tree)
        # each domain will collect items identified by a prefix of the <Lineage> element

        lineageStr = phylum.find("./Lineage").text
        assert(lineageStr)

        group = None
        for g, prefix in groupsOfInterest.items():
            if lineageStr.startswith(prefix):
                group = g
                break

        if group is None:
            print("Skipping: %s  (%s)" % (name, lineageStr))
            continue # This phylum does not appear to belong to one of the groups of interest

        # Store information about this phylum
        phyla[group][name] = {'taxId':taxId, 'parentTaxId': parentTaxId}
        
        count += 1

    print("Found %d phyla" % count)
    return phyla


def standalone():
    report = parseReport(reportPath)
    for k,v in report.items():
        print(k)
        print(len(v))

    print("All bacterial phyla: ")
    print( sorted(report['Bacteria'].keys() ) )

    return 0

if __name__=="__main__":
    import sys
    sys.exit(standalone())

        
