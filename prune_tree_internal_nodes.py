from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces, AttrFace, PhyloTree
import pandas as pd


nmicrobiol201648_s6_PATHd8        = "./data/nmicrobiol201648-s6.txt.nw.PATHd8.out.d8.nw"    # =Tree
nmicrobiol201648_s6_PATHd8_pruned = "./data/nmicrobiol201648-s6.txt.nw.PATHd8.out.d8.nw.pruned.nw"    # =Tree
nodeIdentifiersMappingTable_with_testing_csv = "reference_trees.hug.identifiers.with_testing.csv"
alignmentIdentifiersToTreeIdentifiersMappingTable_csv = "reference_trees.hug.alignment_ids.csv"


treefmt = None

def readTranslationMap():
    treeNodeIdentifiersDf = pd.read_csv(nodeIdentifiersMappingTable_with_testing_csv, dtype={'NodeLabel': 'string', 'DBIdentifier': 'string', 'DBIdentifierType': 'category', 'TaxId': 'int32'} )

    out = {}
    for row in treeNodeIdentifiersDf.itertuples():
        label = row.NodeLabel
        taxId = row.TaxId
        if not taxId is None and taxId != 0:
            out[label] = taxId
            
    return out


def readTranslationMap2():
    treeNodeIdentifiersDf = pd.read_csv(alignmentIdentifiersToTreeIdentifiersMappingTable_csv)
    print(treeNodeIdentifiersDf)

    out = {}

    for row in treeNodeIdentifiersDf.itertuples():
        out[row[2]] = row[1]
            
    return out



translationMap = readTranslationMap()
print(len(translationMap))
map2 = readTranslationMap2()


with open(nmicrobiol201648_s6_PATHd8, "r") as f:
    #with open(itol_newick, "r") as f:
    treefmt = f.read()
    
    # Strip internal support, and parse the resulting string
    tree = PhyloTree(treefmt)

    nodesToPreserve = []
    nodesSet = set()

    count0 = 0
    count1 = 0

    count21 = 0
    count22 = 0

    for node in tree.traverse():
        if not node.is_leaf():
            continue
        #print(node.name)
        if map2[node.name] in translationMap:
            count1 += 1
        else:
            count0 += 1
            print(node.name)

        if node.name in map2:
            count21 += 1
            

    print(count0, count1)
    print(count21, count22)
    #tree.prune(nodesToPreserve, preserve_branch_length=True)
    
    #tree.write(format=1, outfile=nmicrobiol201648_s6_PATHd8_pruned)
    
    

    
