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
    
    

    
