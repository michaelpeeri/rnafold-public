import pandas as pd
import math
from pyvirtualdisplay.display import Display  # use Xvnc to provide a headless X session (required by ete for plotting)
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, RectFace, faces, AttrFace, PhyloTree
from ncbi_taxa import ncbiTaxa
from data_helpers import allSpeciesSource
from ncbi_taxa import ncbiTaxa

#--------------------------------------------------------------
# Configuration
minimalTaxonSize = 9   # Note - this should match the value in tree_traits_analysis_with_taxgroups.r

#--------------------------------------------------------------
# Input file - choose one

# Possible input 1: GLS regression result for individual taxonomic groups and position ranges (created by tree_traits_effects_analysis_with_taxgroups.r)
#csvRegressionEffectsByGroupCsv = "tree_traits_effects_analysis_with_taxgroups.out.dLFE.length.300.csv"
csvRegressionEffectsByGroupCsv = "tree_traits_effects_analysis_with_taxgroups.out.abs(dLFE).length.300.csv"

# Possible input 2: Raw profile value outliers, i.e., profiles with positive dLFE at position 0 (created by find_trait_values_outliers.r)
#csvRegressionEffectsByGroupCsv = "find_trait_values_outliers.out.dLFE.csv"
#--------------------------------------------------------------

baseFontSize = 25         # Scale factor for (most) text
significanceLevel = 1e-2  # p-values smaller than this will be marked as significant
barScale = 500            # Width of 100%-bar
#barScale = 20              # Width of 100%-bar
useOwnXserver = False
#--------------------------------------------------------------


print("Processing all species...")
taxidToLineage = {}
for taxId in allSpeciesSource():    # read taxonomic lineages for all species
    lineage = ncbiTaxa.get_lineage( taxId )
    taxidToLineage[taxId] = lineage


"""
Return a list of "major taxonomic groups" (i.e., groups having at least a minimum number of species)
Input - map taxid -> lineage
Output - collection of pairs (taxId, num_species)
"""
def getMajorTaxonomicGroups(taxidToLineage):
    bigGroups = {}
    for taxId, lineage in taxidToLineage.items():
        for item in lineage:
            if item in bigGroups:
                bigGroups[item] += 1
            else:
                bigGroups[item] = 1

    totalSpecies = len(taxidToLineage)

    out = [x for x in bigGroups.items() if x[1] < totalSpecies and x[1] >= minimalTaxonSize ]  # keep all groups with at least 20 species
    return sorted( out, key=lambda x: -(x[1]) )  # Sort by decreasing size


# fixedHyphenations = {"Gammaproteobacteria":"Gammaproteo-\nbacteria"}
# def hyphenate(s):
#     return(s)  # DEBUG ONLY
#     if s in fixedHyphenations:
#         return fixedHyphenations[s]
#     else:
#         spacepos = s.find(" ")
#         if spacepos > 4:
#             #return s.replace(" ", "\n")
#             return s[:spacepos]  # DEBUG ONLY
#         else:
#             if len(s) < 10:
#                 return s

#             return s[:9]  # DEBUG ONLY

#             slashpos = s.find("/") 
#             if slashpos > 4:
#                 return s.replace("/", "/\n")
            
#             return "%s-\n%s" % (s[:9], s[9:])

class DummyResourceManager(object):
    def __init__(self, *dummy1, **dummy2):  pass
    def __enter__(self): return self
    def __exit__(self, exc_type, exc_val, exc_tb):  pass


def plotRegressionEffectsByGroup():
    majorGroups = getMajorTaxonomicGroups( taxidToLineage )
    majorGroupsTree = ncbiTaxa.get_topology( [x[0] for x in majorGroups] )

    df = pd.read_csv( csvRegressionEffectsByGroupCsv )
    print(df)


    allExplanatoryVars = set(df['ExplanatoryVar'])

    allRanges = list(sorted(set(df['Range'])))

    if useOwnXserver:
        _disp = Display
    else:
        _disp = DummyResourceManager

    with _disp(backend='xvnc') as disp:  # Plotting requires an X session

        # Helper function for plotting a single tree, for the specified trait and with data from all specified ranges
        def plotSingleTree( var, ranges ):
            print("Plotting tree for trait %s (ranges: %s)" % (var, ranges))

            assert(type(ranges)==type(()))

            ts = TreeStyle()
            ts.show_branch_length = False
            ts.show_leaf_name = False
            ts.scale = 130     # scale for branch length
            #ts.rotation = 90
            ts.rotation = 0
            ts.min_leaf_separation = 250


            tree = majorGroupsTree.copy()

            taxidsToKeep = set(df[df['ExplanatoryVar']==var]['TaxGroup'].values)

            def nodeLayoutFunc(node):
                taxid = int(node.name)

                if taxid in taxidsToKeep:
                    taxGroupName = ncbiTaxa.get_taxid_translator([taxid])[taxid]  # There has to be an easier way to look up names...

                    row = None
                    rangeRows = None

                    print(len(ranges))
                    
                    if( len(ranges) == 1 ):
                        row       =  df[(df['ExplanatoryVar'] == var) & (df['TaxGroup'] == taxid) & (df['Range'] == ranges[0]) ]
                        assert(len(row)==len(ranges))
                    elif len(ranges) > 1:
                        row       =  df[(df['ExplanatoryVar'] == var) & (df['TaxGroup'] == taxid) & (df['Range'] == 0) ]
                        assert(len(row)==1)
                        rangeRows =  df[(df['ExplanatoryVar'] == var) & (df['TaxGroup'] == taxid) & (df['Range'].isin(set(ranges))) ]
                    else:
                        assert(False)
                        
                    overallPval = float(row['Pvalue'].values[0])

                    name = TextFace("%s" % taxGroupName, fsize=baseFontSize*2.5 )
                    name.tight_text=True
                    name.margin_left=20
                    name.margin_right=0
                    name.margin_top=40
                    name.margin_bottom=12
                    faces.add_face_to_node(name, node, column=0)

                    #print(rangeRows)

                    # For each range to be included in this plot, add a bar
                    for rangeId in ranges:
                        #print("rangeId = %s" % (rangeId))

                        rowForThisRange = None

                        if len(ranges)==1:
                            rowForThisRange =  row
                        else:
                            rowForThisRange =  rangeRows[rangeRows['Range'] == rangeId]
                            
                        assert(len(rowForThisRange)==1)

                        # Extract p-value and "effect-size" (signed R^2)
                        effectSize = float(rowForThisRange['EffectSize'].values[0])
                        pval       = float(rowForThisRange['Pvalue'].values[0])

                        # Set bar-graph color and significance markers
                        barColor = ""
                        significanceMarker = ""
                        if( pval < significanceLevel ):
                            significanceMarker = " %s" % unichr(0x2731)
                            
                            if effectSize < 0:
                                barColor = "#1133ff"
                            else:
                                barColor = "#ff3311"
                        else:  # not significant
                            if effectSize < 0:
                                barColor = "#b0b0f0"
                            else:
                                barColor = "#ccb090"

                        # Add the minus sign if needed
                        signChar = ""
                        if effectSize < 0:
                            signChar = unichr(0x2212)  # minus sign (more legible than a hypen...)
                                
                        v = RectFace(width=abs(effectSize)*barScale, height=baseFontSize*3.5, fgcolor=barColor, bgcolor=barColor, label={"text":"%s%.2g %s" % (signChar, abs(effectSize), significanceMarker), "fontsize":baseFontSize*1.8, "color":"black" } )
                        #v.rotation = -90
                        v.margin_top=1
                        v.margin_left=30
                        v.margin_right=8
                        v.margin_bottom=12
                        faces.add_face_to_node(v, node, column=0)


                    details = TextFace("N=%d" % row['NumSpecies'], fsize=baseFontSize*1.5 )#, fsize=baseFontSize) #, fstyle="italic")
                    details.background.color = "#dfdfdf"
                    details.margin_left=6
                    details.margin_right=20
                    #details.margin_top=5
                    #details.margin_bottom=0
                    faces.add_face_to_node(details, node, column=1)

                    nstyle = NodeStyle()
                    nstyle["size"] = 0

                    node.set_style(nstyle)


            ts.layout_fn = nodeLayoutFunc

            # Annotate tree nodes
            nodesToKeep = []
            for node in tree.traverse():
                if int(node.name) in taxidsToKeep:
                    nodesToKeep.append(node)

            tree.prune(nodesToKeep)

            rangesStr = ""
            if len(ranges)==1:
                rangesStr = str(ranges[0])
            else:
                rangesStr = "_".join(map(str, ranges))

            tree.render('regressionByTaxgroup_%s_range_%s.pdf' % (var, rangesStr), tree_style=ts, units="mm")
            tree.render('regressionByTaxgroup_%s_range_%s.svg' % (var, rangesStr), tree_style=ts, units="mm")
            # End of plotSingleTree()
            #-----------------------------------------------------------
            

        for var in allExplanatoryVars:
            # Results are organized by taxon and by range. Non-zero ranges are defined in tree_traits_effects_analysis_with_taxgroups.r.
            # The range '0' covers the entire section analyzed (i.e., the pyramid width in the same file).
            # For each trait, we will output one tree for per range, plus an extra tree showing all non-zero ranges.
            # e.g., if there are three ranges, we will output the following plots: [0], [1], [2], [3], [1,2,3].
            # The [0] range is useful for traits for which there is no dependence on the range, and the [1,2,3] is good for traits in
            # which there is such dependendece.
            for ranges in list(map( lambda x: (x,), allRanges)) + [tuple(sorted(set(allRanges)-set((0,)))),]:
                #print("%s %s" % (var, ranges))
                plotSingleTree( var, ranges )
            


def test():
    plotRegressionEffectsByGroup()


test()
