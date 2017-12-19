import pandas as pd
import math
from pyvirtualdisplay.display import Display  # use Xvnc to provide a headless X session (required by ete for plotting)
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, RectFace, faces, AttrFace, PhyloTree
from ncbi_taxa import ncbiTaxa
from data_helpers import allSpeciesSource
from ncbi_entrez import getKingdomForSpecies
from ncbi_taxa import ncbiTaxa

# Configuration
minimalTaxonSize = 9   # Note - this should match the value in tree_traits_analysis_with_taxgroups.r
csvGroupsFilename = "TaxidToKingdom.csv"
#nwFilename = "TaxidToKingdom.nw"
#csvRegressionEffectsByGroupCsv = "tree_traits_effects_analysis_with_taxgroups.out.csv"
csvRegressionEffectsByGroupCsv = "tree_traits_effects_analysis_with_taxgroups.out.save2.csv"



print("Processing all species...")
taxidToKingdom = {}
taxidToLineage = {}
for taxId in allSpeciesSource():    # read kingdoms and taxonomic lineages for all species
    kingdom = getKingdomForSpecies(taxId)
    taxidToKingdom[taxId] = kingdom

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


def writeMajorGroupsTable():
    # Create an empty data-frame for csv output
    df = pd.DataFrame({
    #    'kingdom': pd.Categorical([], categories=['Bacteria', 'Eukaryota', 'Archaea'], ordered=False)  # Do pandas categorical vars ever work?
        'kingdom': pd.Series(dtype="string")
        }, index=pd.Index([], name='tax_id', dtype='int') )

    # Add kingdom data to the data-frame
    for k,v in taxidToKingdom.items():
        df.loc[k, 'kingdom'] = v
        assert(df.loc[k, 'kingdom'] == v)

    # Get list of large taxonomic groups (based on the lineages of all species)
    majorGroups = getMajorTaxonomicGroups(taxidToLineage)

    # Add a binary membership column for each major group
    for groupTaxId, _ in majorGroups:
        groupName = ncbiTaxa.get_taxid_translator( [groupTaxId] )[groupTaxId]
        groupName = "Member_%s_%d" % (groupName.replace(" ", "_").replace("/", "_").replace("-", "_"), groupTaxId)

        groupDf = pd.DataFrame(
            {groupName: pd.Series(dtype='bool')},
            index = pd.Index(df.index.values, name='tax_id', dtype='int')
        )

        for taxId, lineage in taxidToLineage.items():
            isMember = int(groupTaxId in lineage)
            groupDf.loc[taxId, groupName] = isMember

        df = pd.merge( df, groupDf, how='inner', left_index=True, right_index=True )  # Add the new column (is there an easier way to do this?)

    print(df.shape)

    print(df.loc[511145,])

    df.to_csv(csvGroupsFilename)
    print("Wrote %s" % csvGroupsFilename)


fixedHyphenations = {"Gammaproteobacteria":"Gammaproteo-\nbacteria"}
def hyphenate(s):
    if s in fixedHyphenations:
        return fixedHyphenations[s]
    else:
        spacepos = s.find(" ")
        if spacepos > 4:
            return s.replace(" ", "\n")
        else:
            if len(s) < 10:
                return s

            slashpos = s.find("/") 
            if slashpos > 4:
                return s.replace("/", "/\n")
            
            return "%s-\n%s" % (s[:9], s[9:])
    

def plotRegressionEffectsByGroup():
    majorGroups = getMajorTaxonomicGroups( taxidToLineage )
    majorGroupsTree = ncbiTaxa.get_topology( [x[0] for x in majorGroups] )

    df = pd.read_csv( csvRegressionEffectsByGroupCsv )
    print(df)

    baseFontSize = 25

    allExplanatoryVars = set(df['ExplanatoryVar'])

    significanceLevel = 1e-2

    with Display(backend='xvnc') as disp:  # Plotting requires an X session

        for var in allExplanatoryVars:
            ts = TreeStyle()
            ts.show_branch_length = False
            ts.show_leaf_name = False
            ts.scale = 60     # scale for branch length
            #ts.rotation = 90
            ts.rotation = 0
            ts.min_leaf_separation = 5
            

            tree = majorGroupsTree.copy()
            
            taxidsToKeep = set(df[df['ExplanatoryVar']==var]['TaxGroup'].values)

            def nodeLayoutFunc(node):
                taxid = int(node.name)

                if taxid in taxidsToKeep:
                    taxGroupName = ncbiTaxa.get_taxid_translator([taxid])[taxid]  # There has to be an easier way to look up names...
                    
                    row =  df[(df['ExplanatoryVar'] == var) & (df['TaxGroup'] == taxid)]
                    effectSize = float(row['EffectSize'].values[0])
                    pval = float(row['Pvalue'].values[0])

                    significanceMarker = ""
                    if( pval < significanceLevel ):
                        significanceMarker = " (*)" #unichr(0x2731)

                    #name = TextFace("%s\np-val=%.1g%s\nR^2=%.2g\n(N=%d)" % (hyphenate(taxGroupName), pval, significanceMarker, effectSize, row['NumSpecies']) )#, fsize=baseFontSize) #, fstyle="italic")
                    name = TextFace("%s" % hyphenate(taxGroupName), fsize=baseFontSize*1.5 ) #, pval, significanceMarker ) ) #, significanceMarker, effectSize, row['NumSpecies']), fsize=baseFontSize) #, fstyle="italic")
                    #name = TextFace("%s" % hyphenate(taxGroupName), fsize=baseFontSize) #, fstyle="italic")
                    #print("[%s]" % hyphenate(taxGroupName))
                    #name.rotation = -90
                    name.margin_left=20
                    name.margin_right=0
                    name.margin_top=60
                    name.margin_bottom=40
                    #faces.add_face_to_node(name, node, column=1)
                    faces.add_face_to_node(name, node, column=0)

                    details = TextFace("p-val=%.1g%s\nR^2=%.2g\n(N=%d)" % ( pval, significanceMarker, effectSize, row['NumSpecies']), fsize=baseFontSize )#, fsize=baseFontSize) #, fstyle="italic")
                    #details = TextFace("%s\np-val=%.3g%s" % (hyphenate(taxGroupDetails), pval, significanceMarker ) ) #, significanceMarker, effectSize, row['NumSpecies']), fsize=baseFontSize) #, fstyle="italic")
                    #details = TextFace("%s" % hyphenate(taxGroupDetails), fsize=baseFontSize) #, fstyle="italic")
                    #print("[%s]" % hyphenate(taxGroupDetails))
                    #details.rotation = -90
                    details.margin_left=20
                    details.margin_right=0
                    details.margin_top=10
                    details.margin_bottom=90
                    #faces.add_face_to_node(details, node, column=1)
                    faces.add_face_to_node(details, node, column=0)
                    

                        
                    #v = RectFace(width=effectSize*70, height=baseFontSize*2, fgcolor="SteelBlue", bgcolor="SteelBlue", label={"text":"%.2g %s" % (effectSize, significanceMarker), "fontsize":baseFontSize+1, "color":"black" } )
                    #v.rotation = -90
                    #faces.add_face_to_node(v, node, column=1)
                    
                    nstyle = NodeStyle()
                    nstyle["size"] = math.sqrt( effectSize / math.pi ) * 400   # choose radius such that the area is proportional to effect size (since perception is proportional to area, not radius)
                    nstyle["shape"] = "sphere"
                    if( pval < significanceLevel ):
                        nstyle["fgcolor"] = "#1133ff"
                    else:
                        nstyle["fgcolor"] = "#90b0cc"
                        
                    node.set_style(nstyle)
                    # TODO - add text
                    
            
            ts.layout_fn = nodeLayoutFunc
            
            # Annotate tree nodes
            nodesToKeep = []
            for node in tree.traverse():
                if int(node.name) in taxidsToKeep:
                    nodesToKeep.append(node)

            tree.prune(nodesToKeep)


            tree.render('regressionByTaxgroup_%s.pdf' % var, tree_style=ts, units="mm")
            tree.render('regressionByTaxgroup_%s.svg' % var, tree_style=ts, units="mm")
        
    


def test():
    plotRegressionEffectsByGroup()


test()
