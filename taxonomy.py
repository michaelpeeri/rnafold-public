from data_helpers import allSpeciesSource, getSpeciesTemperatureInfo, getSpeciesProperty
import re
from pyvirtualdisplay.display import Display  # use Xvnc to provide a headless X session (required by ete for plotting)
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces, AttrFace, PhyloTree
import numpy as np
from plot_xy import loadProfileData
from mfe_plots import heatmaplotProfiles, getProfileHeatmapTile
from reference_trees import pruneReferenceTree_Nmicrobiol201648
from ncbi_taxa import ncbiTaxa
from cluster_profiles import analyzeProfileClusters

# ---------------------------------------------------------
# configuration
# ---------------------------------------------------------
# Length at which to truncate taxon names
maxDisplayNameLength = 27

# Exclusion and inclusion list
speciesToExclude = frozenset((405948,999415,946362,470, 158189, 1307761, 456481, 505682, 1280, 4932, 508771, 2850, 753081 ))
speciesToInclude = frozenset()  # If non empty, only these species will be included

# Master font scale (points)
fontScale = 14


#
# Collapsed-tree plotting
#
# Distance threshold for deciding how many clusters to use for a given taxon
max_rms_distance_to_merge_clusters = 0.04
# Taxons below this level will be merged
collapedTaxonomicTreeLevel = 4
# Font-size for the taxon name at each level
levelFontSizes = {0: 1.5, 1: 2.1, 2:1.8, 3:1.5}
# Display distance between centroids (to help choose max_rms_distance_to_merge_clusters)
showInterCentroidDistances = False
# Use fixed Y scale (TODO - use adaptive scale, as in the other profile plots)
fixedYscale = 2.9
# For K-means, choose the number of initial points to use
#n_init=10000
n_init=5000

# ---------------------------------------------------------
# Choose plot format:
#
# single-page (long)
kingdomPageAssignment_SinglePage = {2: 0,     # Bacteria
                                    2759: 0,  # Eukaryotes
                                    2157: 0}  # Archaea
# 2-page
kingdomPageAssignment_2page      = {2: 1,     # Bacteria
                                    2759: 0,  # Eukaryotes
                                    2157: 0}  # Archaea

kingdomPageAssignment = kingdomPageAssignment_2page

# ---------------------------------------------------------



def getSpeciesToInclude():
    count = 0;
    taxa = []
    for taxId in allSpeciesSource():

        # Use black list
        if taxId in speciesToExclude:
            continue # Exclude this taxon

        # Use white list (if defined)
        if speciesToInclude:
            if taxId not in speciesToInclude:
                continue # Exclude this taxon

        # Taxon will be included

        taxa.append(taxId)
        count += 1

    print("Using %d species" % count)
    return taxa
    

# Color names: SVG colors. See: http://etetoolkit.org/docs/latest/reference/reference_treeview.html?highlight=imgface#ete3.SVG_COLORS
temperatureRangeToColor = {'Mesophilic': 'LightGreen', 'Hyperthermophilic': 'Crimson', 'Thermophilic': 'Orange', 'Psychrophilic': 'Cyan', 'Unknown': 'White'}

salinityToColor = {'NonHalophilic': 'LemonChiffon', 'Mesophilic': 'LemonChiffon', 'ModerateHalophilic': 'Turquoise', 'ExtremeHalophilic': 'Teal'}

oxygenReqToColor = {'Aerobic': 'CornflowerBlue', 'Anaerobic': 'LightSalmon', 'Facultative': 'Violet', 'Microaerophilic': 'Aquamarine'}

#Counter({'HostAssociated': 43, 'Aquatic': 32, 'Specialized': 25, 'Multiple': 25, 'Terrestrial': 9})
habitatToColor = {'HostAssociated': 'LightPink' , 'Aquatic': 'LightSeaGreen', 'Specialized': 'Gold', 'Multiple': 'Silver', 'Terrestrial': 'Sienna'}

algaeToColor = {'Yes': 'MediumTurquoise', 'No': 'LemonChiffon'}

def nodeLayoutWithTaxonomicNames(node, tileFunc=None):
    level = len(node.get_ancestors())
    if( level==1 or level==2 or level==3):
        if node.name:
            name = AttrFace("name", fsize=fontScale)
            ########faces.add_face_to_node(n, node, 0, position="float")
            faces.add_face_to_node(name, node, column=0)

        tf = TextFace(len(node.get_leaves()), fsize=fontScale*1.5)
        faces.add_face_to_node(tf, node, column=1, position="float")
        
    elif( level==0 ):
        tf = TextFace(len(node.get_leaves()), fsize=fontScale*1.5)
        tf.margin_right = 5
        faces.add_face_to_node(tf, node, column=1, position="float")

    elif node.is_leaf():
        assert(node.name)
        name = AttrFace("name", fsize=fontScale, fstyle="italic")
        faces.add_face_to_node(name, node, column=0)

        #-------------------------------
        # LFE profile
        #
        tile = None
        if not tileFunc is None:
            tile = tileFunc(node.taxId)
            
        if not tile is None:
            profileFace = faces.ImgFace(tile, width=150, height=10, is_url=False) # no support for margin_right?
            faces.add_face_to_node(profileFace, node, column=1, aligned=True)

        #proteinCount = getSpeciesProperty(taxId, 'protein-count')
        #if not proteinCount[0] is None:
        #    proteinCount = int(proteinCount)
                
        #genomeSizeMB = getSpeciesProperty(taxId, 'genome-size-mb')
        #if not genomeSizeMB[0] is None:
        #    genomeSizeMB = float(genomeSizeMB)
                    
        #genomicGC = getSpeciesProperty(taxId, 'gc-content')
        #if not genomicGC[0] is None:
        #    genomeicGC = float(genomicGC)

        #-------------------------------
        # genomicGC
        #
        genomicGC = getSpeciesProperty(node.taxId, 'gc-content')[0]
        if not genomicGC is None:
            genomicGC = float(genomicGC)
            genomicGCFace = faces.RectFace( width= (genomicGC-18.0)/(73.0-18.0)*50 , height=5, fgcolor="SteelBlue", bgcolor="SteelBlue", label={"text":"%.2g"%genomicGC, "fontsize":8, "color":"Black"} )
            genomicGCFace.margin_right = 5
            faces.add_face_to_node(genomicGCFace, node, column=2, aligned=True)
        
            
        #-------------------------------
        # Temperature
        #
        (numericalProp, categoricalProp) = getSpeciesTemperatureInfo(node.taxId)
        if categoricalProp[0] != "Unknown":
            temperatureColor = temperatureRangeToColor[categoricalProp[0]]
            temperatureFace = faces.RectFace( width=25, height=10, fgcolor=temperatureColor, bgcolor=temperatureColor )
            temperatureFace.margin_right = 5
            faces.add_face_to_node(temperatureFace, node, column=3, aligned=True)

        #-------------------------------
        # Salinity
        #
        salinity = getSpeciesProperty(node.taxId, 'salinity')[0]
        if salinity is None:
            salinity = "Unknown"
            
        if salinity != "Unknown":
            salinityColor = salinityToColor[salinity]
            salinityFace = faces.RectFace( width=25, height=10, fgcolor=salinityColor, bgcolor=salinityColor )
            salinityFace.margin_right = 5
            faces.add_face_to_node(salinityFace, node, column=4, aligned=True)


        #-------------------------------
        # Oxygen-req
        #
        oxygenReq = getSpeciesProperty(node.taxId, 'oxygen-req')[0]
        if oxygenReq is None:
            oxygenReq = "Unknown"
            
        if oxygenReq != "Unknown":
            oxygenReqColor = oxygenReqToColor[oxygenReq]
            oxygenReqFace = faces.RectFace( width=25, height=10, fgcolor=oxygenReqColor, bgcolor=oxygenReqColor )
            oxygenReqFace.margin_right = 5
            faces.add_face_to_node(oxygenReqFace, node, column=5, aligned=True)

        #-------------------------------
        # Habitat
        #
        habitat = getSpeciesProperty(node.taxId, 'habitat')[0]
        if habitat is None:
            habitat = "Unknown"
            
        if habitat != "Unknown":
            habitatColor = habitatToColor[habitat]
            habitatFace = faces.RectFace( width=25, height=10, fgcolor=habitatColor, bgcolor=habitatColor )
            habitatFace.margin_right = 5
            faces.add_face_to_node(habitatFace, node, column=6, aligned=True)
            
        #-------------------------------
        # Algae
        #
        algae = getSpeciesProperty(node.taxId, 'algae')[0]
        if algae is None:
            algae = "Unknown"
            
        if algae != "Unknown":
            algaeColor = algaeToColor[algae]
            algaeFace = faces.RectFace( width=25, height=10, fgcolor=algaeColor, bgcolor=algaeColor )
            algaeFace.margin_right = 5
            faces.add_face_to_node(algaeFace, node, column=7, aligned=True)

        #-------------------------------
        # paired fraction
        #
        pairedFraction = getSpeciesProperty(node.taxId, 'paired-mRNA-fraction')[0]
        if not pairedFraction is None:
            pairedFraction = float(pairedFraction)
            pairedFractionFace = faces.RectFace( width= pairedFraction*50 , height=5, fgcolor="SteelBlue", bgcolor="SteelBlue", label={"text":"%.2g"%(pairedFraction*100), "fontsize":8, "color":"Black"} )
            pairedFractionFace.margin_right = 5
            faces.add_face_to_node(pairedFractionFace, node, column=8, aligned=True)
            


"""
Return unary layout function that also receives a tile function (for getting image to plot along each node)
TODO - are there other parameters to layout funcs that should be passed?
"""
def makeNodeLayoutFuncWithTiles( layoutfn, tileFunc=None):
    return lambda node: layoutfn(node, tileFunc=tileFunc)

        




def drawTree(tree, basename, *args, **kw):

    tree.render('%s.pdf' % basename, *args, **kw)
    tree.render('%s.svg' % basename, *args, **kw)
        

def simpleNodeLayoutFunc(node):
    if node.is_leaf():
        if node.name != "n/a":
            l = int(node.name)
            binomicName = ncbiTaxa.get_taxid_translator([l])[l]  # There has to be an easier way to look up names...
            #name = AttrFace("name", fsize=12, fstyle="italic")
            name = TextFace(binomicName, fsize=fontScale*1.5, fstyle="italic")
            faces.add_face_to_node(name, node, column=0)

    else:

        try:
            k = node.testK
            if not k is None:
                #binomicName = ncbiTaxa.get_taxid_translator([k])[k]  # There has to be an easier way to look up names...

                name = TextFace(k, fsize=fontScale*1.5)
                #print("--- %s" % l)
                faces.add_face_to_node(name, node, column=0)

        except AttributeError as e:
            pass
        
        try:
            l = node.testL
            if not l is None:
                binomicName = ncbiTaxa.get_taxid_translator([l])[l]  # There has to be an easier way to look up names...

                name = TextFace(binomicName, fsize=fontScale*1.5)
                #print("--- %s" % l)
                faces.add_face_to_node(name, node, column=1)

        except AttributeError as e:
            pass
        
    


def drawTrees(completeTree, prunedTree, args=None):
    import os
    from glob import glob

    files = []
    if not args is None and args.use_profile_data:
        files = [x for x in glob(args.use_profile_data) if os.path.exists(x)]
    
    print("Loading profile data for %d files..." % len(files))
    (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics) = loadProfileData(files)
            
    tileGenerator = ProfileDataTileGenerator(files, biasProfiles)

    ts = TreeStyle()
    #ts.mode = "c"
    #ts.arc_start = -180
    #ts.arc_span  =  360
    ts.show_branch_length = False
    ts.show_leaf_name = False
    ts.layout_fn = simpleNodeLayoutFunc
    ts.scale = 1000
    
    with Display(backend='xvnc') as disp:  # Plotting requires an X session
        drawTree( completeTree, 'nmicro_s6',        tree_style=ts,  w=300, h=3000, units="mm")
        drawTree( prunedTree,   'nmicro_s6_pruned', tree_style=ts,  w=1500, h=1500, units="mm")

        plotSpeciesOnTaxonomicTree(tileFunc=tileGenerator.getProfileTileFunc(), tree=prunedTree)
        

        # Display is about to close; how to tell tree to disconnect cleanly? (to prevent "Client Killed" message...)
        
    return 0


"""
Plot each profile separataly, and return the filename of a 'tile' that can be included in other graphs.
"""
class ProfileDataTileGenerator(object):
    def __init__(self, files, biasProfiles):
        self._biasProfiles = biasProfiles

        self._yrange = heatmaplotProfiles(biasProfiles)  # plot all profiles, to determine a single color-scale that will fit them all

    def getProfileTile(self, taxid):
        return getProfileHeatmapTile(taxid, self._biasProfiles, self._yrange)

    """
    Return a free function that will return a tile based on the taxid
    """
    def getProfileTileFunc(self):
        return lambda x: self.getProfileTile(x)
                
            
        

def testTileFunc(taxId):
    return "heatmap_profile_taxid_511145.png"




"""
Plot "statistical" tree, with species names and counts 
This tree should illustrate the included species
"""
def plotSpeciesOnTaxonomicTree(tileFunc=None, tree=None):
    taxa = getSpeciesToInclude()

    if tree is None:
        # Get the smallest "taxonomic" (i.e., n-ary) tree that includes all specified species
        # See: http://etetoolkit.org/docs/3.0/tutorial/tutorial_ncbitaxonomy.html
        tree = ncbiTaxa.get_topology(taxa, intermediate_nodes=True)
    
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = makeNodeLayoutFuncWithTiles( nodeLayoutWithTaxonomicNames, tileFunc )
    ts.aligned_header.add_face(TextFace("LFE"), column=1)
    ts.aligned_header.add_face(TextFace("GC%"), column=2)
    ts.aligned_header.add_face(TextFace("Tmp"), column=3)
    ts.aligned_header.add_face(TextFace("Hal"), column=4)
    ts.aligned_header.add_face(TextFace("Oxy"), column=5)

    ts.aligned_header.add_face(TextFace("Hab"), column=6)
    ts.aligned_header.add_face(TextFace("Alg"), column=7)


    ts.aligned_header.add_face(TextFace("Paired"), column=8)

    for cat in ('Psychrophilic', 'Mesophilic', 'Thermophilic', 'Hyperthermophilic'):
        color = temperatureRangeToColor[cat]
        ts.legend.add_face(faces.RectFace( width=100, height=20, fgcolor=color, bgcolor=color, label={"text":cat, "color":"Black", "fontsize":8}), column=1)

    for cat in ('NonHalophilic', 'ModerateHalophilic', 'ExtremeHalophilic'):
        color = salinityToColor[cat]
        ts.legend.add_face(faces.RectFace( width=100, height=20, fgcolor=color, bgcolor=color, label={"text":cat, "color":"Black", "fontsize":8}), column=2)

    for cat in ('Aerobic', 'Facultative', 'Anaerobic', 'Microaerophilic'):
        color = oxygenReqToColor[cat]
        ts.legend.add_face(faces.RectFace( width=100, height=20, fgcolor=color, bgcolor=color, label={"text":cat, "color":"Black", "fontsize":8}), column=3)

    for cat in ('Aquatic', 'Terrestrial', 'HostAssociated', 'Specialized', 'Multiple'):
        color = habitatToColor[cat]
        ts.legend.add_face(faces.RectFace( width=100, height=20, fgcolor=color, bgcolor=color, label={"text":cat, "color":"Black", "fontsize":8}), column=4)

    for cat in ('Yes', 'No'):
        color = algaeToColor[cat]
        ts.legend.add_face(faces.RectFace( width=100, height=20, fgcolor=color, bgcolor=color, label={"text":cat, "color":"Black", "fontsize":8}), column=5)
        

        
    # Annotate tree nodes
    for node in tree.traverse():
        if not node.name:
            continue
        
        taxId = int(node.name)
        binomicName = ncbiTaxa.get_taxid_translator([taxId])[taxId]  # There has to be an easier way to look up names...
        node.add_features(displayName=binomicName, taxId=taxId)
        if( len(binomicName) > maxDisplayNameLength):
            node.name = binomicName[:maxDisplayNameLength-1].rstrip() + '...'
        else:
            node.name = binomicName

    
    n = 0
    nsLight = NodeStyle()
    #nsLight["bgcolor"] = "LightSt"
    nsDark  = NodeStyle()
    nsDark["bgcolor"] = "WhiteSmoke"
    for node in tree.traverse():
        if not node.is_leaf():
            continue
        
        if n%2==0:
            node.set_style(nsDark)
        else:
            node.set_style(nsLight)
        n += 1
        

    with Display(backend='xvnc') as disp:  # Plotting requires an X session

        h = int(round(len(taxa)*0.25, 0))
        tree.render('alltaxa.pdf', tree_style=ts, w=10, h=h, units="mm")
        tree.render('alltaxa.svg', tree_style=ts, w=10, h=h, units="mm")

        print("------------------ Ignore error message ------------------")
        # Display is about to close; how to tell tree to disconnect cleanly? (to prevent "Client Killed" message...)

    return 0


def makeProfilesArray(keys, biasProfiles):
    if( not keys ):
        raise Exception("No keys specified")

    # Filter keys with missing profiles
    keysFound = [key for key in keys if key in biasProfiles]
    if( not keysFound ):
        raise Exception("No profiles found matching requested keys: %s..." % str(keys))

    # Find out the length of the profile
    profileLength = biasProfiles[keysFound[0]].shape[0]  # We will use the length of the first profile, and later verify the following profiles match
    assert(profileLength > 3)

    # Copy the profiles into an array
    out = np.zeros( (len(keysFound), profileLength ) )
    for i, key in enumerate(keysFound):
        out[i,:] = biasProfiles[key]

    return out

def nodeLayoutForCollapsedTree(node, groupMembers, biasProfiles, tileFunc=None):
    level = len(node.get_ancestors())
    taxId = int(node.name)
    nodeName = ncbiTaxa.get_taxid_translator([taxId])[taxId]  # There has to be an easier way to look up names...
    members = groupMembers[taxId]
        
    if level > 0:
        name = TextFace(nodeName, fsize=fontScale*levelFontSizes[level])
        faces.add_face_to_node(name, node, column=0)

        
    if( level >= collapedTaxonomicTreeLevel - 2 ):
        assert(node.name)

        tile = None

        profilesArray = makeProfilesArray( members, biasProfiles )

        #print("//// %d %s" % (level, profilesArray.shape))

        if( profilesArray.shape[0] == 1 ):
            #if not taxId in biasProfiles:
            #    print("TaxId %d not found" % taxId)
                
            tileDummyId = 1900000000+taxId*10

            tile = getProfileHeatmapTile(tileDummyId, {tileDummyId:profilesArray[0]}, (-fixedYscale,fixedYscale) )

            if not tile is None:
                profileFace = faces.ImgFace(tile, width=150, height=fontScale, is_url=False)
                profileFace.margin_top = 4   # separate this taxon's profiles from the next
                faces.add_face_to_node(profileFace, node, column=1 )#, aligned=True)

                countFace = faces.TextFace( "1", fsize=fontScale )
                countFace.background.color = "#dfdfdf"
                countFace.margin_right = 3
                faces.add_face_to_node(countFace, node, column=2 )#, aligned=True)

            else:
                print("No node <- Level: %d profilesArray: %s" % (level, profilesArray.shape))
                
        else:
            # Perform clustering
            centers, labels, dist = analyzeProfileClusters( profilesArray, n_init, max_rms_distance_to_merge_clusters )
            groupCounts = [sum([1 for x in labels if x==i]) for i in range(len(centers))]  # count how many profiles belong to each group
            
            for i in range(centers.shape[0]):
                tileDummyId = 1900000000+taxId*10+i

                tile = getProfileHeatmapTile(tileDummyId, {tileDummyId:centers[i]}, (-fixedYscale,fixedYscale) )

                if not tile is None:
                    profileFace = faces.ImgFace(tile, width=150, height=fontScale, is_url=False)
                    profileFace.margin_bottom = 2   # separate this taxon's profiles from the next
                    if i==0:
                        profileFace.margin_top = 4   # separate this taxon's profiles from the next
                    faces.add_face_to_node(profileFace, node, column=1 )#, aligned=True)

                    countFace = faces.TextFace( groupCounts[i], fsize=fontScale )
                    countFace.background.color = "#dfdfdf"
                    countFace.margin_right = 3
                    faces.add_face_to_node(countFace, node, column=2 )#, aligned=True)

                    if showInterCentroidDistances:
                        distFace = faces.TextFace( "%.3g" % dist, fsize=fontScale )
                        #distFace.background.color = "#dfdfdf"
                        distFace.margin_left = 7
                        distFace.margin_right = 3
                        faces.add_face_to_node(distFace, node, column=3 )#, aligned=True)

                else:
                    print("No node <- Level: %d  Id: %d profilesArray: %s" % (level, i, profilesArray.shape))

    else:
        countFace = faces.TextFace( len(members), fsize=fontScale )
        countFace.background.color = "#dfdfdf"
        countFace.margin_right = 3
        faces.add_face_to_node(countFace, node, column=1 )#, aligned=True)
        
        


def plotCollapsedTaxonomicTree(biasProfiles):
    allTaxa = getSpeciesToInclude()

    # ---------------------------------------------------------
    collapsedTaxa = set()
    groupMembers = {}
    
    for taxid in allTaxa:
        lineage = ncbiTaxa.get_lineage(taxid)

        group = lineage[collapedTaxonomicTreeLevel]

        # Create a mapping containing all leaves found under a grouping taxon
        for ancestor in lineage[:collapedTaxonomicTreeLevel+1]:
            if ancestor in groupMembers:
                groupMembers[ancestor].append( taxid )
            else:
                groupMembers[ancestor] = [taxid]

        collapsedTaxa.add(group)


    def filterTaxaByPage(taxa, pageToInclude):
        out = set()
        for taxid in taxa:
            lineage = ncbiTaxa.get_lineage(taxid)
            taxonPage = kingdomPageAssignment[lineage[2]]
            if taxonPage == pageToInclude:
                out.add(taxid)
        return out
            
    
    with Display(backend='xvnc') as disp:  # Plotting requires an X session
        
        def plotCollapsedTree(taxa, i):
            print("Page %d (%d taxons)" % (i, len(taxa)))
            
            # Get the smallest "taxonomic" (i.e., n-ary) tree that includes all specified species
            # See: http://etetoolkit.org/docs/3.0/tutorial/tutorial_ncbitaxonomy.html
            # In this case, we provide mid-level taxons to get the "collapsed" tree
            tree = ncbiTaxa.get_topology(taxa, intermediate_nodes=True)


            # Plot the tree
            ts = TreeStyle()
            ts.show_leaf_name = False
            ts.layout_fn = lambda node: nodeLayoutForCollapsedTree(node, groupMembers, biasProfiles)


            h = int(round(len(taxa)*0.25, 0))
            tree.render('alltaxa.collapsed.%d.pdf' % i, tree_style=ts, w=7, h=h, units="mm")
            tree.render('alltaxa.collapsed.%d.svg' % i, tree_style=ts, w=7, h=h, units="mm")

        for pageNum in sorted(set(kingdomPageAssignment.values())):
            plotCollapsedTree(filterTaxaByPage(collapsedTaxa, pageNum), pageNum)
    
    print("------------------ Ignore error message ------------------")
    # Display is about to close; how to tell tree to disconnect cleanly? (to prevent "Client Killed" message...)

    return 0


def savePrunedTree(tree):
    tree.write(format=1, outfile="nmicro_s6_pruned_with_names.nw")
    
    tree = tree.copy()


    usedTaxIds = {}


    # Rewrite the tax-id into the name (so it appears in the output tree)
    # Also warn about collisions (duplicate nodes that were mapped to the tax-id)
    for node in tree.traverse():
        if not node.is_leaf():
            continue

        if node.taxId in usedTaxIds:
            print("=="*20)
            print("WARNING: Duplicate found for tax-id %d" % node.taxId)
            print("Old: %s" % usedTaxIds[node.taxId])
            print("New: %s" % node.label)
            print(node.lineageItems)
        else:
            usedTaxIds[node.taxId] = node.label # store to detect future collisions

        # Rewrite the tax-id into the name (so it appears in the output tree)
        node.name = node.taxId 
    
    tree.write(format=1, outfile="nmicro_s6_pruned_with_taxids.nw")


def standalone():
    import argparse
    from glob import glob
    import os
    
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument("--use-profile-data", type=str, default="")
    argsParser.add_argument("--use-tree", choices=("taxonomic", "hug", "taxonomic-collapsed"), default="taxonomic")
    args = argsParser.parse_args()


    if( args.use_tree=="hug" ):
        taxa = getSpeciesToInclude()
        (completeTree, prunedTree) = pruneReferenceTree_Nmicrobiol201648(taxa)
        drawTrees( completeTree, prunedTree, args=args )

        savePrunedTree( prunedTree )

        return 0
    
    elif( args.use_tree=="taxonomic-collapsed" ):

        if( not args.use_profile_data ):
            raise Exception()
        
        files = [x for x in glob(args.use_profile_data) if os.path.exists(x)]
        
        print("Loading profile data for %d files..." % len(files))
        (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics) = loadProfileData(files)
        
        return plotCollapsedTaxonomicTree(biasProfiles)
    
    elif( args.use_tree=="taxonomic" ):
        files = []

        if args.use_profile_data:
            files = [x for x in glob(args.use_profile_data) if os.path.exists(x)]

        print("Loading profile data for %d files..." % len(files))
        (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics) = loadProfileData(files)

        tileGenerator = ProfileDataTileGenerator(files, biasProfiles)
        
        #return plotSpeciesOnTaxonomicTree()
        return plotSpeciesOnTaxonomicTree(tileFunc=tileGenerator.getProfileTileFunc())


if __name__=="__main__":
    import sys
    sys.exit(standalone())


        
    
