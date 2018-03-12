# Note: use 'export OPENBLAS_NUM_THREADS=1'
#       to limit automatica parallelization done by BLAS (during clustering)
from data_helpers import allSpeciesSource, getSpeciesTemperatureInfo, getSpeciesProperty
import re
from pyvirtualdisplay.display import Display  # use Xvnc to provide a headless X session (required by ete for plotting)
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, StaticItemFace, faces, AttrFace, PhyloTree
#from PyQt4 import QtCore, QtGui
import numpy as np
from mfe_plots import getHeatmaplotProfilesValuesRange, getProfileHeatmapTile, getLegendHeatmapTile, loadProfileData, loadPhylosignalProfiles, PCAForProfiles, getHeatmaplotProfilesValuesRange
from reference_trees import pruneReferenceTree_Nmicrobiol201648, pruneTreeByTaxonomy, extendTreeWithSpecies, getTaxidsFromTree
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
fontScale = 16

# Dimensions for drawing profiles
profileDrawingWidth  = 150
profileDrawingHeight = fontScale*1.5

#
profileStepNt  = 10
profileWidthNt = 40

#
# Collapsed-tree plotting
#
# Distance threshold for deciding how many clusters to use for a given taxon
#max_distance_to_split_clusters = 1.30  # previous level (with 4 levels)
max_distance_to_split_clusters = 1.60
# Maximum number of clusters (for each taxon)
max_clusters = 9
# Taxons below this level will be merged
collapedTaxonomicTreeLevel = 3
# Font-size for the taxon name at each level
levelFontSizes = {0: 1.5, 1: 2.1, 2:1.8, 3:1.5}
# Display distance between centroids (to help choose max_distance_to_split_clusters)
showInterCentroidDistances = False
# Use fixed Y scale (TODO - use adaptive scale, as in the other profile plots)
yScale = None
# For K-means, choose the number of initial points to use
#n_init=10000
n_init=10000


# Phylosignal profiles
#phylosignalProfilesFile = "phylosignalLipaOutputFile.csv"

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

def nodeLayoutWithTaxonomicNames(node, tileFunc=None, hideEnvironmentalVars=False):
    nstyle = NodeStyle()
    node.set_style(nstyle)
    
    nstyle["size"] = 0     # hide node symbols
    
    
    level = len(node.get_ancestors())
    if not node.is_leaf():
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

    else: #node.is_leaf():
        assert(node.name)
        name = AttrFace("name", fsize=fontScale, fstyle="italic")
        faces.add_face_to_node(name, node, column=0)

        #-------------------------------
        # LFE profile
        #
        tile = None
        if not tileFunc is None:
            print("Getting tile for {}...".format(node.taxId))
            tile = tileFunc(node.taxId)
            print("...{}".format(tile))
            
        if not tile is None:
            profileFace = faces.ImgFace(tile, width=profileDrawingWidth, height=profileDrawingHeight, is_url=False) # no support for margin_right?
            faces.add_face_to_node(profileFace, node, column=1, aligned=True)
        else:
            print("No tile for node {}".format(node.taxId))

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
        

        if( not hideEnvironmentalVars):
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

        # Hide the lines of dummyTopology nodes (i.e., nodes for which the tree topology is unknown and were added under the top node)
        if "dummyTopology" in node.features and node.dummyTopology:
            nstyle["hz_line_width"] = 0
            nstyle["hz_line_color"] = "#ffffff"
            nstyle["vt_line_width"] = 0
            nstyle["vt_line_color"] = "#ffffff"

"""
Return unary layout function that also receives a tile function (for getting image to plot along each node)
TODO - are there other parameters to layout funcs that should be passed?
"""
def makeNodeLayoutFuncWithTiles( layoutfn, tileFunc=None, hideEnvironmentalVars=False):
    return lambda node: layoutfn(node, tileFunc=tileFunc, hideEnvironmentalVars=hideEnvironmentalVars)

        




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
        

def collectProfilesFromTree(tree, profiles):
    ret = {}
    for taxId in getTaxidsFromTree(tree):
        if taxId in profiles:
            ret[taxId] = profiles[taxId]

    return ret


class DummyResourceManager(object):
    def __init__(self, *dummy1, **dummy2):  pass
    def __enter__(self): return self
    def __exit__(self, exc_type, exc_val, exc_tb):  pass
    

def drawTrees(completeTree, prunedTree, args=None):
    import os
    from glob import glob

    print("drawTrees()")

    files = []
    if not args is None and args.use_profile_data:
        files = [x for x in glob(args.use_profile_data) if os.path.exists(x)]
    
    print("Loading profile data for %d files..." % len(files))
    (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics) = loadProfileData(files)
            
    phylosignalProfiles = None
    if( args.use_phylosignal_data ):
        phylosignalProfiles = loadPhylosignalProfiles( args.use_phylosignal_data )
        
    print("(initializing tile-gen)")
    tileGenerator = ProfileDataTileGenerator(files, biasProfiles, dfProfileCorrs, phylosignalProfiles)

    ts = TreeStyle()
    #ts.mode = "c"
    #ts.arc_start = -180
    #ts.arc_span  =  360
    ts.show_branch_length = False
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.layout_fn = simpleNodeLayoutFunc
    ts.scale = 1500
    #ts.branch_vertical_margin = 10  # Y-axis spacing
    prunedTree.dist = 0.1  # Cut the root distance (since we are displaying a partial tree)

    # PCA
    profilesForPCA = collectProfilesFromTree(prunedTree, biasProfiles )

    PCAForProfiles( profilesForPCA, tileGenerator.getYRange(), profilesYOffsetWorkaround=args.profiles_Y_offset_workaround, profileScale=args.profile_scale, fontSize=args.font_size, overlapAction="hide", highlightSpecies=args.highlight_species )
    return

    if args.X_server:
        _disp = Display
    else:
        _disp = DummyResourceManager
    
    with _disp(backend='xvnc') as disp:  # Plotting requires an X session
        for treeToPlot, title in ( (completeTree, 'nmicro_s6'), (prunedTree, 'nmicro_s6_pruned') ):
            h = int(round(len(treeToPlot)*3.5, 0))
            
            if not args.hide_environmental_vars:
                h = int(round(len(treeToPlot)*6.0, 0))
            print("(drawTree {}) h={}".format(title, h))
            drawTree( treeToPlot, title,        tree_style=ts,  w=100, h=h, units="mm")
            #profilesForPCA = collectProfilesFromTree(treeToPlot, biasProfiles )
            #drawPCAforTree( profilesForPCA, [0]*len(profilesForPCA) )
        
        plotSpeciesOnTaxonomicTree(tileFunc=tileGenerator.getProfileTileFunc(), tree=prunedTree, hideLegend=args.hide_legend, hideEnvironmentalVars=args.hide_environmental_vars, treeScale=1000)
        

        # Display is about to close; how to tell tree to disconnect cleanly? (to prevent "Client Killed" message...)
        
    return 0


"""
Plot each profile separataly, and return the filename of a 'tile' that can be included in other graphs.
"""
class ProfileDataTileGenerator(object):
    def __init__(self, files, biasProfiles, corrData, phylosignalProfiles=None):
        self._biasProfiles = biasProfiles

        self._yrange = getHeatmaplotProfilesValuesRange(biasProfiles)  # plot all profiles once, to determine a single color-scale that will fit them all

        self._phylosignalProfiles = phylosignalProfiles

    def getProfileTile(self, taxid):
        return getProfileHeatmapTile(taxid, self._biasProfiles, self._yrange, phylosignalProfiles=self._phylosignalProfiles)

    """
    Return a free function that will return a tile based on the taxid
    """
    def getProfileTileFunc(self):
        return lambda x: self.getProfileTile(x)
                
    def getYRange(self):
        return self._yrange
        

def testTileFunc(taxId):
    return "heatmap_profile_taxid_511145.png"


"""
Plot "statistical" tree, with species names and counts 
This tree should illustrate the included species
"""
def plotSpeciesOnTaxonomicTree(tileFunc=None, tree=None, phylosignalFile="", hideLegend=False, hideEnvironmentalVars=False, ownXserver=False, treeScale=100):
    taxa = getSpeciesToInclude()

    if tree is None:
        # Get the smallest "taxonomic" (i.e., n-ary) tree that includes all specified species
        # See: http://etetoolkit.org/docs/3.0/tutorial/tutorial_ncbitaxonomy.html
        tree = ncbiTaxa.get_topology(taxa, intermediate_nodes=True)
    
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.scale = treeScale
    #ts.branch_vertical_margin = 10  # Y-axis spacing
    ts.layout_fn = makeNodeLayoutFuncWithTiles( nodeLayoutWithTaxonomicNames, tileFunc, hideEnvironmentalVars=hideEnvironmentalVars )

    # Add header (with a caption above each column)
    ts.aligned_header.add_face(TextFace("LFE"), column=1)
    ts.aligned_header.add_face(TextFace("GC%"), column=2)
    if( not hideEnvironmentalVars):
        ts.aligned_header.add_face(TextFace("Tmp"), column=3)
        ts.aligned_header.add_face(TextFace("Hal"), column=4)
        ts.aligned_header.add_face(TextFace("Oxy"), column=5)
        ts.aligned_header.add_face(TextFace("Hab"), column=6)
        ts.aligned_header.add_face(TextFace("Alg"), column=7)
        ts.aligned_header.add_face(TextFace("Paired"), column=8)

    
    if( not hideLegend ):
        # Draw legend
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
            if( node.is_leaf() ):
                print("Warning: node.name missing for node {}".format(node.taxId))
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

    if ownXserver:
        _disp = Display
    else:
        _disp = DummyResourceManager

    with _disp(backend='xvnc') as disp:  # Plotting requires an X session

        if not hideEnvironmentalVars:
            h = int(round(len(tree)*2.5, 0))
        else:
            h = int(round(len(tree)*4.0, 0))
        
        tree.render('alltaxa.pdf', tree_style=ts, w=100, h=h, units="mm")
        tree.render('alltaxa.svg', tree_style=ts, w=100, h=h, units="mm")

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

    nstyle = NodeStyle()
    nstyle["size"] = 0     # hide node symbols
    node.set_style(nstyle)
        
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

            tile = getProfileHeatmapTile(tileDummyId, {tileDummyId:profilesArray[0]}, yScale )

            if not tile is None:
                profileFace = faces.ImgFace(tile, width=profileDrawingWidth, height=profileDrawingHeight, is_url=False)
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
            centers, labels, dist = analyzeProfileClusters( profilesArray, n_init, max_distance_to_split_clusters, max_clusters )  # use default metric
            groupCounts = [sum([1 for x in labels if x==i]) for i in range(len(centers))]  # count how many profiles belong to each group
            
            for i in range(centers.shape[0]):
                tileDummyId = 1900000000+taxId*10+i

                tile = getProfileHeatmapTile(tileDummyId, {tileDummyId:centers[i]}, yScale )

                if not tile is None:
                    profileFace = faces.ImgFace(tile, width=profileDrawingWidth, height=profileDrawingHeight, is_url=False)
                    profileFace.margin_bottom  = 2   # separate this taxon's profiles from the next
                    if i==0:
                        profileFace.margin_top = 8   # separate this taxon's profiles from the next
                    faces.add_face_to_node(profileFace, node, column=1 )#, aligned=True)

                    countFace = faces.TextFace( groupCounts[i], fsize=fontScale*0.95 )
                    countFace.background.color = "#dfdfdf"
                    countFace.margin_right = 3
                    countFace.margin_top = 0
                    countFace.margin_bottom = 0
                    faces.add_face_to_node(countFace, node, column=2 )#, aligned=True)

                    if showInterCentroidDistances:
                        distFace = faces.TextFace( "%.4g" % dist, fsize=fontScale )
                        #distFace.background.color = "#dfdfdf"
                        distFace.margin_left = 7
                        distFace.margin_right = 3
                        faces.add_face_to_node(distFace, node, column=3 )#, aligned=True)

                else:
                    print("No node <- Level: %d  Id: %d profilesArray: %s" % (level, i, profilesArray.shape))

    else:
        pass
        #countFace = faces.TextFace( len(members), fsize=fontScale*1.5 )
        #countFace.background.color = "#dfdfdf"
        #countFace.margin_right = 3
        #faces.add_face_to_node(countFace, node, column=1 )#, aligned=True)
        
        


def plotCollapsedTaxonomicTree(biasProfiles, ownXserver=False):
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

    global yScale
    yScale = getHeatmaplotProfilesValuesRange( biasProfiles )

    if ownXserver:
        _disp = Display
    else:
        _disp = DummyResourceManager
    
    with _disp(backend='xvnc') as disp:  # Plotting requires an X session
        
        def plotCollapsedTree(taxa, i):
            print("Page %d (%d taxons)" % (i, len(taxa)))
            
            # Get the smallest "taxonomic" (i.e., n-ary) tree that includes all specified species
            # See: http://etetoolkit.org/docs/3.0/tutorial/tutorial_ncbitaxonomy.html
            # In this case, we provide mid-level taxons to get the "collapsed" tree
            tree = ncbiTaxa.get_topology(taxa, intermediate_nodes=True)


            # Configure tree plotting settings
            ts = TreeStyle()
            ts.show_leaf_name = False
            ts.branch_vertical_margin = 6  # Y-axis spacing
            ts.layout_fn = lambda node: nodeLayoutForCollapsedTree(node, groupMembers, biasProfiles)
            ts.show_scale = False
            ts.margin_top = 25
            ts.margin_bottom = 50


            # Create figure legend
            #if i==max(kingdomPageAssignment.values()):  # only on last page
            if i==0:                                     # only on first page

                ts.legend_position = 1 #=top left

                lastProfileStart = (biasProfiles[biasProfiles.keys()[0]].shape[0] - 1) * profileStepNt
                widthFor100nt = profileDrawingWidth * (100. / lastProfileStart)
                scaleFace = faces.RectFace( width=widthFor100nt, height=profileDrawingHeight, fgcolor="black", bgcolor="black" )
                scaleFace.margin_top     = 10
                scaleFace.margin_bottom  = 5
                scaleFace.margin_left    = 20
                scaleFace.margin_right   = 50
                ts.legend.add_face( scaleFace, column=0 )
                scaleFaceText = TextFace("100 bp", fsize=fontScale*1.5 )
                scaleFaceText.margin_top = 0
                ts.legend.add_face( scaleFaceText, column=0 )
                
                widthForWindow = profileDrawingWidth * (float(profileWidthNt) / lastProfileStart)
                windowWidthFace = faces.RectFace( width=widthForWindow, height=profileDrawingHeight, fgcolor="black", bgcolor="black" )
                windowWidthFace.margin_top    = 10
                windowWidthFace.margin_bottom = 5
                windowWidthFace.margin_left   = 20
                windowWidthFace.margin_right  = 50
                ts.legend.add_face( windowWidthFace, column=1 )
                windowWidthFace = TextFace("40 bp", fsize=fontScale*1.5 )
                windowWidthFace.margin_top    = 0
                ts.legend.add_face( windowWidthFace, column=1 )
                
                # Draw color scale
                
                #ts.legend.add_face( TextFace(unichr(0x1d6ab)+"LFE", fsize=12    ), column=0 )
                #ts.legend.add_face( TextFace("\xf0\x9d\x9a\xab"+"LFE", fsize=12    ), column=0 )
                dLFEFace = TextFace("dLFE", fsize=fontScale*2.0    )
                dLFEFace.margin_bottom = 0
                ts.legend.add_face( dLFEFace, column=2 )

                unitsFace = TextFace("kcal/mol/window", fsize=fontScale*1.5 )
                unitsFace.margin_top = 0
                ts.legend.add_face( unitsFace, column=2 )
                
                legendColorScaleTile = getLegendHeatmapTile(yScale)
                legendColorScaleFace = faces.ImgFace(legendColorScaleTile, width=480, height=80, is_url=False) # no support for margin_right?
                ts.legend.add_face( legendColorScaleFace, column=3 )

                # Test rich text plotting, using Qt (this didn't work...)
                #txtItem = QtGui.QGraphicsTextItem( "Hello, World! " + unichr(0x1d6ab)+ "LFE" )
                #txtItem = QtGui.QGraphicsSimpleTextItem()
                #txtItem.setPlainText( QtCore.QString( "Hello, World! " ) )   # <-- This step crashes the interpreter silently
                
                #siFace = StaticItemFace(txtItem)
                #siFace.width = 100
                #siFace.height = 20
                #ts.legend.add_face( siFace, column=5 )

            # Render and save the tree
            if collapedTaxonomicTreeLevel==4:
                h = int( round(  (len(tree)             *4.5), 0) )
            else:
                h = int( round( ((len(tree) + (i==0)*3 )*6.5), 0) )
                
            tree.render('alltaxa.collapsed.%d.pdf' % i, tree_style=ts, w=100, h=h, units="mm")
            tree.render('alltaxa.collapsed.%d.svg' % i, tree_style=ts, w=100, h=h, units="mm")

        # Plot each page (with a different subset of the tree)
        for pageNum in sorted(set(kingdomPageAssignment.values())):
            plotCollapsedTree(filterTaxaByPage(collapsedTaxa, pageNum), pageNum)
    
    print("------------------ Ignore error message ------------------")
    # Display is about to close; how can we disconnect cleanly? (to prevent "Client Killed" message...)

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

    
def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert

def standalone():
    import argparse
    from glob import glob
    import os

    argsParser = argparse.ArgumentParser()
    argsParser.add_argument("--use-profile-data", type=str, default="")
    argsParser.add_argument("--use-tree", choices=("taxonomic", "hug", "taxonomic-collapsed"), default="taxonomic")
    argsParser.add_argument("--use-phylosignal-data", type=str, default="")
    argsParser.add_argument("--limit-taxonomy", type=int, default=None)
    argsParser.add_argument("--include-all-species", action="store_true", default=False)
    argsParser.add_argument("--hide-legend", action="store_true", default=False)
    argsParser.add_argument("--hide-environmental-vars", action="store_true", default=False)
    argsParser.add_argument("--profiles-Y-offset-workaround", type=float, default=0.0)
    argsParser.add_argument("--profile-scale", type=float, default=1.0)
    #argsParser.add_argument("--profile-height", type=float, default=1.0)
    argsParser.add_argument("--font-size", type=float, default=7)
    argsParser.add_argument("--X-server", action="store_true", default=False)
    argsParser.add_argument("--highlight-species", type=parseList(int), default=())



    
    args = argsParser.parse_args()

    phylosignalProfiles = None
    if( args.use_phylosignal_data ):
        phylosignalProfiles = loadPhylosignalProfiles( args.use_phylosignal_data )

    if( args.use_tree=="hug" ):
        taxa = getSpeciesToInclude()
        
        (completeTree, prunedTree) = pruneReferenceTree_Nmicrobiol201648(taxa) # prune complete reference phylogenetic tree to include only dataset species 
        if( not args.limit_taxonomy is None ):  # limit taxonomy to specified taxon
            prunedTree = pruneTreeByTaxonomy( prunedTree, args.limit_taxonomy  )

        if args.include_all_species:
            print("{} =ext=> ".format(len(prunedTree)))
            prunedTree = extendTreeWithSpecies( tree=prunedTree, additionalSpecies=taxa, limitTaxonomy=args.limit_taxonomy )
            print(" =ext=> {}".format(len(prunedTree)))
        
        drawTrees( completeTree, prunedTree, args=args )

        savePrunedTree( prunedTree )

        return 0
    
    elif( args.use_tree=="taxonomic-collapsed" ):

        if( not args.use_profile_data ):
            raise Exception()
        
        files = [x for x in glob(args.use_profile_data) if os.path.exists(x)]

        print("Loading profile data for %d files..." % len(files))
        (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics) = loadProfileData(files)
        
        return plotCollapsedTaxonomicTree(biasProfiles, ownXserver=args.X_server)
    
    elif( args.use_tree=="taxonomic" ):
        files = []

        if args.use_profile_data:
            files = [x for x in glob(args.use_profile_data) if os.path.exists(x)]

        print("Loading profile data for %d files..." % len(files))
        (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics) = loadProfileData(files)

        tileGenerator = ProfileDataTileGenerator(files, biasProfiles, dfProfileCorrs, phylosignalProfiles)
        
        #return plotSpeciesOnTaxonomicTree()
        return plotSpeciesOnTaxonomicTree(tileFunc=tileGenerator.getProfileTileFunc(), hideLegend=args.hide_legend, hideEnvironmentalVars=args.hide_environmental_vars, ownXserver=args.X_server, treeScale=7 )


if __name__=="__main__":
    import sys
    sys.exit(standalone())


        
    
