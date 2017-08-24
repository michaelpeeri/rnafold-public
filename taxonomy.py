from data_helpers import allSpeciesSource, getSpeciesTemperatureInfo, getSpeciesProperty
import re
from pyvirtualdisplay.display import Display  # use Xvnc to provide a headless X session (required by ete for plotting)
#from easyprocess import EasyProcess
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces, AttrFace, PhyloTree
from plot_xy import loadProfileData
from mfe_plots import heatmaplotProfiles, getProfileHeatmapTile
from reference_trees import pruneReferenceTree_Nmicrobiol201648
from ncbi_taxa import ncbiTaxa


# configuration
maxDisplayNameLength = 27
speciesToExclude = frozenset((99999, 9876543))
speciesToInclude = frozenset()

# Test pyvirtualdisplay
#with Display(backend='xvnc') as disp:
#    with EasyProcess('xmessage hello -timeout 5') as proc:
#        proc.wait()  # returns after 5 seconds



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
        name = AttrFace("name", fsize=9)
        #faces.add_face_to_node(n, node, 0, position="float")
        faces.add_face_to_node(name, node, column=0)

        tf = TextFace(len(node.get_leaves()), fsize=12)
        faces.add_face_to_node(tf, node, column=1, position="float")
        
    elif( level==0 ):
        tf = TextFace(len(node.get_leaves()), fsize=12)
        tf.margin_right = 5
        faces.add_face_to_node(tf, node, column=1, position="float")

    elif node.is_leaf():
        name = AttrFace("name", fsize=9, fstyle="italic")
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
        name = AttrFace("name", fsize=12, fstyle="italic")
        faces.add_face_to_node(name, node, column=0)
    


def drawTrees(completeTree, prunedTree):

    ts = TreeStyle()
    #ts.mode = "c"
    #ts.arc_start = -180
    #ts.arc_span  =  360
    ts.show_branch_length = False
    ts.show_leaf_name = False
    ts.layout_fn = simpleNodeLayoutFunc
    ts.scale = 100
    
    with Display(backend='xvnc') as disp:  # Plotting requires an X session
        drawTree( completeTree, 'nmicro_s6',        tree_style=ts,  w=300, h=3000, units="mm")
        drawTree( prunedTree,   'nmicro_s6_pruned', tree_style=ts,  w=300, h=1500, units="mm")

        # Display is about to close; how to tell tree to disconnect cleanly? (to prevent "Client Killed" message...)
        
    return 0



reMatchH5file_find_taxid = re.compile("(\w+_)+taxid_(\d+)(_\w+)+.h5")
class ProfileDataTileGenerator(object):
    def __init__(self, files, biasProfiles, dfProfileCorrs):
        self._biasProfiles = biasProfiles
        self._dfProfileCorrs = dfProfileCorrs

        self._yrange = heatmaplotProfiles(biasProfiles, None, dfProfileCorrs, None)
        
        #self._h5files = {}
        #self.parse_files(files)
        #print(self._h5files)

    #def parse_files(self, files):
    #    for fn in files:
    #        match = reMatchH5file_find_taxid.match(fn)
    #        if match:
    #            taxid = int(match.group(2))
    #            self._h5files[taxid] = fn

    def getProfileTile(self, taxid):
        return getProfileHeatmapTile(taxid, self._biasProfiles, self._dfProfileCorrs, self._yrange)
    
    def getProfileTileFunc(self):
        return lambda x: self.getProfileTile(x)
                
            
        

def testTileFunc(taxId):
    return "heatmap_profile_taxid_511145.png"

"""
Plot "statistical" tree, with species names and counts 
This tree should illustrate the included species
"""
def plotSpeciesOnTaxonomicTree(tileFunc=None):
    taxa = getSpeciesToInclude()

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

        tree.render('alltaxa.pdf', tree_style=ts, w=10, h=48, units="mm")
        tree.render('alltaxa.svg', tree_style=ts, w=10, h=48, units="mm")

        print("------------------ Ignore error message ------------------")
        # Display is about to close; how to tell tree to disconnect cleanly? (to prevent "Client Killed" message...)

    return 0

def standalone():
    import argparse
    from glob import glob
    import os
    
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument("--use-profile-data", type=str, default="")
    argsParser.add_argument("--test", action="store_true", default=False)
    args = argsParser.parse_args()


    if( args.test ):
        taxa = getSpeciesToInclude()
        (completeTree, prunedTree) = pruneReferenceTree_Nmicrobiol201648(taxa)
        drawTrees( completeTree, prunedTree )

        prunedTree.write(format=1, outfile="nmicro_s6_pruned.nw")

        return 0
    else:
        files = []
        if args.use_profile_data:
            files = [x for x in glob(args.use_profile_data) if os.path.exists(x)]

        print("Loading profile data for %d files..." % len(files))
        (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics) = loadProfileData(files)
            
        tileGenerator = ProfileDataTileGenerator(files, biasProfiles, dfProfileCorrs)
        
        #return plotSpeciesOnTaxonomicTree()
        return plotSpeciesOnTaxonomicTree(tileFunc=tileGenerator.getProfileTileFunc())


if __name__=="__main__":
    import sys
    sys.exit(standalone())


        
    
