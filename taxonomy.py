from data_helpers import allSpeciesSource, getSpeciesTemperatureInfo, getSpeciesProperty
import re
from pyvirtualdisplay.display import Display  # use Xvnc to provide a headless X session (required by ete for plotting)
#from easyprocess import EasyProcess
from ete3 import NCBITaxa, Tree, TreeStyle, NodeStyle, TextFace, faces, AttrFace, PhyloTree
from rate_limit import RateLimit
from plot_xy import loadProfileData
from mfe_plots import heatmaplotProfiles, getProfileHeatmapTile


# configuration
maxDisplayNameLength = 27
speciesToExclude = frozenset((99999, 9876543))
speciesToInclude = frozenset()
nmicrobiol201648_s6 = "./data/nmicrobiol201648-s6.txt"
nmicrobiol201648_s8 = "./data/nmicrobiol201648-s8.txt"
itol_newick         = "./data/itol_newick.txt"

# Test pyvirtualdisplay
#with Display(backend='xvnc') as disp:
#    with EasyProcess('xmessage hello -timeout 5') as proc:
#        proc.wait()  # returns after 5 seconds


ncbiTaxa = NCBITaxa()

rl = RateLimit(30)

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
        # paired fraction
        #
        pairedFraction = getSpeciesProperty(node.taxId, 'paired-mRNA-fraction')[0]
        if not pairedFraction is None:
            pairedFraction = float(pairedFraction)
            pairedFractionFace = faces.RectFace( width= pairedFraction*50 , height=5, fgcolor="SteelBlue", bgcolor="SteelBlue", label={"text":"%.2g"%(pairedFraction*100), "fontsize":8, "color":"Black"} )
            pairedFractionFace.margin_right = 5
            faces.add_face_to_node(pairedFractionFace, node, column=6, aligned=True)
            


"""
Return unary layout function that also receives a tile function (for getting image to plot along each node)
TODO - are there other parameters to layout funcs that should be passed?
"""
def makeNodeLayoutFuncWithTiles( layoutfn, tileFunc=None):
    return lambda node: layoutfn(node, tileFunc=tileFunc)

        
reComment = re.compile('[[][^]]*[]]')
def stripTreeInternalSupport(tree):
    return re.sub(reComment, '', tree)


def pruneTree(tree, keepNodes):
    taxIds = set()
    for x in keepNodes:
        if not x.taxId is None:
            taxIds.add(x.taxId)

    for node in tree.traverse(strategy='postorder'):
        childTaxIds = set()
        
        for x in node.get_leaves():
            try:
                if not x.taxId is None:
                    childTaxIds.add(x.taxId)
            except AttributeError as e:
                pass
        #print(len(childTaxIds))

        if not taxIds.intersection(childTaxIds):
            node.detach()
            #print(".")
            

def getContainingSpeciesLevelTaxon(taxId):
    lineage = ncbiTaxa.get_lineage(taxId)
    ranks = ncbiTaxa.get_rank(lineage)

    for taxon in reversed(lineage):
        rank = ranks[taxon]
        if rank == "no rank":
            continue
        elif rank == "subspecies":
            continue
        elif rank == "forma":
            continue
        elif rank == "varietas":
            continue
        elif rank == "species":
            return taxon
        else:
            #print("--"*50)
            #print(taxId)
            #print("Lineage: %s" % lineage)
            #print("Ranks: %s" % ranks)
            #print(ncbiTaxa.get_taxid_translator(lineage))
            return None

taxonItemsToIgnore = frozenset(("Unclassified", "Incertae", "Candidate", "candidate", "Candidatus"))

def testPlotTree():
    treefmt = None
    
    #with open(nmicrobiol201648_s6, "r") as f:
    with open(itol_newick, "r") as f:
        treefmt = f.read()
            
    tree = PhyloTree(stripTreeInternalSupport(treefmt))

    ts = TreeStyle()
    #ts.show_leaf_name = True
    #ts.layout_fn = nodeLayoutWithTaxonomicNames

    unmatched = []

    #xxx = set([1287680, 115713, 203267, 188937, 4781, 187420, 243230, 130081, 227882, 228908, 227377, 224308, 5693, 345663, 208964, 224325, 1116230, 243273, 213585, 64091, 45670, 1069680, 1397361, 280699, 1047168, 284811, 284812, 46234, 418459, 214684, 262768, 243365, 273063, 511145, 176299, 272557, 272558, 402612, 283166, 223926, 163003, 559292, 1041607, 1183438, 2769, 122586, 273116, 593117, 192222, 1574623, 243159, 160490, 212717, 272623, 272631, 272632, 63737, 272634, 1341181, 1125630, 99287, 27923, 400667, 269084, 257314, 96563, 300852, 4927, 381764, 242507, 65357, 104782, 336722, 190304, 882, 347515, 353152, 83332, 93061, 194439, 1223560, 267671, 196164, 1245935, 449447, 420778, 195522, 556484, 5061, 391623, 70601, 85962, 272844, 259536, 272633, 220668, 169963, 295405, 237561, 407035, 997884, 1432061, 1010810, 562, 1010800])
        
    numProcessed = 0
    numMatched = 0
    # Annotate tree leaves
    for node in tree.traverse():
        if node.is_leaf():
            items = node.name.split("_") # "items" (words) are separated by '_'

            #for i,x in enumerate(items):
            #    if x.startswith("Submit") or x.startswith("submit"):
            #        print("Removing submission note on node: %s" % items)
            #        items = items[:i]
            #        break
            
            matchingName = None
            matchingTaxId = None
            #print("---------------------------")
            #print(n)

            # Try to match the species name with the NCBI taxonomy
            # The node names contain concatenated taxon names and we don't know where the species name starts.
            # We will collect candidate strings (based on the criteria below) and choose the candidate that has a match in NCBI.
            candidates = []
            for i in range(len(items)-1, 0, -1): # process suffixes (strings containing last N items), from shortest to longest
                candidate = " ".join(items[i:])

                # Filter candidates that are too short or too long to be a binomial species name (with optional strain description)
                if len(candidate) < 8:
                    continue # skip this suffix
                
                if len(candidate) > 100:
                    break # skip this and longer suffixes

                if items[i] in taxonItemsToIgnore:
                    continue

                if candidate.find("ncertae sedis") > -1:
                    continue
                
                candidates.append( (i,candidate) )

                # For each suffix, also consider its prefixes as candidates. This is bacause the species names sometimes include additional strain identifiers that are not included in NCBI
                if i<=len(items)-2 and items[i] and items[i+1] and items[i][0].isupper() and items[i+1][0].islower(): # require that the first item start with a capital, and the second item start with a lower-case letter (e.g., "Eschericia coli")
                    for j in range(i+1,len(items)-1):
                        candidate = " ".join(items[i:j+1])

                        candidates.append( (i,candidate) )
                        

            if " ".join(items).find("Anabaena") > -1:
                print("="*50)
                print(items)
                print(candidates)
                print("="*50)
                        
            # Sort candidates by priority - from longest to shortest
            candidates = sorted(candidates, key=lambda x:-len(x[1]))

            # Finished collecting candidates; match them in NCBI
            # Note: matching all names in a single call is much faster
            ncbiMatches = ncbiTaxa.get_name_translator([x[1] for x in candidates])

            speciesStartItem = None
            speciesLevelTaxon = None

            # Find the first (i.e., longest) match
            for candidateStartItem, candidateName in candidates:
                if candidateName in ncbiMatches:
                    candidateTaxId = ncbiMatches[candidateName][0] 
                    
                    speciesLevelTaxon = getContainingSpeciesLevelTaxon( candidateTaxId )
                    if speciesLevelTaxon is None: # Make sure this taxon is not above species rank
                        continue

                    # select this candidate
                    matchingName = candidateName
                    matchingTaxId = candidateTaxId
                    speciesStartItem = candidateStartItem # remember at what item the selected species starts
                    break
                    
            # Did we find a match for this leaf?
            if not matchingName is None: # match found
                node.name = matchingName
                node.add_features(taxId = matchingTaxId, lineageItems = items[:speciesStartItem])
                if matchingTaxId != speciesLevelTaxon:
                    node.add_features(speciesLevelTaxon=speciesLevelTaxon)

                #print("--"*20)
                #print(matchingName)
                #print(items[:speciesStartItem])
                numMatched += 1
            else: # no match found
                unmatched.append(" ".join(items))
                node.name = "n/a"

            numProcessed += 1
            if (rl()):
                print("(processed %d matched %d)" % (numProcessed, numMatched))

            #if(numProcessed>1000):
            #    break


    # Save unmatched names to file (for examination)
    with open("unmatched_names.txt", "w") as f:
        f.writelines(["%s\n" % x for x in unmatched])

    
    for node in tree.traverse(strategy='postorder'): # children first
        if node.is_leaf():
            continue

        a = []
        for c in node.children:
            try:
                a.append( c.lineageItems )
            except AttributeError:
                pass


        if a:
            out = None
            
            if len(a)==2:
                out = []
                if a[0][-1]=="Anabaena" or a[1][-1]=="Anabaena":
                    print(">"*50)
                    print(a)
                    print(">"*50)
                    
                for u,v in zip(a[0], a[1]):
                    if u==v:
                        out.append(u)
                    else:
                        break
            elif len(a)==1:
                if a[0][-1]=="Anabaena":
                    print(">"*50)
                    print(a)
                    print(">"*50)
                    
                out = a[0]
            else:
                assert(False)

            if out:
                node.add_features(lineageItems = out, testK = out[-1])

                print(">>> %s" % out[-1])
                if out[-1] == "Nostocaceae":
                    print(node.name)
                    print(a)
                    
            #else:
            #    print("*-"*20)
            #    print(node.name)
            #    print(a)
                

        
    # Now we have our annotated reference phylogenetic tree

    # Get our target list of species to appear on the final tree
    taxa = getSpeciesToInclude()
    allNames = ncbiTaxa.get_taxid_translator(taxa)

    
    #taxa.append(1906157)
    #taxa.append(251221)

    #print(ncbiTaxa.get_rank(taxa))
    
    f = set()
    fnodes = []
    notf = set()
    for x in taxa:
        found = tree.search_nodes(taxId=x)
        if found:
            f.add(x)
            fnodes.append(found[0])
        else:
            containingSpeciesLevelTaxon = getContainingSpeciesLevelTaxon(x)
            if x != containingSpeciesLevelTaxon:
                found = tree.search_nodes(taxId=x)

                if not found:
                    found = tree.search_nodes(speciesLevelTaxon=x)
                    
                if found:
                    f.add(x)
                    fnodes.append(found[0])
        #elif ncbiTaxa.get_rank([x])[x] == 'no rank':
        #    parent = ncbiTaxa.get_lineage(x)[-2]
        #    found = tree.search_nodes(taxId=parent)
        #    if found:
        #        f.add(x)
        #        fnodes.append(found[0])

        if x not in f:
            print("--"*50)
            print("TaxId: %d" % x)
            print("Name: %s" % allNames[x])
            containingSpeciesLevelTaxon = getContainingSpeciesLevelTaxon(x)
            print("Species TaxId: %d %s" % (containingSpeciesLevelTaxon, "" if x==containingSpeciesLevelTaxon else "***"))
            
            lineage = ncbiTaxa.get_lineage(x)
            print("Lineage: %s" % lineage)
            names = ncbiTaxa.get_taxid_translator(lineage)
            print(names)

            for y in reversed(lineage):
                name = names[y]
                res = bool(tree.search_nodes(testK = name))
                print("%s: %s" % (name, res))
            notf.add(x)
            
    print("Found (%d): %s" % (len(f), f))
    print(ncbiTaxa.get_rank(list(f)))

    print("Couldn't find (%d): %s" % (len(notf), notf))
    print(ncbiTaxa.get_taxid_translator(list(notf)).values())
    print(len(fnodes))

    tree2 = tree.copy()

    print("Before pruning: %d" % len(tree2))
    if fnodes:
        #tree2.prune(fnodes, preserve_branch_length=True)
        pruneTree(tree2, fnodes)
    print("After pruning: %d" % len(tree2))
    
    

        
    with Display(backend='xvnc') as disp:  # Plotting requires an X session

        print(len(tree))
        tree.render('nmicro_s6.pdf', tree_style=ts, w=5, h=320, units="mm")
        tree.render('nmicro_s6.svg', tree_style=ts, w=5, h=320, units="mm")

        print(len(tree2))
        tree2.render('nmicro_s6_pruned.pdf', tree_style=ts, w=5, h=20, units="mm")
        tree2.render('nmicro_s6_pruned.svg', tree_style=ts, w=5, h=20, units="mm")

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


    ts.aligned_header.add_face(TextFace("Paired"), column=6)

    for cat in ('Psychrophilic', 'Mesophilic', 'Thermophilic', 'Hyperthermophilic'):
        color = temperatureRangeToColor[cat]
        ts.legend.add_face(faces.RectFace( width=100, height=20, fgcolor=color, bgcolor=color, label={"text":cat, "color":"Black", "fontsize":8}), column=1)

    for cat in ('NonHalophilic', 'ModerateHalophilic', 'ExtremeHalophilic'):
        color = salinityToColor[cat]
        ts.legend.add_face(faces.RectFace( width=100, height=20, fgcolor=color, bgcolor=color, label={"text":cat, "color":"Black", "fontsize":8}), column=2)

    for cat in ('Aerobic', 'Facultative', 'Anaerobic', 'Microaerophilic'):
        color = oxygenReqToColor[cat]
        ts.legend.add_face(faces.RectFace( width=100, height=20, fgcolor=color, bgcolor=color, label={"text":cat, "color":"Black", "fontsize":8}), column=3)

        

        
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

        tree.render('alltaxa.pdf', tree_style=ts, w=9, h=38, units="mm")
        tree.render('alltaxa.svg', tree_style=ts, w=9, h=38, units="mm")

        print("------------------ Ignore error message ------------------")
        # Display is about to close; how to tell tree to disconnect cleanly? (to prevent "Client Killed" message...)

    return 0

def standalone():
    import argparse
    from glob import glob
    import os
    
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument("--use-profile-data", type=str, default="")
    args = argsParser.parse_args()

    files = []
    if args.use_profile_data:
        files = [x for x in glob(args.use_profile_data) if os.path.exists(x)]

    print("Loading profile data for %d files..." % len(files))
    (xdata, ydata, ydata_nativeonly, ydata_shuffledonly, labels, groups, filesUsed, biasProfiles, dfProfileCorrs, summaryStatistics) = loadProfileData(files)

    tileGenerator = ProfileDataTileGenerator(files, biasProfiles, dfProfileCorrs)
    
    #return plotSpeciesOnTaxonomicTree()
    return plotSpeciesOnTaxonomicTree(tileFunc=tileGenerator.getProfileTileFunc())
    #return testPlotTree()


if __name__=="__main__":
    import sys
    sys.exit(standalone())


        
    
