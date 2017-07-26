from data_helpers import allSpeciesSource
import re
from pyvirtualdisplay.display import Display  # use Xvnc to provide a headless X session (required by ete for plotting)
#from easyprocess import EasyProcess
from ete3 import NCBITaxa, Tree, TreeStyle, TextFace, faces, AttrFace, PhyloTree
from rate_limit import RateLimit


# configuration
maxDisplayNameLength = 27
speciesToExclude = frozenset((99999, 9876543))
speciesToInclude = frozenset()
nmicrobiol201648_s6 = "./data/nmicrobiol201648-s6.txt"
nmicrobiol201648_s8 = "./data/nmicrobiol201648-s8.txt"

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
    


def nodeLayoutWithTaxonomicNames(node):
    level = len(node.get_ancestors())
    if( level==1 or level==2 or level==3):
        name = AttrFace("name", fsize=9)
        #faces.add_face_to_node(n, node, 0, position="float")
        faces.add_face_to_node(name, node, column=0)

        tf = TextFace(len(node.get_leaves()), fsize=12)
        faces.add_face_to_node(tf, node, column=1, position="float")
        
    elif( level==0 ):
        tf = TextFace(len(node.get_leaves()), fsize=12)
        faces.add_face_to_node(tf, node, column=1, position="float")

    elif node.is_leaf():
        name = AttrFace("name", fsize=9, fstyle="italic")
        faces.add_face_to_node(name, node, column=0)

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
        print(len(childTaxIds))

        if not taxIds.intersection(childTaxIds):
            node.detach()
            print(".")
            
    


def testPlotTree():
    treefmt = None
    with open(nmicrobiol201648_s6, "r") as f:
        treefmt = f.read()
            
    tree = PhyloTree(stripTreeInternalSupport(treefmt))

    ts = TreeStyle()
    ts.show_leaf_name = True
    #ts.layout_fn = nodeLayoutWithTaxonomicNames

    unmatched = []

    xxx = set([1287680, 115713, 203267, 188937, 4781, 187420, 243230, 130081, 227882, 228908, 227377, 224308, 5693, 345663, 208964, 224325, 1116230, 243273, 213585, 64091, 45670, 1069680, 1397361, 280699, 1047168, 284811, 284812, 46234, 418459, 214684, 262768, 243365, 273063, 511145, 176299, 272557, 272558, 402612, 283166, 223926, 163003, 559292, 1041607, 1183438, 2769, 122586, 273116, 593117, 192222, 1574623, 243159, 160490, 212717, 272623, 272631, 272632, 63737, 272634, 1341181, 1125630, 99287, 27923, 400667, 269084, 257314, 96563, 300852, 4927, 381764, 242507, 65357, 104782, 336722, 190304, 882, 347515, 353152, 83332, 93061, 194439, 1223560, 267671, 196164, 1245935, 449447, 420778, 195522, 556484, 5061, 391623, 70601, 85962, 272844, 259536, 272633, 220668, 169963, 295405, 237561, 407035, 997884, 1432061, 1010810, 562, 1010800])
        
    numProcessed = 0
    numMatched = 0
    # Annotate tree nodes
    for node in tree.traverse():
        if node.is_leaf():
            n = node.name.split("_")
            
            f = None
            taxid = None
            #print("---------------------------")
            #print(n)

            for i in range(len(n)-1, 0, -1):
                candidate = " ".join(n[i:])

                if len(candidate) < 10:
                    continue
                
                if len(candidate) > 70:
                    break

                #print(candidate)
                # try to match full string
                x = ncbiTaxa.get_name_translator([candidate])
                if x:
                    f = candidate
                    taxid = x[candidate][0]
                else:
                    # Failed to find full string, try to match prefixes
                    if i<=len(n)-2 and n[i] and n[i+1] and n[i][0].isupper() and n[i+1][0].islower():
                        #print("----------------")
                        #print(candidate)
                        for j in range(i+1,len(n)-1):
                            cand2 = " ".join(n[i:j+1])
                            #print(cand2)
                            x = ncbiTaxa.get_name_translator([cand2])
                            if x:
                                f = cand2
                                taxid = x[cand2][0]
                    
                # continue to try to find a longer match

            if not f is None:
                node.name = f
                if taxid in xxx:
                    print(f)
                if f.startswith("Escher") or f.startswith("Plasmod") or f.startswith("Aquifex")  or f.startswith("Pythium")  or f.startswith("Dictyo")  or f.startswith("Cornyn"):
                    print(f)
                    print(taxid)
                    
                node.add_features(taxId=taxid)
                numMatched += 1
            else:
                unmatched.append(" ".join(n))
                node.name = "n/a"

            numProcessed += 1
            if (rl()):
                print("(processed %d matched %d)" % (numProcessed, numMatched))

            #if(numProcessed>1000):
            #    break

            #assert(len(n) > len(p) + 3)
            #assert(n[:len(p)]==p)
            #node.name = " ".join(n.split("_")[-3:-1])[:-25]
            #print(n)
            #print(" ".join(n.split("_")[-5:]))
            #print(n[len(p):])

        
        #taxId = int(node.name)
        #binomicName = ncbiTaxa.get_taxid_translator([taxId])[taxId]  # There has to be an easier way to look up names...
        #node.add_features(displayName=binomicName, taxId=taxId)
        #if( len(binomicName) > maxDisplayNameLength):
        #    node.name = binomicName[:maxDisplayNameLength-1].rstrip() + '...'
        #else:
        #    node.name = binomicName

    with open("unmatched_names.txt", "w") as f:
        f.writelines(["%s\n" % x for x in unmatched])

    taxa = getSpeciesToInclude()
    allNames = ncbiTaxa.get_taxid_translator(taxa)

    
    taxa.append(1906157)
    taxa.append(251221)

    print(ncbiTaxa.get_rank(taxa))
    
    f = set()
    fnodes = []
    notf = set()
    for x in taxa:
        found = tree.search_nodes(taxId=x)
        if found:
            f.add(x)
            fnodes.append(found[0])
        else:
            notf.add(x)
    print("Found (%d): %s" % (len(f), f))
    print("Couldn't find (%d): %s" % (len(notf), notf))
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
        tree2.render('nmicro_s6_pruned.pdf', tree_style=ts, w=5, h=40, units="mm")
        tree2.render('nmicro_s6_pruned.svg', tree_style=ts, w=5, h=40, units="mm")

        # Display is about to close; how to tell tree to disconnect cleanly? (to prevent "Client Killed" message...)
        
    return 0
        
        
"""
Plot "statistical" tree, with species names and counts 
This tree should illustrate the included species
"""
def plotSpeciesOnTaxonomicTree():
    taxa = getSpeciesToInclude()

    # Get the smallest "taxonomic" (i.e., n-ary) tree that includes all specified species
    # See: http://etetoolkit.org/docs/3.0/tutorial/tutorial_ncbitaxonomy.html
    tree = ncbiTaxa.get_topology(taxa, intermediate_nodes=True)
    
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = nodeLayoutWithTaxonomicNames

    # Annotate tree nodes
    for node in tree.traverse():
        taxId = int(node.name)
        binomicName = ncbiTaxa.get_taxid_translator([taxId])[taxId]  # There has to be an easier way to look up names...
        node.add_features(displayName=binomicName, taxId=taxId)
        if( len(binomicName) > maxDisplayNameLength):
            node.name = binomicName[:maxDisplayNameLength-1].rstrip() + '...'
        else:
            node.name = binomicName


    with Display(backend='xvnc') as disp:  # Plotting requires an X session

        tree.render('alltaxa.pdf', tree_style=ts, w=9, h=33, units="mm")
        tree.render('alltaxa.svg', tree_style=ts, w=9, h=33, units="mm")

        # Display is about to close; how to tell tree to disconnect cleanly? (to prevent "Client Killed" message...)

    return 0

if __name__=="__main__":
    import sys
    #sys.exit(plotSpeciesOnTaxonomicTree())
    sys.exit( testPlotTree() )


        
    
