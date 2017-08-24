import re
from ete3 import NCBITaxa, Tree, TreeStyle, NodeStyle, TextFace, faces, AttrFace, PhyloTree
from rate_limit import RateLimit
from ncbi_taxa import ncbiTaxa



# configuration
nmicrobiol201648_s6 = "./data/nmicrobiol201648-s6.txt"
nmicrobiol201648_s8 = "./data/nmicrobiol201648-s8.txt"
itol_newick         = "./data/itol_newick.txt"


taxonItemsToIgnore = frozenset(("Unclassified", "Incertae", "Candidate", "candidate", "Candidatus", "incertae", "sedis"))

reNodeSupport = re.compile('[[][^]]*[]]')
def stripTreeInternalSupport(tree):
    return re.sub(reNodeSupport, '', tree)

rl = RateLimit(15)


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
            


        
def pruneReferenceTree_Nmicrobiol201648(taxa):
    treefmt = None
    
    with open(nmicrobiol201648_s6, "r") as f:
    #with open(itol_newick, "r") as f:
        treefmt = f.read()

    # Strip internal support, and parse the resulting string
    tree = PhyloTree(stripTreeInternalSupport(treefmt))

    #ts.show_leaf_name = True
    #ts.layout_fn = nodeLayoutWithTaxonomicNames

    unmatched = []
    matched = []

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
            #
            # See also: FuzzyMatch in R 'evobiR' package (https://cran.r-project.org/web/packages/evobiR/index.html)
            #
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
                
                candidates.append( (i,candidate,-1) )

                # For each suffix, also consider its prefixes as candidates. This is bacause the species names sometimes include additional strain identifiers that are not included in NCBI
                if i<=len(items)-2 and items[i] and items[i+1] and items[i][0].isupper() and items[i+1][0].islower(): # require that the first item start with a capital, and the second item start with a lower-case letter (e.g., "Eschericia coli")
                    for j in range(i+1,len(items)-1):
                        candidate = " ".join(items[i:j+1])

                        candidates.append( (i,candidate,j+1) )
                        

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
            for candidateStartItem, candidateName, _ in candidates:
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
                #print("<%s>" % node.name)

                lineageItems = [x for x in items[:speciesStartItem] if not x in taxonItemsToIgnore]
                
                node.add_features(taxId = matchingTaxId, lineageItems = lineageItems)
                if matchingTaxId != speciesLevelTaxon:
                    node.add_features(speciesLevelTaxon=speciesLevelTaxon)

                matched.append("%s [%d]" % (matchingName, matchingTaxId))


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
        f.writelines(["%s\n" % x for x in sorted(unmatched)])

    # Save unmatched names to file (for examination)
    with open("matched_names.txt", "w") as f:
        f.writelines(["%s\n" % x for x in sorted(matched)])
        

    # Try to annotate non-leaf nodes with common taxonomic group
    for node in tree.traverse(strategy='postorder'): # children first
        if node.is_leaf():
            continue

        a = []
        for c in node.children:
            try:
                a.append( c.lineageItems )
            except AttributeError:
                pass


        print(">>" * 20)
        print(len(a))
        print(a)
        
        if a:
            out = None
            
            if len(a)==2:
                out = []
                #if a[0][-1]=="Anabaena" or a[1][-1]=="Anabaena":
                #    print(">"*50)
                #    print(a)
                #    print(">"*50)
                    
                for u,v in zip(a[0], a[1]):
                    if u==v:
                        out.append(u)
                    else:
                        break
            elif len(a)==1:
                #if a[0][-1]=="Anabaena":
                #    print(">"*50)
                #    print(a)
                #    print(">"*50)
                    
                out = a[0]
            else:
                assert(False)

            if out:
                #print("out = %s" % out)
                node.add_features(lineageItems = out, testK = out[-1])

                print(">>> %s" % out[-1])
                
            #    #if out[-1] == "Nostocaceae":
            #        #print(node.name)
            #        #print(a)
            #        #pass
                    
            #else:
            #    print("*-"*20)
            #    print(node.name)
            #    print(a)
                

        
    # Now we have our annotated reference phylogenetic tree

    # Get our target list of species to appear on the final tree
    #taxa = getSpeciesToInclude()
    allNames = ncbiTaxa.get_taxid_translator(taxa)

    
    #taxa.append(1906157)
    #taxa.append(251221)

    #print(ncbiTaxa.get_rank(taxa))

    for x in (45157,4896,44056):
        print(x)
        print(tree.search_nodes(taxId=x))
        
    
    f = set()
    fnodes = []
    notf = set()
    for x in taxa:
        print("=="*20)
        print("Searching for %d" % x)
        found = tree.search_nodes(taxId=x)
        
        if found:
            f.add(x)
            fnodes.append(found[0])
            print("Exact match found")
        else:
            containingSpeciesLevelTaxon = getContainingSpeciesLevelTaxon(x)
            if x != containingSpeciesLevelTaxon:
                found = tree.search_nodes(taxId=containingSpeciesLevelTaxon)

                if not found:
                    found = tree.search_nodes(speciesLevelTaxon=containingSpeciesLevelTaxon)

                # TODO - CONTINUE HERE
                #if not found:
                #    found = tree.search_nodes(
                    
                if found:
                    f.add(x)
                    fnodes.append(found[0])
                    print("Found")
                else:
                    print("Not found at all...")
                        
        #elif ncbiTaxa.get_rank([x])[x] == 'no rank':
        #    parent = ncbiTaxa.get_lineage(x)[-2]
        #    found = tree.search_nodes(taxId=parent)
        #    if found:
        #        f.add(x)
        #        fnodes.append(found[0])
                

        if x not in f:
            print("--"*50)
            print("TaxId not found: %d" % x)
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

    return (tree, tree2)

    
