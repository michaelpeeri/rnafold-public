import requests, sys
from time import sleep
from ete3 import TreeNode, TreeStyle, NodeStyle, TextFace, faces, AttrFace

class EnsembleTaxonomy(object):
    # API documentation:
    # http://rest.ensembl.org/documentation/info/taxonomy_classification
    server = "http://rest.ensembl.org"
    ext = "/taxonomy/classification/%d?"
    def __init__(self):
        pass

    def getAncestors(self, taxid):
        r = requests.get(
            EnsembleTaxonomy.server + EnsembleTaxonomy.ext % taxid,
            headers={ "Content-Type" : "application/json"})
 
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        # Rate-limit requests
        sleep(0.5)
 
        decoded = r.json()

        ancestors = []

        taxidName = decoded[0]['children'][0]['scientific_name']
        ancestors.append(taxidName)

        for level in decoded:
            ancestors.append(level['name'])

        ancestors.append(decoded[-1]['parent']['name'])

        #print('--------------------------------------------------')
        #print(ancestors)

        return ancestors

et = EnsembleTaxonomy()


#testTaxons = (1041607, 1069680, 1116230, 1125630, 1183438, 1223560, 1287680, 130081, 1341181, 163003, 195522, 208964, 223926, 237561, 243274, 251221, 269084, 2769, 280699, 284811, 295405, 296543, 3055, 336722, 339860, 391623, 400667, 402612, 420247, 420778, 431895, 436017, 44056, 45670, 4781, 511145, 556484,  559292, 633149, 63737, 65357, 695850,  93061, 946362, 96563, 997884)
testTaxons = (1041607, 1069680, 1116230, 1125630, 1183438, 1223560, 1287680, 130081, 1341181, 163003, 195522, 208964, 223926, 237561, 243274, 251221, 269084, 2769, 280699, 284811, 295405, 296543, 3055, 336722, 339860, 391623, 400667, 400682, 402612, 420247, 420778, 431895, 436017, 44056, 45670, 4781, 511145, 556484, 559292, 633149, 63737, 65357, 6669, 695850, 93061, 946362, 96563, 997884)
print("Displaying %d taxons" % len(testTaxons))

def debugSource():
    for i in testTaxons:
        yield i

def stdinSource():
    pass


# class Node(object):
#     def __init__(self, name, parent=None):
#         self._name = name
#         self._parent = parent
#         self._count = 0

#     def __repr__(self):
#         if( self._parent is None ):
#             return "Node(name='%s', count=%d, parent=None)" % (self._name, self._count)
#         else:
#             return "Node(name='%s', count=%d, parent='%s')" % (self._name, self._count, self._parent.name())

#     def increment(self):
#         self._count += 1

#     def name(self):
#         return self._name

#     def parent(self):
#         return self._parent

#     def assignParent(self, parent):
#         if( self._parent is None):
#             self._parent = parent
#         elif self._parent == parent:
#             return
#         else:
#             raise Exception("mismatched parents")


G = TreeNode(name=u'cellular organisms')
nodes = {u'cellular organisms':G}
nstyle = NodeStyle()
nstyle['shape'] = 'circle'
nstyle['size'] = 3


def layout(node):
    #print(node)
    if(len(node.get_ancestors()) < 4):
        print(node.name)
        n = AttrFace("name", fsize=9)
        n.margin_top = 10
        n.margin_bottom = 0
        n.margin_left= 10
        faces.add_face_to_node(n, node, 0, position="float")
        tf = TextFace(len(node.get_leaves()), fsize=12)
        faces.add_face_to_node(tf, node, 0, position="float")
    

for taxid in debugSource():
    ancestors = et.getAncestors(taxid)

    for node, parent in reversed(zip(ancestors, ancestors[1:])):
        #print(node, _parent)
        print(node)

        if parent in nodes:
            #parent = G.search_nodes(name = _parent)
            parentNode = nodes[parent]
        else:
            parentNode = TreeNode(name = parent)
            parentNode.set_style(nstyle)
            #parentNode.add_face(TextFace(_parent), column=0, position="aligned")
            #faces.add_face_to_node(TextFace(parent), parentNode, 0, position="aligned")
            nodes[parent] = parentNode

        #child = G.search_nodes(name = node)
        if node in nodes:
            childNode = nodes[node]
        else:
            childNode = parentNode.add_child(name = node)
            childNode.set_style(nstyle)
            #childNode.add_face(TextFace(node), column=0, position="aligned")
            #faces.add_face_to_node(TextFace(node), childNode, 0, position="aligned")
            nodes[node] = childNode
            
print(G)
#print(nodes)

for n in G.traverse():
    if( len(n.get_ancestors()) >= 4 ):
        n.dist = 0.1
    else:
        n.dist = 1.0

ts = TreeStyle()
ts.mode = 'c'  # circular
ts.show_leaf_name = True
ts.show_branch_length = False
ts.layout_fn = layout
ts.arc_start = -90
ts.arc_span = 180
ts.scale = 120
ts.branch_vertical_margin = 20

G.render("tree.pdf", tree_style=ts)
