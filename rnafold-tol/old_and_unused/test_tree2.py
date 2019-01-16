import requests, sys
from time import sleep

class EnsembleTaxonomy(object):
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

        for level in decoded:
            ancestors.append(level['name'])

        ancestors.append(decoded[-1]['parent']['name'])

        return ancestors

et = EnsembleTaxonomy()


testTaxons = (1041607, 1069680, 1116230, 1125630, 1183438, 1223560, 1287680, 1341181, 163003, 195522, 208964, 223926, 237561, 243274, 251221, 269084, 2769, 284811, 295405, 296543, 3055, 336722, 339860, 391623, 400667, 402612, 420247, 420778, 436017, 44056, 45670, 4781, 511145, 556484, 559292, 633149, 63737, 65357, 695850, 93061, 96563, 997884)

def debugSource():
    for i in testTaxons:
        yield i

def stdinSource():
    pass


class Node(object):
    def __init__(self, name, parent=None):
        self._name = name
        self._parent = parent
        self._count = 0

    def __repr__(self):
        if( self._parent is None ):
            return "Node(name='%s', count=%d, parent=None)" % (self._name, self._count)
        else:
            return "Node(name='%s', count=%d, parent='%s')" % (self._name, self._count, self._parent.name())

    def increment(self):
        self._count += 1

    def name(self):
        return self._name

    def parent(self):
        return self._parent

    def assignParent(self, parent):
        if( self._parent is None):
            self._parent = parent
        elif self._parent == parent:
            return
        else:
            raise Exception("mismatched parents")
        
G = {}

for taxid in debugSource():
    ancestors = et.getAncestors(taxid)

    for node, _parent in zip(ancestors, ancestors[1:]):
        print(node, _parent)

        if not _parent in G:
            G[_parent] = Node(_parent)
            
        if node in G:
            G[node].increment()
            G[node].assignParent(G[_parent])
        else:
            assert(G[_parent] is not None)
            G[node] = Node(node, parent=G[_parent])
            G[node].increment()
            
print(G)

#a.add_node('A')
#>>> a.add_node('B')
#>>> a.add_edge('A','B')
