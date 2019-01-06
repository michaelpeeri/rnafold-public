import requests, sys
from time import sleep
import pygraph.classes.digraph

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

def debugSource():
    for i in (420778, 1069680):
        yield i

G = pygraph.classes.digraph.digraph()

for taxid in debugSource():
    ancestors = et.getAncestors(taxid)

    for node in ancestors:
        if not G.has_node(node):
            G.add_node(node)

    for u,v in zip(ancestors, ancestors[1:]):
        if G.has_edge((u,v)):
            
        else:
            G.add_edge((u,v), attrs=[('count',1)])
            
print(G)

#a.add_node('A')
#>>> a.add_node('B')
#>>> a.add_edge('A','B')
