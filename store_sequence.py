# Input - CDS JSON file for a specific organism
# Build CDS entries from the JSON file
import sys
import urllib2
from urllib import urlencode
from HTMLParser import HTMLParser
import codecs
import redis
import ijson
import config

r = redis.StrictRedis(host=config.host, port=config.port, db=config.db)


#f = open('/home/michael/rnafold/Creinherdtii.json', 'r')
#f = open('/home/michael/rnafold/Ptricornutum.json', 'r')
f = open(sys.argv[1], 'r')

visitedProteinIds = set()

# Parse NCBI taxonomy results; capture the returned taxonomy-id
class NCBIParser(HTMLParser):
    #def __init__(self):
    #    self.taxId = 0
    #    super(NCBIParser, self).__init__(self)

    #def handle_starttag(self, tag, attrs):
    #    self.lasttag = tag

    def handle_data(self, data):
        if(self.get_starttag_text()=='<pre>'):
            #assert(self.taxId == 0)
            self.taxId = int(data)

# Convert binomial species name to an NCBI taxonomy-id, using the NCBI server.
def getTaxonomyId(binomialName):
    response = urllib2.urlopen('http://ncbi.nlm.nih.gov/taxonomy/?%s' % urlencode((('term', binomialName),('report','taxid'),('format','text'))))
    assert(response.getcode()==200)
    #assert(response.info.getheader('Content-Type').startswith('text/html'))
    parser = NCBIParser();
    parser.feed(codecs.decode(response.read()))
    taxId = parser.taxId
    assert(taxId > 0)
    return taxId



species = next(ijson.items(f, 'Species'))
taxId = getTaxonomyId(species)

# Store a species entry for this species
r.set('species:taxid:%d:name' % (taxId,), species)
r.set('species:name:%s:taxid' % (species,), taxId)


f.seek(0) # reset the file for another parsing iteratin; Is there a nicer way to do this?

cdsCount = 0
objects = ijson.items(f, 'CDSlist.item')
for o in objects:
    proteinId = o['attributes']['protein_id']
    # verify there are no duplicates entries
    assert(proteinId not in visitedProteinIds)
    visitedProteinIds.add(proteinId)

    # Build a CDS entry for this protein
    r.set('CDS:taxid:%d:protid:%s:seq' % (taxId, proteinId), o['CDS'])
    r.set('CDS:taxid:%d:protid:%s:num-exons' % (taxId, proteinId), int(o['numExons']))
    r.set('CDS:taxid:%d:protid:%s:strand' % (taxId, proteinId), o['strand'])
    r.set('CDS:taxid:%d:protid:%s:reference' % (taxId, proteinId), o['reference'])
    if('partial' in o['attributes'] and o['attributes']['partial'].lower()=='true' ):
        r.set('CDS:taxid:%d:protid:%s:partial' % (taxId, proteinId), True);

    r.sadd('species:taxid:%d:CDS' % (taxId,), proteinId)
    cdsCount += 1

print("Processed %d CDS entries" % (cdsCount,))
print(r.scard('species:taxid:%d:CDS' % (taxId,)))
