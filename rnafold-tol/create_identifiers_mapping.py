from csv import reader
import codecs
import requests
import urllib
from time import sleep
import xml.etree.ElementTree as ET
from data_helpers import getSpeciesName as getSpeciesName_remote
from pcache import pcache
#from data_helpers import SpeciesCDSSource, CDSHelper
#from genome_model import getGenomeModelFromCache

requestDelaySeconds = 0.4    # must be above 1/3 unless an API key is used or other stuff is taking up time (https://www.ncbi.nlm.nih.gov/books/NBK25497/)

# format: (<EXT_IDENTIFIER>, <TAXID>)
testIdentifiers = \
(('lmo1688', 169963), ('lmo0530', 169963), ('lmo0433', 169963), ('lmo0122', 169963), ('lmo0123', 169963), ('lmo0125', 169963), ('lmo0126', 169963), ('lmo0127', 169963), ('lmo0129', 169963), ('lmo1439', 169963), ('lmo1438', 169963), ('lmo1431', 169963), ('lmo1434', 169963), ('lmo1436', 169963), ('lmo2829', 169963), ('lmo1967', 169963), ('lmo0690', 169963), ('lmo0691', 169963), ('lmo0697', 169963), ('lmo0197', 169963), ('lmo0196', 169963), ('lmo1501', 169963), ('lmo1503', 169963), ('lmo0192', 169963), ('lmo0054', 169963), ('lmo0055', 169963), ('lmo0199', 169963), ('lmo0198', 169963), ('lmo0052', 169963), ('lmo0053', 169963), ('lmo0825', 169963), ('lmo0392', 169963), ('lmo0943', 169963), ('lmo1068', 169963), ('lmo1067', 169963), ('lmo1062', 169963), ('lmo1611', 169963), ('lmo1759', 169963), ('lmo0893', 169963), ('lmo1615', 169963), ('lmo0898', 169963)) + \
(('PA3013', 208964), ('PA4241', 208964), ('PA3019', 208964), ('PA5535', 208964), ('PA5234', 208964), ('PA5530', 208964), ('PA4715', 208964), ('PA5532', 208964), ('PA5539', 208964), ('PA2986', 208964)) + \
(('DVU2669', 882), ('DVU2523', 882), ('DVU2529', 882), ('DVU2281', 882), ('DVU0725', 882), ('DVU3274', 882), ('DVU2885', 882), ('DVU3272', 882), ('DVU1453', 882), ('DVU2448', 882), ('DVU3271', 882), ('DVU2444', 882), ('DVU1682', 882), ('DVU2442', 882), ('DVU1684', 882), ('DVU3279', 882), ('DVU2225', 882), ('DVU3379', 882), ('DVU2667', 882), ('DVU1094', 882), ('DVU1095', 882), ('DVU1096', 882), ('DVU1097', 882), ('DVU0538', 882), ('DVU1091', 882), ('DVU2324', 882), ('DVU1093', 882), ('DVU1397', 882), ('DVU0305', 882), ('DVU0278', 882), ('DVU0302', 882), ('DVU1258', 882), ('DVU1257', 882), ('DVU0627', 882), ('DVU0722', 882), ('DVU1251', 882), ('DVU0242', 882), ('DVU2675', 882), ('DVU1470', 882), ('DVU0900', 882), ('DVU2231', 882), ('DVU2230', 882), ('DVU0712', 882), ('DVU3002', 882), ('DVU1953', 882), ('DVU0716', 882), ('DVU0890', 882), ('DVU0659', 882), ('DVU2970', 882), ('DVU2976', 882), ('DVU1618', 882), ('DVU1619', 882), ('DVU0650', 882), ('DVU1615', 882)) + \
(('VNG1366H', 64091), ('VNG2063G', 64091), ('VNG0883H', 64091), ('VNG2452C', 64091), ('VNG1380H', 64091), ('VNG0447H', 64091), ('VNG1562H', 64091), ('VNG2051G', 64091), ('VNG0487H', 64091), ('VNG1426H', 64091), ('VNG1459H', 64091), ('VNG0450C', 64091), ('VNG0405C', 64091), ('VNG2521H', 64091), ('VNG2372G', 64091), ('VNG0583G', 64091), ('VNG0180G', 64091), ('VNG2158G', 64091), ('VNG1311G', 64091), ('VNG0798H', 64091), ('VNG2569H', 64091), ('VNG0326G', 64091), ('VNG2498H', 64091), ('VNG1781C', 64091), ('VNG0874G', 64091), ('VNG1776G', 64091), ('VNG2343G', 64091), ('VNG0690C', 64091), ('VNG2292H', 64091), ('VNG0563G', 64091), ('VNG1306G', 64091), ('VNG2074H', 64091), ('VNG0700G', 64091), ('VNG0259G', 64091), ('VNG2418G', 64091), ('VNG2445C', 64091), ('VNG0992H', 64091), ('VNG0965C', 64091), ('VNG0512G', 64091), ('VNG0623G', 64091), ('VNG2116C', 64091), ('VNG1140G', 64091), ('VNG1998H', 64091), ('VNG0983C', 64091), ('VNG0736G', 64091), ('VNG2267G', 64091), ('VNG0676C', 64091), ('VNG1153G', 64091), ('VNG1621H', 64091), ('VNG2617G', 64091), ('VNG1299C', 64091), ('VNG1097G', 64091), ('VNG1204G', 64091), ('VNG1656H', 64091), ('VNG0869G', 64091), ('VNG2558G', 64091), ('VNG2036G', 64091), ('VNG1219G', 64091), ('VNG2507G', 64091), ('VNG0340C', 64091), ('VNG0718C', 64091), ('VNG2600G', 64091), ('VNG2207H', 64091), ('VNG0796G', 64091), ('VNG1929G', 64091), ('VNG2675C', 64091), ('VNG2586C', 64091), ('VNG0331H', 64091), ('VNG1834G', 64091), ('VNG1698G', 64091), ('VNG1543G', 64091), ('VNG0901G', 64091), ('VNG1149Cm', 64091), ('VNG0835G', 64091), ('VNG1916H', 64091), ('VNG1412H', 64091), ('VNG1518H', 64091), ('VNG1501G', 64091), ('VNG0069H', 64091), ('VNG1618H', 64091), ('VNG1918C', 64091), ('VNG0823G', 64091), ('VNG0021H', 64091), ('VNG1123Gm', 64091), ('VNG1253C', 64091), ('VNG1162H', 64091), ('VNG1532G', 64091), ('VNG1075G', 64091), ('VNG1133G', 64091))



# Local memoization
localCache = {}
def memoize(tag=""):
    def wrapper(wrappedFunc):
        def cachedCall(*args, **kw):
            global localCache

            m = []
            for arg in args:
                m.append(str(arg))

            for ident in sorted(kw.keys()): # iterate over the keys in well-defined order
                m.append( "{}:{}".format( ident, kw[ident])  )

            lookupKey = "{}:{}".format( "!".join(m), tag )
            
            retVal = localCache.get( lookupKey, None )

            if retVal is None:
                # produce a new value
                retVal = wrappedFunc( *args, **kw )

                # store the new value inside the cache
                localCache[ lookupKey ] = retVal

            return retVal
            
        return cachedCall
            
    return wrapper


# Encode return values as "url-encoded" (e.g., 'Hello world' -> 'Hello+world')
def urlencoded():
    def wrapper(wrappedFunc):
        def cachedCall(*args, **kw):
            return urllib.quote_plus( wrappedFunc( *args, **kw) )
        return cachedCall
    return wrapper
            

#@memoize
#@urlencoded
#def getSpeciesName(taxId):
#    return getSpeciesName_remote(taxId)
getSpeciesName = memoize()(urlencoded()(getSpeciesName_remote))  # Cache values, and also convert to url-encoded format

def findProteinUsingExternalIdentifier( searchTerm, taxId ):
    # Manual use of the API:
    # 1)
    #     curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=VNG1075G+AND+\"Halobacterium+salinarum+NRC-1\"[Organism]"
    # 2) Extract first identifier
    # 3)
    #     curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=10580625&rettype=gp&retmode=xml"

    #requestUri = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={}&rettype=gp&retmode=xml".format( genProtAccession )
    requestUri = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term={}+AND+\"{}\"[Organism]&retmode=xml".format( searchTerm, getSpeciesName(taxId))

    sleep(requestDelaySeconds)
    r = requests.get( requestUri )
    if( r.status_code != 200 ):
        raise Exception("Entrez returned HTTP code {}, content {}".format( r.status_code, r.text))
    
    #return r.text # codecs.decode(r.text, 'utf-8', 'ignore')
    root = None
    try:
        root = ET.fromstring(r.text)
        
    except ET.ParseError as e:
        fn = 'create_identifiers_mapping.py.search_document.xml'
        f = open(fn, 'w')
        f.write(xml)
        f.close()
        print("XML error in report; saved to %s" % fn)
        raise e

    elems = root.findall(".//IdList/Id")
    
    if elems:
        return elems[0].text
    else:
        return None

def getProteinIdentifiers( proteinAccession ):
    ret = []
    
    requestUri = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={}&rettype=gp&retmode=xml".format( proteinAccession )

    sleep(requestDelaySeconds)
    r = requests.get( requestUri )
    if( r.status_code != 200 ):
        #raise Exception("Entrez returned HTTP code {}, content {}".format( r.status_code, r.text))
        return ret
    
    #return r.text # codecs.decode(r.text, 'utf-8', 'ignore')
    root = None
    try:
        root = ET.fromstring(r.text)
        
    except ET.ParseError as e:
        fn = 'create_identifiers_mapping.py.fetch_protein.xml'
        f = open(fn, 'w')
        f.write(xml)
        f.close()
        print("XML error in report; saved to %s" % fn)
        #raise e
        return ret

    elems = root.findall("./GBSeq/GBSeq_locus")
    if elems:
        ret.append( elems[0].text )

    elems = root.findall(".//GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier[GBQualifier_name=\"locus_tag\"]/GBQualifier_value")
    if elems:
        ret.append( elems[0].text )

    return ret

def findExternalProteinIdentifiers( term, taxId ):
    # Find the available identifiers, using two steps:
    protAccession = findProteinUsingExternalIdentifier( term, taxId )
    
    additionalIdentifiers = getProteinIdentifiers( protAccession )

    return additionalIdentifiers

@pcache("entrez-proteins-id")
def createMappingForSpeciesProteins( identifiers, taxId ):
    ret = {}
    print("Creating identifiers mapping for taxid={} (N={})...".format( taxId, len(identifiers)))
    for ident in identifiers:
        ret[ ident ] = findExternalProteinIdentifiers( ident, taxId )
        if len( ret ) % 50 == 49:
            print("{:.4}%".format(float(len(ret))/len(identifiers)*100))
    return ret

def testAll():
    from random import shuffle

    testSet = list(testIdentifiers)
    shuffle( testSet )
    
    for term, taxId in testSet:
        additionalIdentifiers = findExternalProteinIdentifiers( term, taxId )
        print("{} -> {}".format(term, additionalIdentifiers))

    return 0

if __name__=="__main__":
    import sys
    sys.exit( testAll() )

    
        
