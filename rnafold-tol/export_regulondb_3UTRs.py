import csv
import re

f1                 = "./data/regulondb/UTR_5_3_sequence.txt"
id_conversion_file = "./data/Ensembl/Ecoli/identifiers.tsv"

def removeSuffix(ident):  # convert "abcX-1" -> "abcX"
    if ident[-2]=='-' and ident[-1].isdigit():
        return ident[:-2]
    else:
        return ident
    
def getIdentifiersMapping():
    ret = {}
    with open(id_conversion_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            assert(len(row)==3)
            ret[row[1]]               = row[0]
            ret[removeSuffix(row[2])] = row[0]
    return ret

reGeneNameAndPosition = re.compile("""(['\w-]+)[(][^)]+[)]""")
reGeneCoords =          re.compile("""(\S+)[(](\d+)[,](\d+)[)]""")
reUTRCoords  =          re.compile("""(\d+)-(\d+)""")
def readUTRsTable():
    ret = []
    numSkipped = 0
    identMap = getIdentifiersMapping()
    
    with open(f1) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            if len(row)==1 and row[0][0]=="#": continue
            lastGeneName = reGeneNameAndPosition.match( row[6] ).group(1)
            assert( lastGeneName )

            if lastGeneName in identMap:
                translatedName = identMap[lastGeneName]

                #firstGeneInfo = reGeneCoords.match( row[5] ).groups()
                #lastGeneInfo  = reGeneCoords.match( row[6] ).groups()
                #UTRinfo       = reUTRCoords.match(  row[8] ).groups()

                if len(row) > 10:
                    utr3          = row[10]

                ret.append( (translatedName, lastGeneName, utr3) )
            
    print("Skipped: {}".format(numSkipped))
    return ret

if __name__=="__main__":
    import sys
    data = readUTRsTable()

    for line in data:
        print(">{} {}".format( line[0], line[1] ))
        print(line[2])
    
    sys.exit(0)
