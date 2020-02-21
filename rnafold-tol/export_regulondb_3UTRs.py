# RNAFold - Analyze mRNA folding bias (Local Folding Energy) based on randomizations.
# Copyright (C) 2016-2020 Michael Peeri
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
