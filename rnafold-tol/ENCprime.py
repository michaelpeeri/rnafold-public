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
# Annotate ENc prime for species (by running the external program)
from __future__ import print_function
import subprocess
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
import config
from data_helpers import CDSHelper, SpeciesCDSSource, getSpeciesTranslationTable, getSpeciesName, getSpeciesProperty, setSpeciesProperty

ENCprimeCdsCount = "{}/SeqCount".format( config.ENCprimeBasePath )
ENCprimeMain     = "{}/ENCprime".format( config.ENCprimeBasePath )
debugMode = False

def writeSequenceToTempFile_orig(taxId):
    
    allRecords = []
    print("Fetching sequence for taxid={}".format(taxId))
    
    for protId in SpeciesCDSSource(taxId):
        cds = CDSHelper(taxId, protId)
        
        if( cds.length()%3 != 0 ):
            continue
        
        seq = cds.sequence()
        record = SeqRecord( Seq(seq, NucleotideAlphabet), id=protId, description="" )
        allRecords.append( record )
        if( len(allRecords)%1000 == 999 ): print(".")

    fout = NamedTemporaryFile( mode="w", delete=(not debugMode) )
    SeqIO.write( allRecords, fout.name, "fasta")  # write the full sequences into the file

    return (len(allRecords), fout)

def writeSequenceToTempFile(taxId):
    print("Fetching sequence for taxid={}".format(taxId))
    
    allRecords = []
    allCDSs = []
    
    for protId in SpeciesCDSSource(taxId):
        cds = CDSHelper(taxId, protId)
        
        if( cds.length()%3 != 0 ):
            continue
        
        seq = cds.sequence()
        allCDSs.append( seq )

        if( len(allCDSs)%1000 == 999 ): print(".")

    record = SeqRecord( Seq( ''.join(allCDSs), NucleotideAlphabet), id="allCDSs", description="" )
    allRecords.append( record )
    fout = NamedTemporaryFile( mode="w", delete=(not debugMode) )
    SeqIO.write( allRecords, fout.name, "fasta")  # write the full sequences into the file

    return (len(allCDSs), fout)


def createCodonCounts(fastaPath, fastaCount):
    #fout = NamedTemporaryFile(mode='r')

    print("Creating codon counts...")
    #SeqCount -c ExSeqs.fasta 9
    print(" ".join((ENCprimeCdsCount, "-c", fastaPath, str(fastaCount))))
    out = subprocess.check_output((ENCprimeCdsCount, "-c", fastaPath, str(fastaCount)), shell=False)
    #print(out)
    
def createNucleotideCounts(fastaPath, fastaCount):
    #fout = NamedTemporaryFile(mode='r')

    print("Creating codon counts...")
    #SeqCount -c ExSeqs.fasta 9
    print(" ".join((ENCprimeCdsCount, "-n", fastaPath, str(fastaCount))))
    out = subprocess.check_output((ENCprimeCdsCount, "-n", fastaPath, str(fastaCount)), shell=False)
    #print(out)

def createEncPrimeReport(fastaPath, geneticCode):
    fout = NamedTemporaryFile(mode='r', delete=(not debugMode) )

    print("Creating codon counts...")
    #SeqCount -c ExSeqs.fasta 9
    print((ENCprimeMain, "{}.codcnt".format(fastaPath), "{}.acgtfreq".format(fastaPath), str(geneticCode), fout.name, "0", "-q"))
    out = subprocess.check_output((ENCprimeMain, "{}.codcnt".format(fastaPath), "{}.acgtfreq".format(fastaPath), str(geneticCode), fout.name, "0", "-q"), shell=False)
    
    print(out)

    for line in fout.readlines():
        print("= {}".format(line))
        if line.find("Totals: ")==0:
            fields = line.split(" ")
            ENc       = float(fields[1])
            ENc_prime = float(fields[2])
            print("ENc: {}\tENc': {}".format(ENc, ENc_prime))

    return (ENc, ENc_prime)


def calculateENcPrimeForSpecies(taxId, orig=False):
    geneticCode = getSpeciesTranslationTable( taxId )

    if orig:
        cdsCount, fastaFile = writeSequenceToTempFile_orig( taxId )
    else:
        cdsCount, fastaFile = writeSequenceToTempFile( taxId )
        
    createCodonCounts(fastaFile.name, cdsCount)
    createNucleotideCounts(fastaFile.name, cdsCount)
    print("Genomic GC%: {}".format(getSpeciesProperty(taxId, 'gc-content')))
    
    return createEncPrimeReport(fastaFile.name, geneticCode)

def annotateENcPrime(taxId, overwrite=False):
    encPropValue      = getSpeciesProperty(taxId, 'ENc')
    encPrimePropValue = getSpeciesProperty(taxId, 'ENc-prime')

    # TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY #
    ###return (taxId, 0.0, 1.0, False)  # return old values (last value indicates this value are old)
    # TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY #

    
    if (encPropValue[0] is None) or (encPrimePropValue[0] is None) or overwrite:
        ENc, ENc_prime = calculateENcPrimeForSpecies(taxId)
        assert(ENc      <75.0 and ENc      >10.0)   # The actual extreme values for ENc are not clear to me, but let's do a sanity check
        assert(ENc_prime<75.0 and ENc_prime>10.0)
        
        setSpeciesProperty( taxId, 'ENc',       str(ENc),       "ENCprime (custom version)" )
        setSpeciesProperty( taxId, 'ENc-prime', str(ENc_prime), "ENCprime (custom version)" )
        
    else:
        return (taxId, encPropValue[0], encPrimePropValue[0], False)  # return old values (last value indicates this value are old)
        
    return (taxId, ENc, ENc_prime, True)  # return values (last value indicates new values)
    

if __name__=="__main__":
    import sys
    taxId=int(sys.argv[1])
    print("Name: {}".format(getSpeciesName(taxId)))

    print(annotateENcPrime(taxId))
    
    print( calculateENcPrimeForSpecies(taxId, orig=True) )
    print( calculateENcPrimeForSpecies(taxId, orig=False) )
    

