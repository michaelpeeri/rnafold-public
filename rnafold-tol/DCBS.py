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
# Annotate DCBS for species (by running external MATLAB code)
from __future__ import print_function
import subprocess
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
import config
from data_helpers import CDSHelper, SpeciesCDSSource, getSpeciesTranslationTable, getSpeciesName, getSpeciesProperty, setSpeciesProperty

debugMode = True

def writeSequenceToTempFile(taxId):
    
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


def runDCBScodeInMatlab(fastaPath, cdsCount):
    fout = NamedTemporaryFile(mode='r', delete=(not debugMode) )

    print("Calculating DCBS...")
    # /usr/local.cc/bin/matlab -nodisplay -nodesktop -nosplash -singleCompThread -nojvm -r "cds=fastaread('/var/tmp/pbs.9025895.power8.tau.ac.il/tmpIzkenm');csvwrite('/var/tmp/pbs.9025895.power8.tau.ac.il/tmpIzkenm.csv',Compute_DCBS({cds.Sequence}));quit()"
    cmdline = ((config.MatlabPath, "-nodisplay", "-nodesktop", "-nosplash", "-singleCompThread", "-nojvm", "-r",
                "cds=fastaread('{}');csvwrite('{}',Compute_DCBS({{cds.Sequence}}));quit()".format(fastaPath, fout.name)))
    print(" ".join(cmdline))
    out = subprocess.call(cmdline, shell=False)

    lastVal = None
    lineNum = 0
    for line in fout.readlines():
        lastVal = float(line)
        lineNum += 1

    if lineNum != cdsCount+1:
        raise Exception("Expected to find DCBS values for {} proteins; found {} instead".format(cdsCount, lineNum))

    return lastVal
        

def calcDCBS(taxId):
    cdsCount, fastaFile = writeSequenceToTempFile( taxId )
    return runDCBScodeInMatlab( fastaFile.name, cdsCount )
        
def run():
    import sys
    taxId = int(sys.argv[1])
    print(calcDCBS(taxId))


def annotateDCBS(taxId, overwrite=False):
    dcbsPropValue      = getSpeciesProperty(taxId, 'DCBS-geomean')

    # TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY #
    ## return (taxId, 0.0, False)  # return old values (last value indicates this value are old)
    # TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY #

    
    if (dcbsPropValue[0] is None) or overwrite:
        DCBS = float(calcDCBS(taxId))
        #assert(ENc_prime<75.0 and ENc_prime>10.0)
        
        setSpeciesProperty( taxId, 'DCBS-geomean',   str(DCBS),       "DCBS (matlab, Renana)")
        
    else:
        return (taxId, dcbsPropValue[0], False)  # return old values (last value indicates this value are old)
        
    return (taxId, DCBS, True)  # return values (last value indicates new values)

