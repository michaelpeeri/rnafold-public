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
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
import _distributed
from data_helpers import allSpeciesSource, CDSHelper, SpeciesCDSSource, getSpeciesProperty, setSpeciesProperty
from codonw import readCodonw

debugMode = False


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

    record = SeqRecord( Seq(''.join(allCDSs), NucleotideAlphabet), id="allCDSs", description="" )
    allRecords.append( record )

    fout = NamedTemporaryFile( mode="w", delete=(not debugMode) )
    SeqIO.write( allRecords, fout.name, "fasta")  # write the full sequences into the file

    return (len(allRecords), fout)


#def annotateGenomicCUB():#
#    fFullSeqs = NamedTemporaryFile(mode="w")         # create a temporary file
#    SeqIO.write( fullSeqs, fFullSeqs.name, "fasta")  # write the full sequences into the file
#    dfCodonw = readCodonw( fFullSeqs.name )          # run codonw and get the gene-level results


def calculateGenomeLevelCUBmeasures(taxId):
    (cdsCount, cdsfile) = writeSequenceToTempFile(taxId)
    cubDf = readCodonw( cdsfile.name )          # run codonw and get the gene-level results
    return cubDf

    
def annotateCUBmeasures(taxId, overwrite=False):
    caiPropValue      = getSpeciesProperty(taxId, 'genomic-CAI')

    # TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY #
    ###return (taxId, 0.0, 1.0, False)  # return old values (last value indicates this value are old)
    # TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY ### TESTING ONLY #

    
    if (caiPropValue[0] is None) or overwrite:
        cubDf = calculateGenomeLevelCUBmeasures(taxId)

        print(cubDf)
        #print(cubDf.index)
        #print(cubDf.iloc[0].at['CAI'])
        
        CAI = cubDf.iloc[0].at['CAI']
        CBI = cubDf.iloc[0].at['CBI']
        Fop = cubDf.iloc[0].at['Fop']
        Nc  = cubDf.iloc[0].at['Nc']

        assert(CAI <  1.0 and CAI>  0.0)
        assert(CBI <  1.0 and CBI> -0.5)
        assert(Fop <  1.0 and Fop>  0.0)
        assert(Nc  < 75.0 and Nc > 10.0)   # The actual extreme values for ENc are not clear to me, but let's do a sanity check
        print(CAI, CBI, Fop, Nc)
        
        setSpeciesProperty( taxId, 'genomic-CAI',       "{:.4}".format(CAI),       "codonw 1.4.4" )
        setSpeciesProperty( taxId, 'genomic-CBI',       "{:.4}".format(CBI),       "codonw 1.4.4" )
        setSpeciesProperty( taxId, 'genomic-Fop',       "{:.4}".format(Fop),       "codonw 1.4.4" )
        setSpeciesProperty( taxId, 'genomic-Nc-codonw', "{:.4}".format(Nc),        "codonw 1.4.4" )
        
    else:
        return (taxId, caiPropValue[0], False)  # return old values (last value indicates this value are old)
        
    return (taxId, CAI, True)  # return values (last value indicates new values)
    

def runDistributed():
    import _distributed
    import dask

    scheduler = _distributed.open()
    delayedCalls = []
    
    for taxId in allSpeciesSource():
        call = dask.delayed( annotateCUBmeasures )(taxId)
        delayedCalls.append( call )

    print("Starting %d calls..." % len(delayedCalls))
    
    futures = scheduler.compute(delayedCalls) # submit all delayed calculations; obtain futures immediately

    try:
        _distributed.progress(futures) # wait for all calculations to complete
    except Exception as e:
        print(E)
    print("\n")

    print("Waiting for all tasks to complete...")
    _distributed.wait(futures)

    results = {}
    errorsCount = 0
    newValuesCount = 0
    oldValuesCount = 0
    for f in futures:
        try:
            (taxId, CAI, isFreshValue) = scheduler.gather(f)
            results[taxId] = CAI
            if isFreshValue:
                newValuesCount += 1
            else:
                oldValuesCount += 1
            
        except Exception as e:
            print(e)
            errorsCount += 1
            
    print("Finished %d species with %d errors" % (len(results), errorsCount))
    print("{} new values; {} old values".format(newValuesCount, oldValuesCount))
    return results
            
        
print(runDistributed())

# DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #
#annotateCUBmeasures(511145)
#annotateCUBmeasures(1198115)
#annotateCUBmeasures(312017)
#annotateCUBmeasures(999415)
# DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #
