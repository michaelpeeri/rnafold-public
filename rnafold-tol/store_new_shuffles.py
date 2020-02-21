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
import mysql_rnafold as db
from data_helpers import CDSHelper, AddShuffledSequences, SpeciesCDSSource, getSpeciesTranslationTable
from codon_randomization import SynonymousCodonPermutingRandomization, VerticalRandomizationCache
import logging

def storeRandomizedSequences(cds, seqs, shuffleIds, shuffleType):
    assert( len(seqs)==len(shuffleIds) )
    addseqs = AddShuffledSequences( cds.getTaxId(), cds.getProtId() )

    newSequenceIds = []
    for seq, shuffleId in zip(seqs, shuffleIds):
        newSequenceId = addseqs.addSequence( shuffleId, seq, shuffleType )
        newSequenceIds.append(newSequenceId)

    return newSequenceIds


# Return randomized sequences for the specified native sequence object
def createRandomizedSeqs(cds, newShuffleIds, shuffleType=db.Sources.ShuffleCDSv2_python):
   
    shuffler = SynonymousCodonPermutingRandomization(cds.getTranslationTable())

    nativeSeq = cds.sequence()
    #print(nativeSeq[:10])

    newShuffles = []
    for shuffleId in newShuffleIds:
        totalPermutationsCount, identity, newseq = None, None, None

        try:
            totalPermutationsCount, identity, newseq = shuffler.randomize(nativeSeq)
        except Exception as e:
            print(e)
            raise
        
        assert((identity <= 1.0) and (identity > 0.0))
        
        if(identity > 0.95):
            print("Warning: Identity of randomized sequence is high - %.3g%% (length=%d nt, total permutations=%.2g)" % (identity*100.0, len(newseq), totalPermutationsCount))
        
        if(totalPermutationsCount < 500):
            raise Exception("Low number of possible permutations %.2g (length=%d nt, identity=%.3g%%)" % (totalPermutationsCount, len(newseq), identity*100.0))
        newShuffles.append( newseq )
    
    return newShuffles


_caches = {}
    
def getRandomizedSequenceCacheForVerticalPermutations(taxId):
    global _caches

    if (taxId, db.Sources.ShuffleCDS_vertical_permutation_1nt) in _caches:
        cache = _caches[(taxId, db.Sources.ShuffleCDS_vertical_permutation_1nt)]
        
    else:
        # read all native sequences
        protIds = []
        cdss = []
        for protId in SpeciesCDSSource(taxId):
            cds = CDSHelper(taxId, protId)
            
            if( cds.length()%3 != 0 ):
                continue
            
            seq = cds.sequence()
            
            protIds.append(protId)
            cdss.append(seq)
            
        geneticCode = getSpeciesTranslationTable( taxId )
        scpr = SynonymousCodonPermutingRandomization( geneticCode ) 
        randomizer = lambda cdss: scpr.verticalPermutation( cdss )
        cache = VerticalRandomizationCache(shuffleType=db.Sources.ShuffleCDS_vertical_permutation_1nt,
                                           taxId=taxId,
                                           nativeSeqsMap=dict(zip(protIds, cdss)),
                                           geneticCode=geneticCode,
                                           randomizer=randomizer )
        _caches[(taxId, db.Sources.ShuffleCDS_vertical_permutation_1nt)] = cache
        print(_caches.keys())

        
    return cache


# Generate and store new randomized sequences for the specified native protein
def storeNewShuffles(taxId, protId, newShuffleIds, shuffleType=db.Sources.ShuffleCDSv2_python, dontStore=False):
    
    cds = CDSHelper(taxId, protId)
    print(protId)
    
    if shuffleType == db.Sources.ShuffleCDSv2_python:
        return storeRandomizedSequences(cds,
                                        createRandomizedSeqs(cds, newShuffleIds, shuffleType),
                                        newShuffleIds,
                                        shuffleType
        )
    
    elif shuffleType == db.Sources.ShuffleCDS_vertical_permutation_1nt:
        cache = getRandomizedSequenceCacheForVerticalPermutations( taxId )

        seqs = map( lambda shuffleId: cache.getShuffledSeq( protId, shuffleId ), newShuffleIds )
        print(seqs)

        if dontStore: return seqs
        
        return storeRandomizedSequences(cds,
                                        seqs,
                                        newShuffleIds,
                                        shuffleType)
                                
    else:
        raise Exception("Unsupported shuffleType={}".format(shuffleType))
    


if __name__=="__main__":
    import sys
    taxId = int(sys.argv[1])
    storeNewShuffles( taxId, "dummy_protid", range(20), db.Sources.ShuffleCDS_vertical_permutation_1nt, dontStore=True ) # force the creation of 20 sets of randomized sequences for the specified species




    
