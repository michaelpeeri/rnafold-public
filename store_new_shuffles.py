import mysql_rnafold as db
from data_helpers import CDSHelper, AddShuffledSequences
from codon_randomization import SynonymousCodonPermutingRandomization

    
def storeRandomizedSequences(cds, seqs, shuffleIds):
    assert(len(seqs)==len(shuffleIds))
    addseqs = AddShuffledSequences(cds.getTaxId(), cds.getProtId(), db.Sources.ShuffleCDSv2_python)

    newSequenceIds = []
    for seq, shuffleId in zip(seqs, shuffleIds):
        newSequenceId = addseqs.addSequence( shuffleId, seq )
        newSequenceIds.append(newSequenceId)

    return newSequenceIds



def createRandomizedSeqs(cds, newShuffleIds):
   
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
            raise Exception("Identity of randomized sequence is too high - %.3g%% (length=%d nt, total permutations=%.2g)" % (identity*100.0, len(newseq), totalPermutationsCount))
        
        if(totalPermutationsCount < 200):
            raise Exception("Low number of possible permutations %.2g (length=%d nt, identity=%.3g%%)" % (totalPermutationsCount, len(newseq), identity*100.0))
        newShuffles.append( newseq )
    
    return newShuffles


def storeNewShuffles(taxId, protId, newShuffleIds):
    cds = CDSHelper(taxId, protId)
    return storeRandomizedSequences(cds,
                                    createRandomizedSeqs(cds, newShuffleIds),
                                    newShuffleIds)
    
