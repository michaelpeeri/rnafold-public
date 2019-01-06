import argparse
from hashlib import md5
from collections import Counter
from rate_limit import RateLimit
from data_helpers import CDSHelper, SpeciesCDSSource
import mysql_rnafold as db
from codon_randomization import splitCodons, printCounterAsHistogram


rl = RateLimit(15)


def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert


argsParser = argparse.ArgumentParser()
argsParser.add_argument("--species", type=parseList(int), required=True)
argsParser.add_argument("--shuffle-type", type=str, default="11")
args = argsParser.parse_args()

species = args.species
shuffleType = args.shuffle_type
#defaultMappingType = (db.Sources.ShuffleCDSv2_matlab, db.Sources.ShuffleCDSv2_python)
shuffleTypesMapping = {""                   :db.Sources.ShuffleCDSv2_python,
                       "11"                 :db.Sources.ShuffleCDSv2_python,
                       "ShuffleCDSv2_matlab":db.Sources.ShuffleCDSv2_python,
                       "ShuffleCDSv2_python":db.Sources.ShuffleCDSv2_python,
                       "12"                 :db.Sources.ShuffleCDS_vertical_permutation_1nt,
                       "ShuffleCDS_vertical_permutation_1nt"
                                            :db.Sources.ShuffleCDS_vertical_permutation_1nt }
shuffleType=shuffleTypesMapping[args.shuffle_type]
maxShuffles=20
maxCodons=2000

def getNthCodon(n, seq):
    start = n*3
    end = start+3
    if len(seq) < end:
        return None
    return seq[start:end]

warnings = Counter()

numUniqueShuffles = Counter()

for taxId in species:
    proteinsDone = 0
    
    #nativeColumns = [[] for x in range(maxCodons)]
    #shuffledColumns = [[[] for x in range(maxCodons)] for y in range(maxShuffles)]
    allNativeSeqs = {}
    allShuffledSeqs = {}
    
    for protId in SpeciesCDSSource(taxId):
        cds = CDSHelper(taxId, protId)
        warnings.update(("total-cds",))

        allIds = cds.shuffledSeqIds(shuffleType=shuffleType)[:maxShuffles]

        nativeSeq = cds.sequence()
        if( len(nativeSeq)%3 != 0 ):
            warnings.update(("has-broken-codons",))
            continue
        
        nativeCodons = Counter( splitCodons(nativeSeq) )

        hasMismatchedCodons = False
        allNativeSeqs[protId] = nativeSeq
        hashesForShuffles = set()
        

        #for i, c in enumerate(splitCodons(nativeSeq)[:maxCodons]):
        #    nativeColumns[i].append(c)

        shuffledSeqs = []

        for shuffleId in range(len(allIds)):
            shuffledSeq = cds.getShuffledSeq(shuffleId, shuffleType)
            shuffledCodons = Counter( splitCodons(shuffledSeq) )

            hashesForShuffles.add( md5(shuffledSeq).hexdigest() )


            if shuffledCodons != nativeCodons:
                warnings.update(("num-horizontal-codon-mismatch",))
                hasMismatchedCodons = True

            shuffledSeqs.append(shuffledSeq)

            #for i, c in enumerate(splitCodons(shuffledSeq)[:maxCodons]):
            #    shuffledColumns[shuffleId][i].append(c)

        numUniqueShuffles.update( (len(hashesForShuffles),) )
        if len(hashesForShuffles) != len(allIds):
            warnings.update(("has-unexpected-num-shuffles",))
            
        #print(nativeColumns)
        #print(shuffledColumns)
        allShuffledSeqs[protId] = shuffledSeqs
                

        if hasMismatchedCodons:
            warnings.update(("has-horizontal-codon-mismatch",))
            
        proteinsDone += 1

        # DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #
        #if proteinsDone > 0:
        #    break
        # DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #### DEBUG ONLY #

        if(rl()):
            print("{} done -- {}".format(proteinsDone, warnings) )
            
    print(len(allNativeSeqs))
    print(len(allShuffledSeqs))

    codon = 0
    while (True):  # Iterate over each codon, until no more sequences get to that position
        foundAnyCodons = False
        nativeCodons = []
        shuffledCodons = []
        #print("codon: {}".format(codon))

        for protId, nativeSeq in allNativeSeqs.items():
            nativeNthCodon = getNthCodon(codon, nativeSeq)
            
            if nativeNthCodon is None:
                continue    # this sequence doesn't reach this codon position

            nativeCodons.append(nativeNthCodon)
            foundAnyCodons = True


            # iterate over all shuffled seqs to collect codons at this position (for all shuffles)
            for shuffleId, shuffledSeq in enumerate( allShuffledSeqs[protId] ):
                shuffledNthCodon = getNthCodon(codon, shuffledSeq)
                assert(not shuffledNthCodon is None)

                if shuffleId >= len(shuffledCodons):
                    shuffledCodons.append([])
                shuffledCodons[shuffleId].append( shuffledNthCodon )

        hasVerticalCodonMismatches = False
        
        # iterate over all found shuffle-ids to compare frequencies
        nativeFreqs = Counter(nativeCodons)
        for shuffleId, shuffledCodons in enumerate( shuffledCodons ):
            shuffledFreqs = Counter(shuffledCodons)
            if nativeFreqs != shuffledFreqs:
                warnings.update(("num-vertical-codon-mismatch",))
                hasVerticalCodonMismatches = True

        if hasVerticalCodonMismatches:
            warnings.update(("has-vertical-codon-mismatch",))
                

        if not foundAnyCodons: break
        codon += 1
                

            
        

print("{} done -- {}".format(proteinsDone, warnings) )

printCounterAsHistogram( numUniqueShuffles, stops=tuple(range(1,21)) )
