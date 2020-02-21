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
from itertools import compress, chain, imap
from copy import copy
from math import factorial
import random
from Bio.Seq import Seq
import Bio.Data
from local_cache import LocalStringsCache

def splitCodons(seq):
    assert(len(seq)%3==0)
    return [seq[i*3:i*3+3] for i in range(len(seq)/3)]

"""
Convert a forward codon table to a backwards codon table (AA->codons)
Note: make_back_table from Bio.Data.CodonTable only returns the first codon
Note: the returned table does not include include entries for start/stop codons
"""
def make_full_back_table(forward_table):
    ret = {}
    for codon,aa in forward_table.items():
        codon = codon.lower()
        if aa is None:
            continue
        
        if aa in ret:
            ret[aa].append(codon)
        else:
            ret[aa] = [codon]
    return ret

def applyCodonsPermutation(mutableSeq, newCodons, positions):
    for newCodon, pos in zip(newCodons, positions):
        mutableSeq[pos] = newCodon

negate = lambda p: imap( lambda x: not x, p ) # return a sequence where each element is negated

class SynonymousCodonPermutingRandomization(object):
    """
    Args:
    skipIdenticalSequences - refuse returning sequences identical to the source (even though they are valid permutation). This means the sampling from the population of random perumtations is done without replacement; This can be thought of as either distorting the distribution or avoiding distrortion; The difference is only noticeble for sequences which have a very small number of permutations (compared to the number actually generated).
    """
    def __init__(self, geneticCode=1, skipIdenticalSequences=True):
        self._geneticCode = geneticCode
        self._code = Bio.Data.CodonTable.unambiguous_dna_by_id[geneticCode]
        self._back_table = make_full_back_table(self._code.forward_table)
        self._skipIdenticalSequences = skipIdenticalSequences

    def randomize(self, nucleotideSeq):
        seq = Seq(nucleotideSeq.lower())
        codonsSeq = splitCodons(nucleotideSeq.lower())
        origSeqTranslation = seq.translate(table=self._code)

        randomizedCodonsSeq = copy(codonsSeq)
        permutationsCount = 1
        for aa,aaCodons in self._back_table.items(): # Iterate over each AA and its codons
            aaCodonsSet = frozenset(aaCodons)
            isSynonymousPosition = map( lambda x: x in aaCodonsSet, codonsSeq)
            synonymousPositions = filter(lambda x:x != -1, [i if v else -1 for v,i in zip(isSynonymousPosition, range(len(isSynonymousPosition)))])

            if(sum([int(x) for x in isSynonymousPosition]) < 2): # nothing to permute
                continue

            pool = list(compress(codonsSeq, isSynonymousPosition))

            # Calculate the total number of permutations, which is the product
            # of the numbers of permutations for each AA.
            # The permutations for each AA are a multiset permutation:
            # Source: https://en.wikipedia.org/wiki/Permutation#Permutations_of_multisets
            permutationsCountDenom = 1
            for c in aaCodons:
                countOfC = sum( [1 if x==c else 0 for x in codonsSeq] )
                permutationsCountDenom *= factorial(countOfC)
            permutationsCount *= factorial(len(pool)) / permutationsCountDenom

            # Permute synonymous codons for the current aa
            random.shuffle(pool)
            applyCodonsPermutation( randomizedCodonsSeq, pool, synonymousPositions)

        # Final checks
        # refuse providing permutations if the number of possible permutations is very small
        if(permutationsCount < 50):
            print("Warning: sequence only has %d possible permutations" % permutationsCount)
            #raise Exception("Sequence only has %d possible permutations" % permutationsCount)

        # calculate nucleotide identity to the original sequence
        identityCount = sum([x==y for x,y in zip(''.join(randomizedCodonsSeq), seq)])
        assert(identityCount >= 0 and identityCount <= len(codonsSeq)*3)
        identity = float(identityCount)/len(seq)

        # refuse returning a sequence identical to the source (even though it is a valid permutation)
        #if( self._skipIdenticalSequences and identity > 0.9999 ):
        #    return self.randomize(nucleotideSeq)

        # warn if the permutation only effects a small number of positions
        # note: identities are normally over 67%
        if( identity > 0.9 ):
            print("Warning: shuffled sequence has %.1f%% identity to original sequence" % (identity*100))


        resultingSeq = ''.join(randomizedCodonsSeq)
        assert(Seq(resultingSeq).translate(table=self._code) == origSeqTranslation)  # translation was not maintained by randomization!
        return (permutationsCount, identity, resultingSeq)

    """
    Permute a sequence, allowing substitutions only in positions specified by a mask. Positions can include an arbitrary list of codons, specified by positions

    nucleotideSeq - nucleotide coding sequence (using the species genetic code). Length must be divisible by 3.
                    Start and end codons are not treated specially.
    codonMask - List of logical values of length nucleotideSeq/3 (not nucleotideSeq!)
                True  = randomize this codon
                False = keep this codon
    """
    def randomizeWithMask(self, nucleotideSeq, codonMask):
        assert( len(nucleotideSeq)%3 == 0 )   # length must be divisible by 3
        nucleotideMask = list(chain(*zip(codonMask,codonMask,codonMask))) # repeat each element of codonMask 3 times
        assert( len(nucleotideMask) == len(nucleotideSeq) )  

        originalMaskedNucleotides   = ''.join(compress( nucleotideSeq,        nucleotideMask  ))
        unmaskedNucleotides         = ''.join(compress( nucleotideSeq, negate(nucleotideMask) ))
        assert( len(originalMaskedNucleotides)+len(unmaskedNucleotides) == len(nucleotideSeq) )   # each nucleotide must be included in exactly one subset
        maskedFraction = float(sum(nucleotideMask))/len(nucleotideSeq)  # the fraction of masked nucleotides (will be used to calculate %identity)

        if len(originalMaskedNucleotides) == 0:   # Masked (randomized area) covers nothing - just return the original sequence
            return (1, 100.0, nucleotideSeq)

        origSeqTranslation = Seq( nucleotideSeq ).translate(table=self._code)
        
        (permutationsCount, identity, randomizedMaskedNucleotides) = self.randomize( originalMaskedNucleotides )

        idxMasked   = iter(range(len(originalMaskedNucleotides)))
        idxUnmasked = iter(range(len(unmaskedNucleotides)))

        resultingSeq = ''.join(imap( lambda isMasked: randomizedMaskedNucleotides[next(idxMasked)] if isMasked else unmaskedNucleotides[next(idxUnmasked)],
                            nucleotideMask ))
        assert( len(resultingSeq) == len(nucleotideSeq) )
        #for i in range(0, len(nucleotideSeq), 40):
        #    print("")
        #    print( ''.join(map(lambda x: '+' if x else ' ', nucleotideMask[i:i+40] ) ) )
        #    print( nucleotideSeq[i:i+40] )
        #    print(           ret[i:i+40] )

        assert( all( map( lambda x: True if x[2] else x[0]==x[1],  zip( resultingSeq, nucleotideSeq, nucleotideMask ) ) ) ) # all unmasked nucleotide must remain unchanged
        assert(Seq(resultingSeq).translate(table=self._code) == origSeqTranslation)  # translation was not maintained by randomization!
        identity = identity*maskedFraction + 1.0*(1-maskedFraction)
        
        return (permutationsCount, identity, resultingSeq)
        
    def verticalPermutation( self, cdss ):
        assert(all([x%3==0 for x in map(len, cdss)]))  # length of all CDSs must be divisible by 3
        longestLengthNt = max(map(len, cdss)) # find the length of the longest CDS

        cdsCodons = map( splitCodons, cdss )

        identities = []
        
        for codonPos in range(longestLengthNt/3):
            print("--- Doing codon {} ---".format(codonPos))

            # Collect all codons in position codonPos
            poolCodons = []
            for seq in cdsCodons:
                if codonPos < len(seq):
                    codon = seq[codonPos]
                    assert(len(codon)==3)
                    poolCodons.append( codon )
            assert(len(poolCodons) >= 1)

            # Permute all codons in this position
            pool = ''.join(poolCodons)
            totalPermutationsCountForSeq, identity, shuffledPool = self.randomize(pool)

            # Update the original sequences with the new permutation
            nextPoolCodonPos = 0
            for seqIdx, seq in enumerate(cdsCodons):
                if codonPos < len(seq):
                    newCodon = shuffledPool[nextPoolCodonPos:nextPoolCodonPos+3]
                    #if len(newCodon) != 3:
                    #    print(codonPos)
                    #    print(len(newCodon))
                    #    print(shuffledPool)
                    assert(len(newCodon)==3)
                    seq[codonPos] = newCodon
                    assert(cdsCodons[seqIdx][codonPos] == newCodon)
                    nextPoolCodonPos += 3

            identities.append( identity )

        ret = list( map( lambda x: ''.join(x), cdsCodons ) )

        # Check all translations are unaltered
        for (u,v) in zip(cdss, ret):
            assert( len(u) > 0 )
            assert( len(u) == len(v) )
            xlation1 = Seq( u ).translate(table=self._code)
            xlation2 = Seq( v ).translate(table=self._code)
            if( xlation1 != xlation2 ):
                print("Translation mismatch:")
                print(u)
                print(v)
                print(xlation1)
                print(xlation2)
            assert(xlation1==xlation2)

        return (ret, identities)
            


class VerticalRandomizationCache(object):
    _id_sets = "x//sets"
    _id_native_seqs_written = "x//native_ok"
    #_id_genetic_code = "x//genetic_code"

    
        
    def _storeNativeSeqs( self, nativeSeqsMap ):
        if self._cache.get_value( VerticalRandomizationCache._id_native_seqs_written ) == "1":  # native seqs already written; nothing to do
            return

        # write all seqs
        for protId, seq in nativeSeqsMap.items():
            entry_key = self._make_entry_key( protId, -1 )
            self._cache.insert_value(entry_key, seq)

        self._cache.insert_value( VerticalRandomizationCache._id_native_seqs_written, "1" )  # make native seqs already written flag

    def _make_entry_key(self, protId, shuffleId ):
        assert(type(shuffleId) == type(0))
        assert(shuffleId >= -1)
        return "s//{}//{}".format(shuffleId, protId)
    
        
    def _nativeSeqsSource( self ):
        if self._cache.get_value( VerticalRandomizationCache._id_native_seqs_written ) != "1":  # were native seqs already written?
            raise Exception("Native seqs not available!")

        nativePrefix = "s//-1//"
        prefixLen = len(nativePrefix)
        for key, val in self._cache.all_matching_values_source( nativePrefix ):
            assert( key[:prefixLen] == nativePrefix )
            yield (key[prefixLen:], val)  # yield clean protIds without the prefix
        
    
    def __init__(self, shuffleType, taxId, nativeSeqsMap, geneticCode, randomizer):
        self._cache = LocalStringsCache("codon_randomization_cache_type_{}_taxid_{}".format( shuffleType, taxId ) )
        self._availableShuffles = self._checkAvailableSets()
        assert(type(geneticCode)==type(0))
        assert(geneticCode > 0)
        self._geneticCode = geneticCode
        self._storeNativeSeqs( nativeSeqsMap )
        self._randomizer = randomizer


    def _checkAvailableSets(self):
        val = self._cache.get_value(VerticalRandomizationCache._id_sets)
        if val is None:
            self._cache.insert_value(VerticalRandomizationCache._id_sets, "")
            return set()
        elif val=="":
            return set()
        else:
            return set(map(int, val.split(",")))

    def _updateAvailableSets(self):
        self._cache.update_value( VerticalRandomizationCache._id_sets, ",".join( map( str, self._availableShuffles ) ) )

    #def _storeGeneticCode(self):
    #    val = self._cache.get_value(VerticalRandomizationCache._id_genetic_code)
    #    if val is None:
    #        self._cache.insert_value(VerticalRandomizationCache._id_genetic_code, str(self._geneticCode) )
    #    else:
    #        assert( val == str(self._geneticCode) )

    def _storeNewShuffleSet( self, shuffleId ):
        nativeSeqs = []
        protIds = []
        for protId, seq in self._nativeSeqsSource():
            protIds.append(protId)
            nativeSeqs.append(seq)
        
        shuffledSeqs, _ = self._randomizer( nativeSeqs )
        assert( len(shuffledSeqs) == len(nativeSeqs) )
        assert( map(len, shuffledSeqs) == map(len, nativeSeqs) )

        for protId, shuffledSeq in zip( protIds, shuffledSeqs ):
            entry_key = self._make_entry_key( protId, shuffleId )
            self._cache.insert_value(entry_key, shuffledSeq)

        # Add the newly created shuffle to the list of available shuffles
        self._availableShuffles.add( shuffleId )
        self._updateAvailableSets()


    def getShuffledSeq(self, protId, shuffleId ):
        entry_key = self._make_entry_key( protId, shuffleId )
        
        if shuffleId in self._availableShuffles or shuffleId==-1:  # is this shuffleId already stored?
            return self._cache.get_value( entry_key )

        self._storeNewShuffleSet( shuffleId )
        
        assert(shuffleId in self._availableShuffles)  # now shuffleId must be stored...
        return self._cache.get_value( entry_key )

        

#map( lambda x: x[1] if x[0] else x[2], zip( [0,1,1,0,1,0,0,1,0,1,0,1,1,1,0,1,0,1,0,1], [8]*20, [1]*20 ) )

def testCountRandomizations(N=50000):
    c = SynonymousCodonPermutingRandomization()
    s = set()
    testSeq  = 'atcccgcgcccacctaataacacagcgatcagaccgctcagaccacatatacgatcggactcg'
    #testSeq  = 'atcccgcgcccacctaataacacagcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcgatcgcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcgatcgcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcgatcgcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcg'
    totalPermutationsCount = None
    for i in range(N):
        totalPermutationsCount, identity, perm = c.randomize(testSeq)
        s.add( perm )
    print('-----------------')
    print('Iters: %d; Found %d unique randomizations' % (N, len(s)))
    print('Sequence length: %d codons' % (len(testSeq)/3))
    print('Total randomizations: %g' % float(totalPermutationsCount))
    return 0

def testMaskedRandomization(N=10000):
    c = SynonymousCodonPermutingRandomization()
    s = set()
    totalPermutationsCount = None
    #
    #testSeq =  'atcgcggcagcagcggcggcggcggcggca'
    #testMask = [  0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
    #testMask = [  0, 0, 0, 0, 0, 0,  1, 1, 1, 1]
    #
    testSeq  = 'atcccgcgcccacctaataacacagcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcgatcgcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcgatcgcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcgatcgcgatcagcgaaccacatatacgatcggaaagccctacgcgagagcactcg'
    #testMask = [  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    testMask = [  0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #
    #testSeq  =       'atcccgcgcccacctaataacacagcgatcagaccgctcagaccacatatacgatcggactcg'
    #testMask =       [  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    #
    #testSeq  =  'xxxxxxatcccgcgcccacctaataacacagcgatcagaagaagaagaagaccgctcagaccacatatacgatcggactcgcccccccccccc'
    #testMask =  [  0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
    #
    #testSeq  = 'atctaataaatcccgcgcccacctaataacacagcgatcagaagaagaagaagaccgctcagaccacatatacgatcggactcgtatacgatcgga'
    #testMask = [  0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]

    assert(len(testSeq)/3==len(testMask))
    for i in range(N):
        totalPermutationsCount, identity, perm = c.randomizeWithMask(testSeq, testMask)
        s.add( perm )
    print('-----------------')
    print('Iters: %d; Found %d randomizations' % (N, len(s)))
    print('Sequence length: %d codons' % (len(testSeq)/3))
    print('Total randomizations: %g' % float(totalPermutationsCount))
    return 0


def testVerticalShuffling():
    from nucleic_compress import randseq
    
    geneticCode = 1
    numSeqs = 505
    fixedSeqLength = 621
    
    protIds = ["PROT{}.1".format(x) for x in range(numSeqs) ]
    seqs = list(map( lambda _: randseq(fixedSeqLength), protIds ))


    scpr = SynonymousCodonPermutingRandomization( geneticCode )
    randomizer = lambda cdss: scpr.verticalPermutation( cdss )
    cache = VerticalRandomizationCache(shuffleType=9999,
                                       taxId=555,
                                       nativeSeqsMap=dict(zip(protIds, seqs)),
                                       geneticCode=geneticCode,
                                       randomizer=randomizer )

    print(">native")
    s = cache.getShuffledSeq( "PROT9.1", -1 )
    for x0 in range(0, 80, len(s)): print( s[x0:x+80] )
    for i in range(10):
        s = cache.getShuffledSeq( "PROT9.1", i )
        print(">s{}".format(i))
        for x0 in range(0, 80, len(s)): print( s[x0:x+80] )

    return 0


def printCounterAsHistogram(counter, stops=(0,1,5,20,50,100,200,300,400,500,600,700,800,900,1000,1500)):
    levels = [] # [0]*(len(stops)+2)
    ranges = []
    #currLevel = 0

    def getCounterElementByRange( s0, s1 ):
        return sum([x[1] for x in counter.items() if (x[0]>=s0 and x[0]<s1) ] )

    if counter:
        levels.append( getCounterElementByRange( min( min(counter.keys()), stops[0]-1), stops[0] ) )
    ranges.append( ("min", str(stops[0])) )
    #currLevel += 1

    for s0,s1 in zip(stops,stops[1:]):
        assert(s0<s1)
        levels.append( getCounterElementByRange( s0, s1 ) )
        ranges.append( (str(s0), str(s1) ) )
        #currLevel += 1

    if counter:
        levels.append( getCounterElementByRange( stops[-1], max( max(counter.keys()), stops[-1]+1) ) )
    ranges.append( (str(stops[-1]), "max" ) )
    #currLevel += 1

    for level, currRange in zip( levels, ranges ):
        if counter:
            print("{: >5}-{: <5} {:7} {}".format( currRange[0], currRange[1], level, "="*int(round(float(level)/sum(counter.values())*200)) ) )
        else:
            print("{: >5}-{: <5} {:7}".format(    currRange[0], currRange[1], level ) )
        


def testAll():
    #ret = testCountRandomizations(N=1000)
    #if ret: return ret

    #ret = testMaskedRandomization(N=20000)
    #if ret: return ret

    ret = testVerticalShuffling()
    if ret: return ret
    
    return 0


if __name__=="__main__":
    import sys
    sys.exit(testAll())
