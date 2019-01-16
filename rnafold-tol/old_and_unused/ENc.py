from __future__ import print_function
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import NucleotideAlphabet
from Bio import SeqIO
import Bio.Data
#import numpy as np
#import pandas as pd
#import matplotlib
#matplotlib.use("cairo")
#import matplotlib.pyplot as plt
#plt.style.use('ggplot') # Use the ggplot style


fastafile = sys.argv[1]

translationTable = int(sys.argv[2])
codonTable = Bio.Data.CodonTable.unambiguous_dna_by_id[translationTable]
codonTableDict = codonTable.forward_table
assert(len(codonTableDict) <= 64 and len(codonTableDict) >= 57)
reverseTable = {}
for u,v in codonTableDict.iteritems():
    if v in reverseTable:
        reverseTable[v].append(u)
    else:
        reverseTable[v] = [u]
reverseTable["start"] = codonTable.start_codons
reverseTable["stop"] = codonTable.stop_codons
assert(len(reverseTable)==20+2)
#print("Genetic code: %s" % reverseTable)

debugMode = False
if( len(sys.argv) > 3 ):
    if( sys.argv[3]=='-v' ):
        debugMode = True

# def calcFgroups():
#     Fgroups = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[]}
#     for aa,codons in reverseTable.iteritems():
#         if aa=="start" or aa=="stop":
#             continue

#         Fgroups[len(codons)].append(aa)

#     assert(sum([len(x) for x in Fgroups.values()])==20)
        
#     return Fgroups

def getAAs():
    for aa in reverseTable.keys():
        if aa=="start" or aa=="stop":
            continue
        yield aa

def calcENc(seq):
    FforGroups = {2:[], 3:[], 4:[], 5:[], 6:[]}

    #codonseq =
    codonCounts = {}

    for codon in getCodonSeq(seq[:-3]):
        if codon in codonCounts:
            codonCounts[codon] += 1
        else:
            codonCounts[codon] = 1
    assert(sum(codonCounts.values()) == (len(seq)/3)-1)

    totalCodonsUsed = len(codonCounts)

    #if( totalCodonsUsed < 40 or len(seq) < 300 ):
    #    return None

    # Iterate over all AAs belonging to groups F2 and up (for each aa in F1, Nc must be 1)
    for aa in [aa for aa in getAAs() if len(reverseTable[aa]) > 1]:
        Fgroup = len(reverseTable[aa])
        assert(Fgroup > 1)

        # since we belong to Fx, there are x possible codons
        ni = [0]*Fgroup
            
        for i, codon in enumerate(reverseTable[aa]):
            if codon in codonCounts:
                ni[i] += codonCounts[codon]
        n = float(sum(ni))
        
        if( debugMode ):
            print("AA: %s Counts: %s" % (aa, ", ".join(map(str, ni))))

        if n > 1:
            G = 0.0
            for i in range(Fgroup):
                pi = float(ni[i])/n
                assert(pi>=0.0 and pi<=1.0)
                #print(i, pi)
                G += pi**2

            #print("G=%.3g" % G)

            assert(G >= 0.0 and G <= 1.0)
            #H = 1.0 - G
            #print("G=%.2g H=%.2g" % (G, H))

            Fhat = ((n * G) - 1.0) / (n-1)  # Wright (1990) Eq. (1);  Nei, Tajima (1981) Eq. (7)
            #Fhat = G  # TODO add ref.
        elif n==1:
            Fhat = 1.0
        elif n==0:
            Fhat = None

        #if( Fhat!=0.0 ): 
        #    print("Nc=%.3g" % (1/Fhat))
        #else:
        #    print("Nc=Inf")

        if( Fhat==0.0 ):
            Fhat = G

        if( debugMode ):
            print(Fhat)


        assert(n==0 or Fhat <= 1.0 or Fhat is None)
            
        #print("%.4g" % Fhat)
        if( not Fhat is None ):
            FforGroups[Fgroup].append(Fhat)

    #print("Fhat: %s" % FforGroups)

    Ngroups = [0.0]*7
    Ngroups[1] = 2.0 # TODO - make this generic

    # Calculate the Nc for each codon degeneracy class (which is the cardinality of that class, divided by the mean F (homozygosity))
    cardinalities = {2:9, 3:1, 4:5, 5:0, 6:3}  # TODO - make this generic
    missing = []
    for F, group in FforGroups.iteritems():
        if( debugMode ):
            print(F, group)
            
        if not group:
            Ngroups[F] = 0.0
            #assert(cardinalities[F]==0)
            missing.append(F)
            continue

        meanF = float(sum(group)) / len(group)
        nominator = float(cardinalities[F])

        if( meanF == 0.0 ):
            #print("zero for: %d, %s" % (F, group))
            missing.append(F)
            continue

        Ngroups[F] = nominator / meanF
        
        #print(F, nominator, meanF)
        
    if( debugMode ):
        print(Ngroups)

    # Handle groups that are completely unrepresented. In most cases, this means refusing to give a score.
    for F in missing:
        if( F==3 ):
            Ngroups[F] = (Ngroups[F-1] + Ngroups[F+1]) * 0.5
        #elif( F==2 ):
        #    #Ngroups[F] = filter(lambda x:x>0.0, Ngroups)[0]
        #elif( F==6 ):
        #    #Ngroups[F] = filter(lambda x:x>0.0, Ngroups)[-1]
        else:
            #print("%s is missing" % F)
            if( cardinalities[F] > 0 ):
                return None


    Nc = sum(Ngroups)
    #print(Nc)
    #assert(Nc>=21 and Nc<=(64 - len(reverseTable['stop'])))
    #if( Nc > totalCodonsUsed ):
    #    print("Effective: %f  real: %d" % (Nc, totalCodonsUsed))
    #assert(Nc <= totalCodonsUsed)
    return min(61.0, Nc)


#Fgroups = calcFgroups()

def getCodonSeq(seq):
    for i in range(len(seq)/3):
        yield seq[i*3:(i*3+3)]

with open(sys.argv[1], "r") as f:
    for entry in SeqIO.parse(f, 'fasta'):

        seq = str(entry.seq).upper()

        if( len(seq) % 3 != 0 ):
            print("Warning: skipping sequence %s because it does not contain a whole number of codons..." % entry.id)
            continue

        #translation = seq.seq.translate(table=codonTable)

        ENc = calcENc(seq[:-3]) # exclude the stop codon

                

            
            #print(aa, Fgroup)

        padding = " "*(20-len(entry.id))

        if( not ENc is None ):
            print("%s %s Nc = %2.3f" % (entry.id, padding, ENc))
        else:
            print("%s %s Nc = None" % (entry.id, padding))
    
        
