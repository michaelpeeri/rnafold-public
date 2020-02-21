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
#
from __future__ import print_function
import sys
from collections import Counter
from math import log10
from random import randint
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import NucleotideAlphabet
from Bio import SeqIO
import Bio.Data
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style

translationTable = int(sys.argv[2])
codonTable = Bio.Data.CodonTable.unambiguous_dna_by_id[translationTable] # Sorry, this uses undocumented API :-(
codonTableDict = codonTable.forward_table
print(codonTable.start_codons)
print(codonTable.stop_codons)
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
print(reverseTable)


def getCodons(seq):
    if(len(seq) % 3 != 0):
        print("Warning: sequence does not contain complete codons, skipping...")
        return
        
    for n in range(0, len(seq)-3, 3):  # Skip the stop codon
        codon = seq[n:(n+3)]
        if( 'N' in codon ):
            continue
        yield codon


# def enumPossibleCodingSequences(translated):
#     if len(translated)==0:
#         return

#     if len(translated)==1:
#         for codon in reverseTable[translated[0]]:
#             yield codon

#     for codon in reverseTable[translated[0]]:
#         for i in enumPossibleCodingSequences(translated[1:]):
#             yield codon + i

#def xcount(iterable):
#    x = 0
#    for _ in iterable:
#        x += 1
#    return x
            
# print("------------------------")
# s = "YKMML"
# print(s)
# print(xcount(enumPossibleCodingSequences(s)))
# s = "YKMMLAVVVILLA"
# print(s)
# print(xcount(enumPossibleCodingSequences(s)))
# s = "YKMMLAVVVILLAWVRYKMMLAVVVILLAWVRYKMMLAVVVILLAWVRYKMMLAVVVILLAWVR"
# print(s)
# print(xcount(enumPossibleCodingSequences(s)))



#def codonPMF(aas):

#    p0 = [1.0]

    



def frequenciesToProbabilities(freqs):
    total = float(sum(codonFrequency.values()))
    out = {}

    for u, v in freqs.iteritems():
        out[u] = float(v) / total
    return out

codonFrequency = Counter()
with open(sys.argv[1], "r") as f:
    for seq in SeqIO.parse(f, 'fasta'):
        codonFrequency.update(getCodons(str(seq.seq)))
assert(len(codonFrequency) <= 64)

Pcodon = frequenciesToProbabilities(codonFrequency)
assert(abs( sum(Pcodon.values()) - 1.0) < 1e-6)
print("Pcodon: %s" % str(Pcodon))

t = {}
sigmaP = 0.0
for aa, codons in reverseTable.iteritems():
    if(aa == "start" or aa == "stop"):
        continue
    totalP = sum([Pcodon[codon] for codon in codons])

    gc = 0.0
    sp = 0.0
    for codon in codons:
        p = Pcodon[codon] / totalP
        sp += p
        #gc += float(codon.count("G") + codon.count("C")) / 3.0 * p
        gc += float(codon.count("G") + codon.count("C")) / 3.0 * p
    #overall = gc / len(codons)
    overall = gc
    assert(abs(sp-1.0) < 1e-6)
        
    sigmaP += totalP
    #print(aa, totalP)
    t[aa] =overall
assert(len(t)==20)
assert(abs(sigmaP-1.0) < 1e-6)
print(t)











def combine(a, b):
    out = [0.0] * (len(a)+len(b))

    for n1,p1 in enumerate(a):
        for n2,p2 in enumerate(b):
            out[n1+n2] += p1*p2
    return out

def codonPMF(aa):
    codons = reverseTable[aa]
    ret = [0]*4 # = [ E(g+c=0), E(g+c=1), E(g+c=2), E(g+c=3) ]
    total = sum( [Pcodon[codon] for codon in codons] )
    for codon in codons:
        gccount = codon.count("G") + codon.count("C")
        weight = Pcodon[codon] / total
        #weight = 1./len(codons)
        ret[gccount] += weight
    return ret
    
AA_PMF = {}
for a in reverseTable.keys():
    if( a != "stop" and a != "start" ):
        AA_PMF[a] = codonPMF(a)
print(AA_PMF)

def getProteinPMF(aas, verbose=False):
    pmf = [1.0]

    if verbose:
        print("-"*20)
        print("translation: %s" % aas)
        print("AA,E0,E1,E2,E3")
    
    for aa in aas:
        pmf = combine(pmf, AA_PMF[aa])
        if verbose:
            print("%s,%s" % (aa, ",".join(map(lambda x:"%.3g"%x, AA_PMF[aa]))))

    if verbose:
        for n,p in enumerate(pmf):
            print("%d,%.5g" % (n, p))

    return pmf

def test(aas):
    print("-"*20)
    print(aas)
    pmf = getProteinPMF(aas, True)
    print(pmf)
    assert(abs(sum(pmf) - 1.0) < 1e-8)
    
#test("W")
#test("WW")
#test("YKMML")
test("YKMMLAVVVILLA")
#test("YKMMLAVVVILLAWVRYKMMLAVVVILLAWVRYKMMLAVVVILLAWVRYKMMLAVVVILLAWVR")
#test("AAAAAAAAAA")
#test("TTTTTTTTTT")
#test("YYYYYYYYYY")


def getPercentile( value, pmf ):
    pmf = np.array(pmf)
    assert(abs(sum(pmf)-1.0) < 1e-6)

    cmf = np.cumsum(pmf)
    if value<len(pmf):
        psmaller = cmf[value]
    else:
        psmaller = 1.0

    if( value > 0 ):
        plarger = 1.0 - cmf[value-1]
    else:
        plarger = 1.0

    if value<len(pmf) and value>=0:
        pmfval = pmf[value]
    else:
        pmfval = 0.0
        
    #print(value, psmaller, pmfval, plarger)
    assert(abs(psmaller + plarger - pmfval - 1.0) < 1e-6 )

    if psmaller > plarger:
        direction = -1
    elif plarger > psmaller:
        direction = 1
    else:
        direction = 0
    #return (min(psmaller, plarger), direction)
    if( min(psmaller, plarger) > 1e-200 ):
        return log10(min(psmaller, plarger)) * direction
    else:
        return log10(1e-200) * direction


print("************************")
print(getPercentile( 0,  [0.1]*10 ))
print(getPercentile( 1,  [0.1]*10 ))
print(getPercentile( 2,  [0.1]*10 ))
print(getPercentile( 3,  [0.1]*10 ))
print(getPercentile( 4,  [0.1]*10 ))
print(getPercentile( 5,  [0.1]*10 ))
print(getPercentile( 6,  [0.1]*10 ))
print(getPercentile( 7,  [0.1]*10 ))
print(getPercentile( 8,  [0.1]*10 ))
print(getPercentile( 9,  [0.1]*10 ))
print(getPercentile( 10,  [0.1]*10 ))




# -------------------------------------------------------------------
# Calculate mean neutral GC% "tendency" for each CDS
# (i.e., the expected GC% based on the translation and the genomic codon composition)
meanNeutralGC = []
gccontent = []
tendencypval = []
with open(sys.argv[1], "r") as f:
    for seq in SeqIO.parse(f, 'fasta'):
        logp = 0.0

        if( len(seq) % 3 != 0 ):
            print("Warning: skipping sequence %s because it does not contain a whole number of codons..." % seq.id)
            continue

        translation = seq.seq.translate(table=codonTable)


        # Calculate the actual GC%
        s = seq.seq.upper()
        gcCount = s.count("C") + s.count("G")
        atCount = s.count("A") + s.count("T")

        gc = float(gcCount) / (gcCount+atCount)
        gccontent.append(gc)
        

        # Calculate the mean tendency
        tend = sum( [ t[aa] for aa in translation[:-1] ] ) / (len(translation)-1)
        meanNeutralGC.append(tend)

        gc_pmf = getProteinPMF(translation[:-1], False) # randint(1,50)==1)

        pc = getPercentile( gcCount, gc_pmf )
        tendencypval.append(pc)

        print(gc, tend, pc)

        
        for codon in getCodons(str(seq.seq)):
            logp += log10(Pcodon[codon])
        #print(logp)

        # DEBUG ONLY - REMOVE THIS - DEBUG ONLY - REMOVE THIS - DEBUG ONLY - REMOVE THIS - DEBUG ONLY - REMOVE THIS
        if( len(gccontent) > 500 ):
            break
        # DEBUG ONLY - REMOVE THIS - DEBUG ONLY - REMOVE THIS - DEBUG ONLY - REMOVE THIS - DEBUG ONLY - REMOVE THIS


fig, ax = plt.subplots()
plt.scatter( gccontent, meanNeutralGC, s=10)
l1 = np.arange(min(gccontent), max(gccontent), 0.01)
plt.plot( l1, l1, c="red")
plt.xlabel('Actual CDS GC%')
plt.ylabel('Mean CDS neutral GC% tendency')
plt.savefig('%s.gc-vs-tendency.pdf' % sys.argv[1])
plt.savefig('%s.gc-vs-tendency.svg' % sys.argv[1])



ttt = meanNeutralGC[:]
ttt.sort()
print(ttt[:10])
print(ttt[-10:])
df = pd.DataFrame(ttt)
        
fig, ax = plt.subplots()

histBins = np.arange(0.14, 0.85+0.01, 0.01)
hist, _ = np.histogram( ttt, bins = histBins )

plt.bar( np.arange(0.14, 0.85, 0.01), hist, width=0.01)

plt.annotate( s='n=%d' % len(ttt), xy=(0.6, 130) )
plt.annotate( s='low 25%% <= %.3g' % df.quantile(0.25), xy=(0.6, 115) )
plt.annotate( s='mean=%.3g' % np.array(ttt).mean(), xy=(0.6, 100) )
plt.annotate( s='median=%.3g' % np.median(np.array(ttt)), xy=(0.6, 85) )
plt.annotate( s='high 25%% >= %.3g' % df.quantile(0.75), xy=(0.6, 70) )
#ax.legend(fontsize='small', loc='upper left')

#plt.grid()
plt.title('Mean neutral GC% tendency')
plt.xlabel('mean GC%')
plt.ylabel('# CDSs')

plt.savefig('%s.gc-tendency.pdf' % sys.argv[1])
plt.savefig('%s.gc-tendency.svg' % sys.argv[1])



fig, ax = plt.subplots()
plt.scatter( gccontent, tendencypval , s=10)

plt.ylim([max(min(tendencypval), -50), min(max(tendencypval), 50)])

#l1 = np.arange(min(gccontent), max(gccontent), 0.01)
#plt.plot( l1, l1, c="red")
plt.xlabel('Actual CDS GC%')
plt.ylabel('Signed log(p-value) for neutral model')
plt.savefig('%s.gc-vs-tendency-pval.pdf' % sys.argv[1])
plt.savefig('%s.gc-vs-tendency-pval.svg' % sys.argv[1])






fig, ax = plt.subplots()

hist, histBins = np.histogram( tendencypval, bins = 100, normed=True )

plt.bar( histBins[:-1], hist, width=float(histBins[-1]-histBins[0])/100 )

#plt.annotate( s='n=%d' % len(ttt), xy=(0.6, 130) )
#plt.annotate( s='low 25%% <= %.3g' % df.quantile(0.25), xy=(0.6, 115) )
#plt.annotate( s='mean=%.3g' % np.array(ttt).mean(), xy=(0.6, 100) )
#plt.annotate( s='median=%.3g' % np.median(np.array(ttt)), xy=(0.6, 85) )
#plt.annotate( s='high 25%% >= %.3g' % df.quantile(0.75), xy=(0.6, 70) )
#ax.legend(fontsize='small', loc='upper left')

#plt.grid()
plt.title('Neutral model log(p-value)')
plt.xlabel('Neutral model log(p-value), signed')
plt.ylabel('Frequency')

plt.savefig('%s.p-val-hist.pdf' % sys.argv[1])
plt.savefig('%s.p-val-hist.svg' % sys.argv[1])
