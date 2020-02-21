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
from binascii import crc32
from Bio import SeqIO

fastafile1 = sys.argv[1]
fastafile2 = sys.argv[2]

def fastaSource(filename):
    for seq in SeqIO.parse(open(filename, 'r'), 'fasta'):
        yield seq

def getCrc(seq):
    return crc32(str(seq.seq).lower()) & 0xffffffff

hashes = {}
for seq in fastaSource(fastafile1):
    crc = getCrc(seq)
    hashes[seq.id] = crc


countIdentical = 0
countDifferent = 0
countNotFound = 0

for seq in fastaSource(fastafile2):
    crc2 = getCrc(seq)

    if( not seq.id in hashes):
        countNotFound += 1
        continue

    crc1 = hashes[seq.id]
    if( crc1==crc2 ):
        countIdentical += 1
    else:
        countDifferent += 1

print("Identical: %d  Different: %d  Not found: %d" % (countIdentical, countDifferent, countNotFound))
