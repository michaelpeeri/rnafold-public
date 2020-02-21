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
#import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from CAI import CAI

# def processFiles( args ):
#     a = [ 'CAI' ]
#     a.extend( ['-s', args.sequence  ] )
#     a.extend( ['-r', args.reference ] )
#     if not args.genetic_code is None:
#         a.extend( ['-g', str(args.genetic_code) ] )

#     out = subprocess.check_output(a, shell=False)

#     print(out)
#     return 0

def getCodingSequences(fn):
    ret = {}
    with open(fn, 'r') as f:
        for seq in SeqIO.parse( f, 'fasta' ):
            ret[seq.id] = seq.seq
    return ret
            

def processFiles( args ):
    seqs =       getCodingSequences( args.sequence )
    refs = list( getCodingSequences( args.reference ).values())
    

    for geneId, seq in seqs.items():
        if not args.genetic_code is None:
            cai = CAI( seq, reference=refs, genetic_code=args.genetic_code )
        else:
            cai = CAI( seq, reference=refs)
        print( "{}\t{}".format( geneId, cai ))
    

if __name__=="__main__":
    import sys
    import argparse
    
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument( "-s", "--sequence",     type=str )
    argsParser.add_argument( "-r", "--reference",    type=str )
    argsParser.add_argument( "-g", "--genetic-code", type=int, default=None )
    args = argsParser.parse_args()
    
    sys.exit( processFiles(args) )

