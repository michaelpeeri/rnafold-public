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

