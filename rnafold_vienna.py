import re
import subprocess

# ----------------------------------------------------------
# configuration
vienna_rnafold_path = "~/anaconda2_node/bin/RNAfold"
#vienna_rnafold_path = "~/anaconda2/bin/RNAfold"
# ----------------------------------------------------------


# Parse vienna 'rnafold' output
# Example:
#
#
#[line 1]:   ccagucgaccagacuauauacaaccuacgcguaucgcgcga
#[line 2]:   ..((((.....))))............(((((...))))). ( -9.50)
#
reMFEScore = re.compile(".*\n.*[(]\s*([\d.-]+)[)]\n")   # parse the score (group 1), ignore the structure
reMFEScoreWithStructure = re.compile(".*\n([().]+)\s+[(]\s*([\d.-]+)[)]\n")  # return the score (group 2) and structure (group 1)

def RNAfold_direct(seq):
    out = subprocess.check_output("echo %s | %s --noPS" % (seq, vienna_rnafold_path), shell=True)
    score = float(reMFEScore.match(out).group(1))
    assert(score<=0.0)
    return score


"""
Parse the 2nd-ary structure string representation returned by vienna,
e.g.: (((.....)))..........((.(((((...))))).)).'
"""
def parseStructure(struct):
    stack = []

    ret = []

    for pos, s in enumerate(struct):
        if s=='(':
            stack.append(pos)
        elif s==')':
            opposingPos = stack.pop()
            assert(not opposingPos is None)
            assert(opposingPos < pos)  # this will obviously always be true...
            #print("%d-%d" % (opposingPos, pos))
            ret.append( (opposingPos, pos) )
        elif s=='.':
            pass
        else:
            raise Exception("Structure representation contains unsupported character '%s'" % s)
    assert(len(stack)==0)  # after parsing, all pairs must have been matched...

    return sorted(ret)
    
    

def RNAfoldWithStructure(seq):
    out = subprocess.check_output("echo %s | %s --noPS" % (seq, vienna_rnafold_path), shell=True)
    match = reMFEScoreWithStructure.match(out)
    score = float(match.group(2))
    assert(score<=0.0)

    struct = match.group(1)
    parsedPairs = parseStructure(struct)

    return (parsedPairs, score)
    

def testsMain():
    print(RNAfoldWithStructure("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"))
    print(RNAfoldWithStructure("gucgaccagacuauauacaaccuacgcguaucgcgcgaagc"))
    print(RNAfoldWithStructure("cgcgcgcgcgcgcgcgcgcgcgcgcgcgcgcgcgcgcgcgc"))
    print(RNAfoldWithStructure("cgcgcgcgcgcgcgcgcgcgcatatatatatatatatatat"))
    print(RNAfoldWithStructure("cgcgcgcgctatatatatatatatatatatatagcgcgcgc"))
    print(RNAfoldWithStructure("atatatatacgcgcgcgcgcgcgcgcgcgcgcgtatatata"))
    return 0
    
if __name__=="__main__":
    import sys
    sys.exit(testsMain())
    
    
