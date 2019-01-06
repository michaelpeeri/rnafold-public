

def testGroupsAC_AG_AT(sequence):
    yield "ac"
    yield "ag"
    yield "at"

def genGroupMembers(sequence, group):
    for start in range(len(sequence)-1):
        s =sequence[start:start+2]
        if s == group:
            yield((start, s))
    
def generatePermutations(sequence, groups):
    while(True):
        perm = list(sequence)
        for group in groups(sequence):
            members = list(genGroupMembers(sequence, group))
            print(group, members)

            rev = group[1] + group[0]

            for m in members:
                perm[m[0]:m[0]+2] = rev.upper()
                
        yield ''.join(perm)

        
def tests_main():

    sequence = "agctagctagccgcgagagacgatcgatcgatgctagctagctagctagctagctaccggcgcgagaatatatatagcgcatgcgagagagctagt"
    
    c = 0
    for i in generatePermutations(sequence, testGroupsAC_AG_AT):
        print(i)
        c += 1
        if( c > 20 ):
            break
        
        
        

if __name__=="__main__":
    import sys
    sys.exit(tests_main())
