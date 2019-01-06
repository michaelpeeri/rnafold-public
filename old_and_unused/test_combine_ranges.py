tests = (("1", (0,1,2,4,6,7,8,9,10,13,14,16,17,20,21,22,23,24,25,26,27,28,29,30,32)),
         ("2", (0,2,4,5,7,8,9,12,13,14,15,18,19,20,21,22,25,26,27,28,29,30,31,34,35)))

def makeRangesList(l):
    return [(x,1) for x in l]

def combineRanges(l):
        

for label, testData in tests:
    print(label)
    print(makeRangesList(testData))
    
