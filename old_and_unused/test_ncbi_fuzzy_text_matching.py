from reference_trees import FuzzyTaxonMatcher



matcher = FuzzyTaxonMatcher()

def doTest(query, scope, verbose=False):
    print("Q: %s" % query)
    match = matcher.match(query, scope, verbose)

    if match is None:
        print("R: ---")
        return None
    else:
        (taxId, matchedName, similarity) = match
        print("R: %s (%.2g%%)" % (matchedName, similarity*100))
        return matchedName

def test(query, scope, expectedResult = None):
    print("=="*50)
    print("=="*50)
    actualResult = doTest(query, scope, verbose=False)
    if actualResult != expectedResult:
        print("TEST FAILED")
        actualResult = doTest(query, scope, verbose=True)
        print(actualResult)
        print(expectedResult)
        return False
    else:
        return True

test("Curtissbacteria GWA1 OP11 40 16 partial", "Candidatus Curtissbacteria", "Candidatus Curtissbacteria bacterium GW2011_GWA1_40_16")
test("Woesebacteria GWD2 OP11 40 19", "Candidatus Woesebacteria", "Candidatus Woesebacteria bacterium GW2011_GWD2_40_19")
test("Woesebacteria GWA1 OP11 41 7", "Candidatus Woesebacteria", "Candidatus Woesebacteria bacterium GW2011_GWA1_41_7")
test("Woesebacteria GWA1 OP11 44 23", "Candidatus Woesebacteria", "Candidatus Woesebacteria bacterium GW2011_GWA1_44_23")
test("Woesebacteria 31 71", "Candidatus Woesebacteria", "Candidatus Woesebacteria bacterium RIFOXYA1_FULL_31_71")
test("Woesebacteria 34 12", "Candidatus Woesebacteria", "Candidatus Woesebacteria bacterium RBG_16_34_12") 
test("Woesebacteria 34 1x", "Candidatus Woesebacteria")
test("Woesebacteria 3x 12", "Candidatus Woesebacteria")
test("Woesebacteria 33 12", "Candidatus Woesebacteria") # ambiguious match
test("Woesebacteria 35 12", "Candidatus Woesebacteria")
test("Woesebacteria 33 11", "Candidatus Woesebacteria", "Candidatus Woesebacteria bacterium RIFOXYD1_FULL_33_11")
test("Woesebacteria 34 4", "Candidatus Woesebacteria")
test("Woesebacteria 42 9", "Candidatus Woesebacteria")

test("Nomurabacteria 36 19","Candidatus Nomurabacteria") # ambiguious match
test("Nomurabacteria 39 10","Candidatus Nomurabacteria")
test("Nomurabacteria 39 11","Candidatus Nomurabacteria")
test("Nomurabacteria 38 10","Candidatus Nomurabacteria")
test("Nomurabacteria 31 12","Candidatus Nomurabacteria", "Candidatus Nomurabacteria bacterium CG1_02_31_12")
test("Nomurabacteria 31 13","Candidatus Nomurabacteria")
test("Nomurabacteria 30 12","Candidatus Nomurabacteria")




