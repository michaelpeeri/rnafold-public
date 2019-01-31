import json
from random import sample
from collections import Counter
from intervaltree import IntervalTree, Interval
from gff import createGffDb


trueVals = frozenset(('true', 'True', 'TRUE', '1'))
class GenomeMoleculeModel(object):

    _gf_molecule = 1
    _gf_start    = 4
    _gf_end      = 5
    _gf_strand   = 7
    _gf_props    = 9
    
    def __init__(self, sequenceFile, gff, name, isLinear=None, props={}):
        self.sequenceFile = sequenceFile
        if type(gff)==type(""):
            self.gffFile = gff
            # TODO load file
        else:
            self.gff = gff
        
        self.name = name
        self.props = props
        
        if not isLinear is None: # get isLinear from props; allow isLinear to override
            self.isLinear = isLinear
        else:
            val = self.props.get('Is_circular', ['false'])
            self.isLinear = val[0] in trueVals

        print(self.props)
        self.features = IntervalTree()
        self._loadFeatures()

    def _loadFeatures(self):
        #print("mRNAs: {}".format( sum( 1 for _ in self.gff.region( seqid=self.name, featuretype='mRNA' ))))

        for feature in self.gff.region( seqid=self.name, featuretype='CDS' ):
            ft = feature.astuple()
            #print(ft)
            assert( ft[GenomeMoleculeModel._gf_molecule]==self.name )
            
            # ('transcript:ENST00000641515', '1', 'havana', 'mRNA', 65419, 71585, '.', '+', '.', '{"ID":["transcript:ENST00000641515"],"Parent":["gene:ENSG00000186092"],"Name":["OR4F5-202"],"biotype":["protein_coding"],"tag":["basic"],"transcript_id":["ENST00000641515"],"version":["2"]}', '[]', 4681)
            # ('transcript:ENST00000335137', '1', 'ensembl', 'mRNA', 69055, 70108, '.', '+', '.', '{"ID":["transcript:ENST00000335137"],"Parent":["gene:ENSG00000186092"],"Name":["OR4F5-201"],"biotype":["protein_coding"],"ccdsid":["CCDS30547.1"],"tag":["basic"],"transcript_id":["ENST00000335137"],"transcript_support_level":["NA (assigned to previous version 3)"],"version":["4"]}', '[]', 4681)
            # ('transcript:ENST00000426406', '1', 'ensembl_havana', 'mRNA', 450703, 451697, '.', '-', '.', '{"ID":["transcript:ENST00000426406"],"Parent":["gene:ENSG00000284733"],"Name":["OR4F29-201"],"biotype":["protein_coding"],"ccdsid":["CCDS72675.1"],"tag":["basic"],"transcript_id":["ENST00000426406"],"transcript_support_level":["NA (assigned to previous version 2)"],"version":["3"]}', '[]', 4684)
            # ('transcript:ENST00000332831', '1', 'ensembl_havana', 'mRNA', 685679, 686673, '.', '-', '.', '{"ID":["transcript:ENST00000332831"],"Parent":["gene:ENSG00000284662"],"Name":["OR4F16-201"],"biotype":["protein_coding"],"ccdsid":["CCDS41221.1"],"tag":["basic"],"transcript_id":["ENST00000332831"],"transcript_support_level":["NA (assigned to previous version 3)"],"version":["4"]}', '[]', 4686)

            sprops = ft[GenomeMoleculeModel._gf_props] 
            props = json.loads( sprops )
            strand = ft[GenomeMoleculeModel._gf_strand]
            assert( strand=='+' or strand=='-' )
            data = {'strand':ft[GenomeMoleculeModel._gf_strand], 'props':sprops}
            start = ft[GenomeMoleculeModel._gf_start]
            end = ft[GenomeMoleculeModel._gf_end]

            assert(end >= start)
            self.features.addi( start, end+1, data )  # end is inclusive in GFF3 but not inclusive in intervaltree

    def findFeatureByStart( self, start ):
        val = self.features[start:start+1]
        if val:
            return list(val)[0]
        else:
            return None
        
    def findNextFeature(self, start, debug=False):
        # TODO handle circular chromosomes
        width = 100
        while (True):

            candidates = sorted( self.features[start:start+width] )

            if debug:
                print("Trying with width={} ({}->{})...".format(width, start, start+width))
                print("candidates: {}".format(candidates))

            if candidates:
                return candidates[0]

            width *= 5
            if width > 2000000:
                return None
            
    def findPrevFeature(self, start, debug=False):
        # TODO handle circular chromosomes
        width = 100
        while (True):
            
            candidates = sorted( self.features[start-width:start] )

            if debug:
                print("Trying with width={} ({}<-{})...".format(width, start-width, start))
                print("candidates: {}".format(candidates))

            if candidates:
                return candidates[-1]

            width *= 5
            if width > 2000000:
                return None

    def find3PrimeFlankingRegion( self, feature, debug=False ):
        
        if feature.data['strand']=='+':
            print("+")
            next = self.findNextFeature( feature.end, debug=debug )
            if next is None:
                return None
            
            print("[[{}----->{}]]    [{}--->{}]".format( feature.begin, feature.end, next.begin, next.end ))
            ret = (feature.end, next.begin)
        elif feature.data['strand']=='-':
            print("-")
            prev = self.findPrevFeature( feature.begin-1, debug=debug )
            if prev is None:
                return None
            
            print("[{}<---{}]    [[{}<------{}]]".format( prev.begin, prev.end, feature.begin, feature.end ))
            ret = (prev.end, feature.begin)
        else:
            raise ValueError("Unknown strand '{}".format(feature.data['strand']))

        print(ret)
        if ( ret[1] <= ret[0] ):
            return None
        
        return ret
        

        
class GenomeModel(object):

    _gf_mol_name  = 1
    _gf_mol_props = 9

    def __init__(self, sequenceFile, gffFile, isLinear, variant):
        self.sequenceFile = sequenceFile
        self.gffFile = gffFile
        self.isLinear = isLinear
        self.variant = variant

        print("Loading gff3 file...")
        self.gff = createGffDb(gffFile, variant)
        print("Done.")

        chromosomeFeatures = list(self.gff.features_of_type( "chromosome" ))
        moleculeNames = [c.astuple()[GenomeModel._gf_mol_name] for c in chromosomeFeatures]
        moleculeProps = [json.loads(c.astuple()[GenomeModel._gf_mol_props]) for c in chromosomeFeatures]
        #self.molecules = [Genome(name,props) for (name,props) in zip(moleculeNames, moleculeProps)]
        #assert(self.molecules)
        #print( [json.loads(c.astuple()[9]) for c in self.gff.features_of_type( "chromosome" )] )

        self.moleculeModels = [GenomeMoleculeModel(sequenceFile, self.gff, name, props=props) for (name,props) in zip(moleculeNames, moleculeProps)]

        print([x.name for x in self.moleculeModels])
        print([len(x.features) for x in self.moleculeModels])

def displayInterval(tree, begin, end):
    print(tree[begin-10:end+11])

def testForwardIteration(mol):
    allFeatures = sorted( mol.features.items() )
    print(len(allFeatures))

    iteratedFeatures = []
    node = mol.findNextFeature(0)
    iteratedFeatures.append(node)
    
    while (True):
        node = mol.findNextFeature(node.end+1)
        
        if node is None:
            break
        
        iteratedFeatures.append(node)
        
    print(len(iteratedFeatures))

    missingFeatures = frozenset(allFeatures) - frozenset(iteratedFeatures)
    print(missingFeatures)
    for f in missingFeatures:
        print("--> missing:")
        displayInterval( mol.features, f.begin, f.end )
        

def testBackwardIteration(mol):
    allFeatures = list( reversed( sorted( mol.features.items() )))
    print(len(allFeatures))

    last = allFeatures[0]

    iteratedFeatures = []
    node = mol.findPrevFeature(last.end+100)
    iteratedFeatures.append(node)
    
    while (True):
        node = mol.findPrevFeature(node.begin)
        
        if node is None:
            break
        
        iteratedFeatures.append(node)
    print(len(iteratedFeatures))

    missingFeatures = frozenset(allFeatures) - frozenset(iteratedFeatures)
    print(missingFeatures)
    for f in missingFeatures:
        print("--> missing:")
        displayInterval( mol.features, f.begin, f.end )

def test3primeFlankingRegions(mol):

    stats = Counter()

    successSet = set()
    failSet = set()
    over500set = set()
    over10set = set()
    over40set = set()
    over100set = set()
    over500set = set()
    
    for feat in mol.features.items():
        val = mol.find3PrimeFlankingRegion( feat )
        if val is None:
            failSet.add(feat)
            continue
        
        (begin, end) = val
        assert( end > begin )
        stats.update( (end-begin,) )
        successSet.add(feat)
        
        if end-begin >= 10:
            over10set.add(feat)
            if end-begin >= 40:
                over40set.add(feat)
                if end-begin >= 100:
                    over100set.add(feat)
                    if end-begin >= 500:
                        over500set.add(feat)

    print("Found {} flanking regions".format(sum(stats)))
    print(len(successSet))
    print(stats.most_common(10))

    print("Failed (sample):")
    for f in sample( failSet, 10 ):
        print(f)

    print("Over500 (sample):")
    for f in sample( over500set, 10 ):
        print(f)

        
    print("Count (len >  10): {}".format( len( over10set  )))
    print("Count (len >  40): {}".format( len( over40set  )))
    print("Count (len > 100): {}".format( len( over100set )))
    print("Count (len > 500): {}".format( len( over500set )))
    

    # Test cases:
    # Overlap over start codon
    print( mol.findFeatureByStart(550439) )
    mol.find3PrimeFlankingRegion( mol.findFeatureByStart(550439), debug=True )
    print(mol.findNextFeature(550493))

    print("---"*10)
    print(mol.features[59000:60000])
    print("---"*10)
    

    f59687 = mol.findFeatureByStart(59687)
    mol.find3PrimeFlankingRegion( f59687, debug=True )
    assert( f59687 in over100set )
    assert( f59687 not in over500set )
    
        

def testAll():
    # gm1 = GenomeModel(
    #     sequenceFile='/tamir1/mich1/cellfold/data/Ensembl/Homo.sapiens/Homo_sapiens.GRCh38.dna_rm.toplevel.fa.gz',
    #     gffFile='/tamir1/mich1/cellfold/data/Ensembl/Homo.sapiens/Homo_sapiens.GRCh38.95.gff3.gz',
    #     isLinear=True,
    #     variant="Ensembl" )
                               
    gm2 = GenomeModel(
        sequenceFile='/tamir1/mich1/data/Ensembl/Ecoli/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_rm.toplevel.fa.gz',
        gffFile='/tamir1/mich1/termfold/data/Ensembl/Ecoli/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz',
        isLinear=True,
        variant="Ensembl" )

    print("Forward:")
    testForwardIteration(  gm2.moleculeModels[0] )
    print("Backward:")
    testBackwardIteration( gm2.moleculeModels[0] )

    test3primeFlankingRegions( gm2.moleculeModels[0] )
    
    return 0

    

if __name__=="__main__":
    import sys
    sys.exit( testAll() )
