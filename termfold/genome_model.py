import json
from random import sample
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pyfaidx
from intervaltree import IntervalTree, Interval
from data_helpers import getSpeciesGenomeSequenceFile, getSpeciesGenomeAnnotationsFile, getSpeciesGenomeAnnotationsVariant
from gff import createGffDb


trueVals = frozenset(('true', 'True', 'TRUE', '1'))
class GenomeMoleculeModel(object):

    _gf_molecule = 1
    _gf_start    = 4
    _gf_end      = 5
    _gf_strand   = 7
    _gf_props    = 9
    
    def __init__(self, sequenceFile : str, gff, name : str, isLinear =None, props : dict ={}):
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

    def findFeatureByStart( self, start:int ):
        val = self.features[start:start+1]
        if val:
            return list(val)[0]
        else:
            return None
        
    def findNextFeature(self, start:int, debug:bool=False):
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
            
    def findPrevFeature(self, start:int, debug:bool=False):
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

    def find3PrimeFlankingRegion( self, feature, debug:bool =False ):
        
        if feature.data['strand']=='+':
            next = self.findNextFeature( feature.end, debug=debug )
            if next is None:
                return None

            if debug:
                print("+")
                print("[[{}----->{}]]    [{}--->{}]".format( feature.begin, feature.end, next.begin, next.end ))
                
            ret = (self.name, feature.end, next.begin)
            
        elif feature.data['strand']=='-':
            prev = self.findPrevFeature( feature.begin-1, debug=debug )
            if prev is None:
                return None

            if debug:
                print("-")
                print("[{}<---{}]    [[{}<------{}]]".format( prev.begin, prev.end, feature.begin, feature.end ))
                
            ret = (self.name, prev.end, feature.begin)
            
        else:
            raise ValueError("Unknown strand '{}".format(feature.data['strand']))

        if debug:
            print(ret)
            
        if ( ret[2] <= ret[1] ):
            return None
        
        return ret
        

class GenomeModel(object):

    _gf_mol_name  = 1
    _gf_mol_props = 9

    def __init__(self, sequenceFile:str, gffFile:str, isLinear =None, variant:str ="Ensembl", geneticCode:int =1):
        self.sequenceFile = sequenceFile
        self.gffFile = gffFile
        self.variant = variant
        self.geneticCode = geneticCode

        print("Loading gff3 file...")
        self.gff = createGffDb(gffFile, variant)
        print("Done.")

        chromosomeFeatures = list(self.gff.features_of_type( "chromosome" ))
        moleculeNames = [c.astuple()[GenomeModel._gf_mol_name] for c in chromosomeFeatures]
        moleculeProps = [json.loads(c.astuple()[GenomeModel._gf_mol_props]) for c in chromosomeFeatures]
        #self.molecules = [Genome(name,props) for (name,props) in zip(moleculeNames, moleculeProps)]
        #assert(self.molecules)
        #print( [json.loads(c.astuple()[9]) for c in self.gff.features_of_type( "chromosome" )] )

        self.moleculeModels = [GenomeMoleculeModel(sequenceFile, self.gff, name, isLinear=isLinear, props=props) for (name,props) in zip(moleculeNames, moleculeProps)]

        print([x.name for x in self.moleculeModels])
        print([len(x.features) for x in self.moleculeModels])

        self.fasta = pyfaidx.Fasta(sequenceFile)

    def getRegion(self, region:tuple, rc=False):
        (chromosome, begin, end) = region
        return self.fasta.get_seq(chromosome, begin, end, rc=rc)


class GenomeModelsCache(object):
    def __init__(self):
        self._cache = {}

    def __getitem__(self, taxId):
        val = self._cache.get(taxId, None)

        if not val is None:
            return val
        else:
            return self._init_item(taxId)
        
    def _init_item(self, taxId):
        genomeSeqFile      = getSpeciesGenomeSequenceFile( taxId )
        genomeAnnotFile    = getSpeciesGenomeAnnotationsFile( taxId )
        genomeAnnotVariant = getSpeciesGenomeAnnotationsVariant( taxId )
        geneticCode        = getSpeciesTranslationTable( taxId )
        if genomeSeqFile is None or genomeAnnotFile is None or geneticCode is None:
            raise ValueError("No supporting annotations for taxId={}".format(taxId))

        gm = GenomeModel(
            sequenceFile=genomeSeqFile,
            gffFile=genomeAnnotFile,
            isLinear=False,
            variant=genomeAnnotVariant,
            geneticCode=geneticCode )

        return gm

_genome_models_cache = GenomeModelsCache()

def getGenomeModelFromCache( taxId ):
    return _genome_models_cache[taxId]

    
def displayInterval(tree, begin:int, end:int):
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
    from collections import Counter

    stats = Counter()

    successSet = set()
    failSet = set()
    over500set = set()
    over10set = set()
    over40set = set()
    over100set = set()
    over250set = set()
    over500set = set()
    
    for feat in mol.features.items():
        val = mol.find3PrimeFlankingRegion( feat )
        if val is None:
            failSet.add(feat)
            continue
        
        (chromosome, begin, end) = val
        assert( end > begin )
        stats.update( (end-begin,) )
        successSet.add(feat)
        
        if end-begin >= 10:
            over10set.add(feat)
            if end-begin >= 40:
                over40set.add(feat)
                if end-begin >= 100:
                    over100set.add(feat)
                    if end-begin >= 250:
                        over250set.add(feat)
                        if end-begin >= 500:
                            over500set.add(feat)

    print("Found {} flanking regions".format(sum(stats)))
    print(len(successSet))
    print(stats.most_common(10))

    print("Failed (sample):")
    for f in sample( failSet, 10 ):
        print(f)

    print("Over100 (sample):")
    for f in sample( over100set, 20 ):
        print(f)

        
    print("Count (len >  10): {}".format( len( over10set  )))
    print("Count (len >  40): {}".format( len( over40set  )))
    print("Count (len > 100): {}".format( len( over100set )))
    print("Count (len > 250): {}".format( len( over250set )))
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
    #assert( f59687 not in over500set )  - fails because some CDSs in the middle of this UTR are missing

def CDS3PrimeFlankingRegionSource( mol, minLength:int=10, debug=False ):
    for feat in mol.features.items():
        region = mol.find3PrimeFlankingRegion( feat, debug=debug )
        if not region is None:
            length = region[2]-region[1]
            if length >= minLength:
                yield( feat, region )
    

def CDSWith3PrimeSequencesSource( genomeModel, minLength:int=10, debug=False ):
    numRegions = 0
    
    for mol in genomeModel.moleculeModels:
        for feat, region in CDS3PrimeFlankingRegionSource( mol, minLength=minLength, debug=debug ):
            if debug:
                print("-----------"*3)
                print(feat)
                print(region)
            
            if feat.data['strand']=='+':
                #begin = region[1]-30
                begin = feat.begin
                end = region[2]
                seq = genomeModel.getRegion( (region[0], begin, end) )
            
            elif feat.data['strand']=='-':
                begin = region[1]  -1
                #end = region[2]+30 -1
                end = feat.end -1
                seq = genomeModel.getRegion( (region[0], begin, end), rc=True )

            else:
                raise ValueError("Unknown strand '{}".format(feature.data['strand']))

            #print(seq[:42])
            ntSeq = Seq( seq.seq, generic_dna )
            if debug:
                print(ntSeq)
                print(len(seq.seq))
            
            #if feat.data['strand']=='-':
            #    ntSeq = ntSeq.reverse_complement()

            aaSeq = ntSeq.translate(table=genomeModel.geneticCode)
            numRegions += 1

            stopCodonPos = feat.end - feat.begin - 3

            if debug:
                print(aaSeq)
                print(aaSeq[((stopCodonPos-6)//3):((stopCodonPos+9)//3):])

            aaAtStopCodonPosition = aaSeq[stopCodonPos//3]
            if not (aaAtStopCodonPosition=="*" or aaAtStopCodonPosition=="X"):
                print("Warning: stop codon not found...")
            
                
            yield (feat, region, ntSeq, stopCodonPos)
            


    print("Found {} regions".format(numRegions))

def writeFlankingSeqToFasta( genomeModel, minLength:int=10, debug=False ):
    #from Bio import SeqIO
    #from Bio.SeqRecord import SeqRecord
    #from Bio.Seq import Seq
    #from Bio.Alphabet import generic_dna

    numRegions = 0

    for (feat,region,ntSeq,stopCodonPos) in CDSWith3PrimeSequencesSource( genomeModel, minLength=minLength, debug=debug ):

            #print(seq[:42])
            if debug:
                print(ntSeq)
                print(len(ntSeq))
            
            #if feat.data['strand']=='-':
            #    ntSeq = ntSeq.reverse_complement()

            aaSeq = ntSeq.translate(table=genomeModel.geneticCode)
            numRegions += 1
            
            if debug:
                print(aaSeq)
            
            #recs.append( SeqRecord( Seq(seq), id=seq.name, name=seq.fancy_name, description=orig.description) )


    print("Found {} regions".format(numRegions))
    #with open(args.output_fasta, "w") as outfile:
    #    out = SeqIO.write( outRecords, outfile, "fasta")

        

def testAll():
    # gm1 = GenomeModel(
    #     sequenceFile='/tamir1/mich1/cellfold/data/Ensembl/Homo.sapiens/Homo_sapiens.GRCh38.dna_rm.toplevel.fa.gz',
    #     gffFile='/tamir1/mich1/cellfold/data/Ensembl/Homo.sapiens/Homo_sapiens.GRCh38.95.gff3.gz',
    #     isLinear=True,
    #     variant="Ensembl",
    #     geneticCode=1 )
                               
    gm2 = GenomeModel(
        sequenceFile='/tamir1/mich1/termfold/data/Ensembl/Ecoli/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_rm.toplevel.fa', # bgzip compression is not supported yet!
        gffFile=     '/tamir1/mich1/termfold/data/Ensembl/Ecoli/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz',
        isLinear=True,
        variant="Ensembl",
        geneticCode=11 )

    print("Forward:")
    testForwardIteration(  gm2.moleculeModels[0] )
    print("Backward:")
    testBackwardIteration( gm2.moleculeModels[0] )

    test3primeFlankingRegions( gm2.moleculeModels[0] )

    #writeFlankingSeqToFasta( gm2, minLength=1, debug=True )
    writeFlankingSeqToFasta( gm2, minLength=1, debug=True )
        
    
    return 0

    

if __name__=="__main__":
    import sys
    sys.exit( testAll() )
