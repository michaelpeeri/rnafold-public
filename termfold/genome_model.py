import json
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
            if end > start:
                self.features.addi( start, end+1, data )  # end is inclusive in GFF3 but not inclusive in intervaltree
            #print(data)
            
            #if len(self.features)>10:
            #    break
        
        
        
class GenomeModel(object):

    _gf_mol_name = 1
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


def testAll():
    gm1 = GenomeModel(
        sequenceFile='/tamir1/mich1/cellfold/data/Ensembl/Homo.sapiens/Homo_sapiens.GRCh38.dna_rm.toplevel.fa.gz',
        gffFile='/tamir1/mich1/cellfold/data/Ensembl/Homo.sapiens/Homo_sapiens.GRCh38.95.gff3.gz',
        isLinear=True,
        variant="Ensembl" )
                               
    gm2 = GenomeModel(
        sequenceFile='/tamir1/mich1/data/Ensembl/Ecoli/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_rm.toplevel.fa.gz',
        gffFile='/tamir1/mich1/termfold/data/Ensembl/Ecoli/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz',
        isLinear=True,
        variant="Ensembl" )
    return 0

    

if __name__=="__main__":
    import sys
    sys.exit( testAll() )
