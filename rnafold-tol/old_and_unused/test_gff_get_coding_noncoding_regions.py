from __future__ import print_function
from intervaltree import IntervalTree    # https://pypi.python.org/pypi/intervaltree
import gff

print("Loading gff file...")
db = gff.createGffDb("data/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff")
print("Done.")

yeastChromosomes = ('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrmt')

# Collect all Yeast coding genes; skip 'Dubious' ORFs
includedGenes = 0
skippedGenes = 0
errorCount = 0

for area in yeastChromosomes:
    codingRegions = IntervalTree()
    nonCodingRegions = IntervalTree()

    transcribedRegions = IntervalTree()
    nonTranscribedRegions = IntervalTree()
    selectedChromosome = (area, 0, 100000000)

    chromosomeBoundaries = None
    for chromosome in db.features_of_type('chromosome', limit=selectedChromosome, completely_within=False ):
        chromosomeBoundaries = (chromosome.start, chromosome.end)

    assert(chromosomeBoundaries[0] < chromosomeBoundaries[1] - 1000)
    assert(chromosomeBoundaries[1] <= selectedChromosome[2])   # Make sure our arbitrary limit is suitable

    nonCodingRegions.addi( chromosomeBoundaries[0], chromosomeBoundaries[1], area)
    nonTranscribedRegions.addi( chromosomeBoundaries[0], chromosomeBoundaries[1], area)

    # -----------------------------------------------------
    # Collect coding regions
    # -----------------------------------------------------
    # struture for blocked_reading_frames:              blocked_reading_frame -> CDS
    # TODO - filter pseudogenes !!!!
    for cds in db.features_of_type('CDS', limit=selectedChromosome, completely_within=False):
        assert(cds.featuretype == 'CDS')

        geneId = None

        for feat in db.parents(cds.id):
            if feat.featuretype == 'mRNA':
                pass
            elif feat.featuretype == 'gene':
                geneId = feat.id
            elif feat.featuretype == 'transposable_element_gene':
                geneId = feat.id
            elif feat.featuretype == 'tRNA_gene':
                print("---Warning!----")
                continue
            elif feat.featuretype == 'ncRNA_gene':
                print("---Warning!----")
                continue
            elif feat.featuretype == 'blocked_reading_frame':
                geneId = feat.id
            else:
                print(feat.featuretype)
                assert(False)

        if( geneId==None ):
            print("Warning (coding)")
            print(cds.id)
            print(cds.attributes)
            continue  # Skip this CDS
        assert( geneId != None )

        if( cds.end - cds.start < 1 ):
            continue  # Skip this CDS

        codingRegions.addi( cds.start, cds.end, geneId )
        nonCodingRegions.chop( cds.start, cds.end )

    # -----------------------------------------------------
    # Collect transcribed regions
    # -----------------------------------------------------
    # Note: standard coding genes have the general structure: gene -> mRNA -> CDS
    #       non-coding genes have the general structure       XXXX_gene -> noncoding_exon   (XXXX can be tRNA, rRNA, snRNA, snoRNA, ncRNA)
    for mRNA in db.features_of_type(('mRNA', 'noncoding_exon'), limit=selectedChromosome, completely_within=False):
        geneId = None
        for gene in db.parents(mRNA.id):
            if gene.featuretype == 'gene':
                geneId = gene.id
            elif gene.featuretype == 'transposable_element_gene' or gene.featuretype == 'LTR_retrotransposon' or gene.featuretype=='tRNA_gene' or gene.featuretype=='ncRNA_gene' or gene.featuretype=='snoRNA_gene' or gene.featuretype=='rRNA_gene' or gene.featuretype=='snRNA_gene' or gene.featuretype=='telomerase_RNA_gene':
                # TODO - What to do about transposable element genes?!
                geneId = gene.id
            else:
                print(":", gene.featuretype)
        
        if( geneId==None ):
            print("Warning")
            # TODO - What to do about transposable element genes?!
            continue

        if( mRNA.end - mRNA.start < 1 ):
            continue

        transcribedRegions.addi( mRNA.start, mRNA.end, geneId )
        nonTranscribedRegions.chop( mRNA.start, mRNA.end )

    # Add transposable elements as transcribed regions
    # transposable elements have the structure: transposable_element -> transposable_elemenet_gene -> CDS   (with no mRNA annotation).
    #       struture for blocked_reading_frames:              blocked_reading_frame -> CDS
    for feat in db.features_of_type(('transposable_element_gene', 'blocked_reading_frame'), limit=selectedChromosome, completely_within=False):
        geneId = feat.id
        assert(geneId != None)

        if( feat.end - feat.start < 1 ):
            continue

        transcribedRegions.addi( feat.start, feat.end, geneId )
        nonTranscribedRegions.chop( feat.start, feat.end )
    
    print("================ [", area, "] ================")
    a = []
    b = []
    for i in range(1,chromosomeBoundaries[1], chromosomeBoundaries[1]/143):
        a.append("#" if codingRegions.overlaps(i) else "-")
        b.append("#" if nonCodingRegions.overlaps(i) else "-")
        assert(a[-1] != b[-1])
    #print("".join(b))

    #print(area, len(transcribedRegions), len(nonTranscribedRegions))
    c = []
    d = []
    for i in range(1,chromosomeBoundaries[1], chromosomeBoundaries[1]/143):
        c.append("#" if transcribedRegions.overlaps(i) else "-")
        d.append("#" if nonTranscribedRegions.overlaps(i) else "-")
        assert(c[-1] != d[-1])
        # assert( non-coding | transcribed )
        if( not ((a[len(c)-1]=="-") | (c[-1]=="#")) ):
            print("--------------------------------------------")
            print("Found coding-region anomaly")
            print(area, i)
            print("Coding: ", a[len(c)-1])
            print(list(codingRegions[i-2500:i+2500]))
            print("Transcribed: ", c[-1])
            print(list(transcribedRegions[i-2500:i+2500]))
            errorCount += 1

        #assert((a[len(c)-1]=="-") | (c[-1]=="#"))
    print("Transcribed: 5'", "".join(c), "3'")
    #print("".join(d))

    print("Coding:      5'", "".join(a), "3'")


print("%d genes included, %d skipped" % (includedGenes, skippedGenes))
print("Found %d errors" % errorCount)
