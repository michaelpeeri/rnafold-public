# For RNAfold SI [running under TERMfold environment!]
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from data_helpers import SpeciesCDSSource, CDSHelper
from genome_model import getGenomeModelFromCache

taxId = 511145
problematicGenes = frozenset(('b2536', 'b0556', 'b1553' ))

outputData  = 'three_prime_intergenic_distances.csv'
outputFasta = 'three_prime_regions.fna'

gm = getGenomeModelFromCache( taxId )

def removeGenePrefix(ident):
    if ident[:5]=="gene:":
        return ident[5:]
    else:
        return ident
    

alreadyProcessedGenes = {}
totalProteinsProcessed = 0
totalSkipped = 0

seqsForWriting=[]

for protId in SpeciesCDSSource(taxId):
    cds = CDSHelper( taxId, protId )
    totalProteinsProcessed += 1

    #feature = gm.findFeatureById( protId )
    geneId = cds.getGeneId()
                
    #flanking3UTRRegionLengthNt = cds.flankingRegion3UtrLength()
    
    feature = gm.findFeatureById( protId )
    #feature = cds.getMatchingFeatureFromGenomeModel()
    #print(feature)
    strand = feature[1].data['strand']

    if strand=='+':
        otherFeature = gm.moleculeModels[ feature[0] ].find5PrimeFlankingRegion( feature[1] )
        
        if otherFeature is None:
            totalSkipped += 1
            continue

        assert( otherFeature['downstream-feature'].begin <= otherFeature['downstream-feature'].end)
        flanking3UTRRegionLengthNt = otherFeature['curr-feature'].begin       -  otherFeature['downstream-feature'].end

        threePrimeUTRCoords = (feature[1].begin-20, feature[1].begin+2, False) # include the first 3 nucleotides of the CDS
        
    else:
        otherFeature = gm.moleculeModels[ feature[0] ].find5PrimeFlankingRegion( feature[1] )

        if otherFeature is None:
            totalSkipped += 1
            continue
        
        assert( otherFeature['downstream-feature'].begin <= otherFeature['downstream-feature'].end)
        flanking3UTRRegionLengthNt = otherFeature['downstream-feature'].begin - otherFeature['curr-feature'].end

        threePrimeUTRCoords = (feature[1].end-3, feature[1].end+20, True) # include the first 3 nucleotides of the CDS

    #assert( flanking3UTRRegionLengthNt > -50 )  # extreme overlaps are not possible

    threePrimeUTR = gm.moleculeModels[ feature[0] ].getSequence( *threePrimeUTRCoords )

    print("{},{},{},{},{}".format( protId, geneId, strand, flanking3UTRRegionLengthNt, threePrimeUTR.seq ))

    seqsForWriting.append( SeqRecord( Seq(threePrimeUTR.seq[:-3], NucleotideAlphabet), id=protId) )
        
    continue

    #if feature[1].data['strand']==downstreamFeature['downstream-feature'].data['strand']:
    

    # downstreamFeature = gm.moleculeModels[ feature[0] ].find3PrimeFlankingRegion( feature[1] )

    # intergenicInfo = []
    # initiationGeneId = None

    # if feature[1].data['strand']==downstreamFeature['downstream-feature'].data['strand']:

    #     thisGeneId = removeGenePrefix( downstreamFeature['curr-feature'].data['gene-id'] )
    #     assert(thisGeneId==geneId)

    #     nextGeneId = removeGenePrefix( downstreamFeature['downstream-feature'].data['gene-id'] )

    #     initiationGeneId = nextGeneId

    #     intergenicInfo = (nextGeneId,
    #                       flanking3UTRRegionLengthNt,
    #                       geneId)

    #     # print("{},{},{},{}".format(protId,
    #     #                            nextGeneId,
    #     #                            flanking3UTRRegionLengthNt,
    #     #                            geneId))

    # else:
    #     upstreamFeature = gm.moleculeModels[ feature[0] ].find5PrimeFlankingRegion( feature[1] )
    #     #usGeneId = removeGenePrefix( upstreamFeature['downstream-feature'].data['gene-id'] )
    #     usProps = json.loads( upstreamFeature['downstream-feature'].data['props'] )
    #     usGeneId = removeGenePrefix( upstreamFeature['downstream-feature'].data['gene-id'] )
    #     usProtId = usProps['protein_id'][0]
        
    #     usCDS = CDSHelper( taxId, usProtId )
        
    #     if feature[1].data['strand']=='+':
    #         usflanking3UTRRegionLengthNt = upstreamFeature['curr-feature'].begin       - upstreamFeature['downstream-feature'].end
    #     else:
    #         usflanking3UTRRegionLengthNt = upstreamFeature['downstream-feature'].begin - upstreamFeature['curr-feature'].end


    #     initiationGeneId = geneId

    #     intergenicInfo = (geneId,
    #                       usflanking3UTRRegionLengthNt,
    #                       usGeneId)
        
    #     # print("{},{},{},{}".format(protId,
    #     #                            geneId,
    #     #                            usflanking3UTRRegionLengthNt,
    #     #                            usGeneId))

    # initiationGene = gm.findFeatureById( initiationGeneId )
    # strand = initiationGene[1].data['strand']
    # if strand=='+':
    #     threePrimeUTRCoords = (initiationGene[1].begin-20, initiationGene[1].begin+2, False)
    # else:
    #     threePrimeUTRCoords = (initiationGene[1].end-3, initiationGene[1].end+20, True)

    # threePrimeUTR = gm.moleculeModels[ initiationGene[0] ].getSequence( *threePrimeUTRCoords )

    # dataForOutput = intergenicInfo + threePrimeUTRCoords + (threePrimeUTR.seq,)

    # if initiationGeneId in problematicGenes:
    #     continue
    
    # if initiationGeneId in alreadyProcessedGenes:
    #     #raise Exception("Gene {} already processed".format( initiationGeneId ))
    #     if( dataForOutput != alreadyProcessedGenes[initiationGeneId] ):
    #         print("Old --> {},{},{},{},{},{},{}".format( *(alreadyProcessedGenes[initiationGeneId]) ))
    #         print("New --> {},{},{},{},{},{},{}".format( *dataForOutput ))
    #     totalSkipped += 1
    #     #pass
    
    # else:
    #     alreadyProcessedGenes[ initiationGeneId ] = dataForOutput + (protId,)
    #     #print("{},{},{},{},{},{},{}".format(*dataForOutput))


SeqIO.write( seqsForWriting, outputFasta, "fasta")  # write the full sequences into the file

print("Processed {} coding sequences".format( totalProteinsProcessed ))
print("Skipped {} coding sequences".format( totalSkipped ))

