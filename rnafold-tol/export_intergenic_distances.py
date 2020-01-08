# For RNAfold SI [running under TERMfold environment!]
import json
import subprocess
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from data_helpers import SpeciesCDSSource, CDSHelper
from genome_model import getGenomeModelFromCache
import config


outputData  = 'three_prime_intergenic_distances_{}.csv'
outputFasta = 'three_prime_regions_{}.fna'


# def removeGenePrefix(ident):
#     if ident[:5]=="gene:":
#         return ident[5:]
#     else:
#         return ident


reResultLine = re.compile("(\w+)\s[<]unknown description[>]\s([-]?\d+[.]\d+)")
def calculateaSDEnergies(seqs, args, taxId):
    fastafn = outputFasta.format(taxId)
    SeqIO.write( seqs, fastafn, "fasta")  # write the full sequences into the file

    "matlab -nodisplay -nodesktop -nosplash -singleCompThread -r \"calculate_aSD('{}');\"".format( fastafn )
    cmdline = (config.MatlabPath, "-nodisplay", "-nodesktop", "-nosplash", "-singleCompThread", "-batch",
                "\"calculate_aSD('{}');quit()\"".format(fastafn))
    #cmdline = (config.MatlabPath, "-nodisplay", "-nodesktop", "-nosplash", "-singleCompThread", "-batch",
    #            "\"disp('{}');quit()\"".format(fastafn))
    print(" ".join(cmdline))
    out = subprocess.check_output(" ".join(cmdline), shell=True)
    #print(out)

    ret = {}
    lineCount = 0
    for line in str(out, encoding="ascii").split('\n'):
        if not line: continue
        lineCount += 1
        m = reResultLine.match(line)
        if not m:
            print(line)
            continue
        
        protId = m.group(1)
        assert(len(protId) > 2)
        aSD = float(m.group(2))
        assert(aSD <= 0.0)
        ret[protId] = aSD
    print("processed {} lines".format(lineCount))
        
    return ret
                    

    
def processGenome(args, taxId):

    alreadyProcessedGenes = {}
    totalProteinsProcessed = 0
    totalSkipped = 0

    seqsForWriting=[]
    recordsForWriting={}
    
    gm = getGenomeModelFromCache( taxId )

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

        threePrimeUTR = gm.moleculeModels[ feature[0] ].getSequence( *threePrimeUTRCoords )

        if flanking3UTRRegionLengthNt < -50:
            print("Warning: found gene with apparent long overlap: {},{},{},{},{}".format( protId, geneId, strand, flanking3UTRRegionLengthNt, threePrimeUTR.seq ))
            #totalSkipped += 1
            #continue

        if threePrimeUTR.seq[-2:] != 'TG':
            print("Warning: skipping gene with start codon at the correct place: {},{},{},{},{}".format( protId, geneId, strand, flanking3UTRRegionLengthNt, threePrimeUTR.seq ))
            totalSkipped += 1
            continue

        # All done - emit the output
        #fout.write("{},{},{},{},{}".format( protId, geneId, strand, flanking3UTRRegionLengthNt, threePrimeUTR.seq ))
        recordsForWriting[protId] = (geneId, strand, flanking3UTRRegionLengthNt, threePrimeUTR.seq )

        seqsForWriting.append( SeqRecord( Seq(threePrimeUTR.seq[:-3], NucleotideAlphabet), id=protId) )

    aSD = calculateaSDEnergies( seqsForWriting, args, taxId )
    print(len(aSD))

    with open( outputData.format(taxId), 'wt') as fout:
        for protId, record in recordsForWriting.items():
            aSDval = aSD.get(protId, None)
            vals = (protId,) + record + (aSDval,)
            fout.write("{},{},{},{},{},{}\n".format( *vals ))
    

    print("Processed {} coding sequences for taxid {}".format( totalProteinsProcessed, taxId ))
    print("Skipped {} coding sequences".format( totalSkipped ))



def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert

if __name__=="__main__":
    import sys
    import argparse

    argsParser = argparse.ArgumentParser()
    argsParser.add_argument( "--taxid",           type=parseList(int) )
    args = argsParser.parse_args()

    for taxId in args.taxid:
        assert(int(taxId))
        #assert(taxId in speciesConfig)
        #config = speciesConfig[ taxId ]
        processGenome(args, taxId)
    
    sys.exit(0)
