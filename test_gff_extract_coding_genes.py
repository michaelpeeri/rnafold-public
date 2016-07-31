from __future__ import print_function
import argparse
import gff



def parseOption(possibleValues, name):
    def checkOption(value):
        if value in possibleValues:
            return value
        else:
            raise argparse.ArgumentTypeError("Unknown %s '%s', allowed values: %s" % (name, value, ",".join(possibleValues)))
    return checkOption

argsParser = argparse.ArgumentParser()
#argsParser.add_argument("--taxid", type=int)
argsParser.add_argument("--gff")
argsParser.add_argument("--variant", type=parseOption(set(("yeastgenome","NCBI")), "variant"))
argsParser.add_argument("--output-fasta")
args = argsParser.parse_args()


print("Loading gff file...")
db = gff.createGffDb(args.gff)
print("Done.")


# Collect all Yeast coding genes; skip 'Dubious' ORFs
includedGenes = 0
skippedGenes = 0

codingGenes = set()

prev = 0

for area in ('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI'):     #, 'chrmt'):
    codingGenesInArea = 0
    for cds in db.features_of_type('CDS', limit=(area, 0, 100000000), completely_within=False):
        for mRNA in db.parents(cds.id):
            for gene in db.parents(mRNA.id):
                if gene.featuretype == 'gene':
                    if args.variant=="yeastgenome.org" and gene.attributes['orf_classification'][0] == 'Dubious':
                        skippedGenes += 1
                        continue

                    #if gene.id in codingGenes:
                    #    print("%s already in!" % (gene.id,))
                    # Note - manually tested that duplicates have identical coordinate

                    codingGenes.add(gene.id)
                    codingGenesInArea += 1

    includedGenes += codingGenesInArea
    print(codingGenesInArea, len(codingGenes)-prev)
    prev = len(codingGenes)
    
    print(area, codingGenesInArea)
    #print(list(codingGenes)[1:10])

print("%d genes included, %d skipped" % (includedGenes, skippedGenes))

print(len(codingGenes))

# Write file containing selected gene-ids:
if args.variant=="yeastgenome.org":
    with open("yeast_coding_genes.txt", "w") as out:
        for geneId in codingGenes:
            out.write("%s\n" % geneId)
