from __future__ import print_function
import gff

print("Loading gff file...")
db = gff.createGffDb("data/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff")
print("Done.")


# Collect all Yeast coding genes; skip 'Dubious' ORFs
includedGenes = 0
skippedGenes = 0

for area in ('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrmt'):
    codingGenes = set()
    for cds in db.features_of_type('CDS', limit=(area, 0, 100000000), completely_within=False):
        for mRNA in db.parents(cds.id):
            for gene in db.parents(mRNA.id):
                if gene.featuretype == 'gene':
                    if gene.attributes['orf_classification'][0] == 'Dubious':
                        skippedGenes += 1
                        continue

                    codingGenes.add(gene.id)
    
    print(area, len(codingGenes))
    includedGenes += len(codingGenes)
    print(list(codingGenes)[1:10])

print("%d genes included, %d skipped" % (includedGenes, skippedGenes))
