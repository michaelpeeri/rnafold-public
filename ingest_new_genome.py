import os
import sys
import re
import subprocess
from ensembl_ftp import EnsemblFTP
from data_helpers import r, getSpeciesName, speciesNameKey, speciesTaxIdKey, speciesTranslationTableKey, speciesCDSList
from config import ensembl_data_dir

# ~/anaconda2/bin/python2 test_gff_extract_coding_genes.py --gff ./data/Ensembl/Twhipplei/Tropheryma_whipplei_str_twist.ASM748v1.36.gff3.gz --variant Ensembl --output-gene-ids ./data/Ensembl/coding_genes_list.Twhi.list --transl-table 11
def processGff3(gff3Filename, genesListFilename, args):
    print("Processing gff3 file...")
    x = subprocess.check_output((sys.executable, "test_gff_extract_coding_genes.py",
                                 "--gff", gff3Filename,
                                 "--variant", "Ensembl",
                                 "--output-gene-ids", genesListFilename,
                                 "--transl-table", str(args.nuclear_genetic_code) ), shell=False)


# ~/anaconda2/bin/python2 store_seqs_from_fasta.py --taxid 203267 --input ./data/Ensembl/Twhipplei/Tropheryma_whipplei_str_twist.ASM748v1.cds.all.fa.gz --variant Ensembl --type cds --output-fasta ./data/Ensembl/Twhipplei/Tropheryma_whipplei_str_twist.ASM748v1.cds.all.fa.gz.filtered.fna --gene-ids-file ./data/Ensembl/coding_genes_list.Twhi.list --dry-run
reSequenceOkLine = re.compile("Inserting (\S+) [(]sequence \S+[)]...")
reEntryCountLine = re.compile("Processed (\d+) CDS entries")
reTotalCountLine = re.compile("[(]out of \d+ CDS entries for this species[)]")
reConnectionRecycling = re.compile(".*sqlalchemy.pool.QueuePool.*exceeded timeout; recycling")
reSkippingAmbiguousSeq = re.compile("Skipping record (\S+), containing non-nucleotide or ambiguous symbols.*")
reWarningEntriesSkipped = re.compile("Warning: (\d+) entries skipped and (\d+) entries not found")
def loadCDSSequences(cdsSequencesFilename, genesListFilename, args, dryRun=True):
    print("Reading CDS sequences file...")

    arguments = (sys.executable, "store_seqs_from_fasta.py",
                 "--taxid", str(args.taxid),
                 "--input", cdsSequencesFilename,
                 "--variant", "Ensembl",
                 "--type", "cds",
                 "--output-fasta", "%s.filtered.fna" % cdsSequencesFilename,
                 "--gene-ids-file", genesListFilename)
    if dryRun:
        arguments = arguments + ("--dry-run",)

    report = subprocess.check_output(arguments, shell=False)
    
    unknownLinesCount = 0
    loadedGenesCount = 0
    loadedGenes = []
    skippedGenes = []
    for line in report.splitlines():
        # Note: Order cases from most to least common (for performance)
        match = reSequenceOkLine.match(line)
        if match:
            loadedGenes.append(match.group(1))
            continue
        
        elif reEntryCountLine.match(line):
            loadedGenesCount = int(reEntryCountLine.match(line).group(1))
            print(line)
            continue
        
        elif reTotalCountLine.match(line):
            print(line)
            continue
        
        elif reSkippingAmbiguousSeq.match(line):
            skippedGenes.append( reSkippingAmbiguousSeq.match(line).group(1) )
            print(line)
            continue
        
        elif reConnectionRecycling.match(line):
            continue
        
        elif reWarningEntriesSkipped.match(line):
            print(line)
            continue
        
        else:
            unknownLinesCount+= 1  # Count "unexpected" lines
            print(line)
            continue
        
    if unknownLinesCount:
        return None
    else:
        return (loadedGenesCount, loadedGenes, skippedGenes)
    
def ingestGenome(args):

    # Sanity test 1 -- genes list file doesn't already exists (if it does, short name may have been reused...)
    genesListFilename = "%s/coding_genes_list.%s.list" % (ensembl_data_dir, args.short_name)
    assert(not os.path.exists(genesListFilename))

    # Sanity test 2 -- this species doesn't already exist in the DB
    if r.exists(speciesNameKey % args.taxid):
        if r.exists(speciesCDSList % args.taxid):
            raise Exception("Species with taxid=%d already exists as '%s'" % (args.taxid, getSpeciesName(args.taxid)))
        else:
            pass # Name defined but no CDS list; maybe an earlier run of this script was interrupted?
    

    # Sanity tests passed
    
    # Step 1 - get files from Ensembl FTP
    
    ftp = EnsemblFTP(args.local_name, args.remote_name, release=args.release, subsection=args.subsection)
    (genomefn, cdsfn, gff3fn) = ftp.fetchAll()
    ftp.close()

    assert( os.path.exists(genomefn) and os.path.isfile(genomefn) )
    assert( os.path.exists(cdsfn)    and os.path.isfile(cdsfn) )
    assert( os.path.exists(gff3fn)   and os.path.isfile(gff3fn) )

    # Step 2 - parse GFF3 file to yield list of acceptable CDS genes
    
    processGff3(gff3fn, genesListFilename, args)

    numGenesReturnedFromGff3 = None
    with open(genesListFilename, "r") as f:
        numGenesReturnedFromGff3 = len(f.readlines())
    if numGenesReturnedFromGff3 < 400:
        raise Exception("Processing gff3 file only yielded %d results; aborting" % numGenesReturnedFromGff3)
    print("%d protein-coding genes found in gff3" % numGenesReturnedFromGff3)

    # Step 3 - Add required annotations for this species to redis DB
    
    # TODO - add redis items here
    #redis-cli -h power5 -a rnafold set "species:taxid:203267:name" "Tropheryma whipplei str. Twist"
    r.set(speciesNameKey % args.taxid, args.full_name)

    #redis-cli -h power5 -a rnafold set "species:name:Tropheryma whipplei str. Twist:taxid" "203267"
    r.set(speciesTaxIdKey % args.full_name, args.taxid)

    #redis-cli -h power5 -a rnafold set "species:taxid:203267:genomic-transl-table"  "11"
    r.set(speciesTranslationTableKey % args.taxid, args.nuclear_genetic_code)

    # Step 4 - Load CDS sequences to DB
    print("Doing trial run for gene loading...")
    if not loadCDSSequences(cdsfn, genesListFilename, args, dryRun=True) is None:
        print("Trial run succeeded.")
        print("Performing actual gene loading...")
        (cdsLoadedCount, cdsIds, skippedGenes) = loadCDSSequences(cdsfn, genesListFilename, args, dryRun=False)
        print("Loaded %d CDS genes..." % cdsLoadedCount)
        if skippedGenes:
            print("Skipped genes: %s" % skippedGenes)
    else:
        print("Dry-run loading may have encountered errors; aborting without actual load...")
        return -1

    # Step 5 - Generate randomized sequences
    # TODO
    
    return 0


def standaloneRun():
    import argparse
    import sys
    
    argsParser = argparse.ArgumentParser()
    #argsParser.add_argument("--verbose", type=int, default=0)
    argsParser.add_argument("--local-name", type=str, required=True)
    argsParser.add_argument("--remote-name", type=str, required=True)
    argsParser.add_argument("--release", type=int, default=36)
    argsParser.add_argument("--subsection", type=str, default="bacteria_3_collection")
    argsParser.add_argument("--nuclear-genetic-code", type=int, required=True)
    argsParser.add_argument("--taxid", type=int, required=True)
    argsParser.add_argument("--full-name", type=str, required=True)
    argsParser.add_argument("--short-name", type=str, required=True)
    args = argsParser.parse_args()

    sys.exit(ingestGenome(args))
    
    
if __name__=="__main__":
    standaloneRun()
