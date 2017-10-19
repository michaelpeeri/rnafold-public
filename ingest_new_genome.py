import os
import sys
import re
import subprocess
from ensembl_ftp import EnsemblFTP
from data_helpers import r, getSpeciesName, speciesNameKey, speciesTaxIdKey, speciesTranslationTableKey, speciesCDSList
from config import ensembl_data_dir





# ~/anaconda2/bin/python2 gff3_extract_coding_genes.py --gff ./data/Ensembl/Twhipplei/Tropheryma_whipplei_str_twist.ASM748v1.36.gff3.gz --variant Ensembl --output-gene-ids ./data/Ensembl/coding_genes_list.Twhi.list --transl-table 11
def processGff3(gff3Filename, genesListFilename, args):
    print("Processing gff3 file...")
    output = []
    try:
        output = subprocess.check_output((sys.executable, "gff3_extract_coding_genes.py",
                                          "--gff3", gff3Filename,
                                          "--variant", args.variant, #"Ensembl",
                                          "--output-gene-ids", genesListFilename,
                                          "--transl-table", str(args.nuclear_genetic_code) ), shell=False)
    except subprocess.CalledProcessError as e:
        print(e)
        print("--" * 20)
        print("Process output: ")
        print("--" * 20)
        for line in output:
            print(line)
        print("--" * 20)
        print(e)


# ~/anaconda2/bin/python2 store_seqs_from_fasta.py --taxid 203267 --input ./data/Ensembl/Twhipplei/Tropheryma_whipplei_str_twist.ASM748v1.cds.all.fa.gz --variant Ensembl --type cds --output-fasta ./data/Ensembl/Twhipplei/Tropheryma_whipplei_str_twist.ASM748v1.cds.all.fa.gz.filtered.fna --gene-ids-file ./data/Ensembl/coding_genes_list.Twhi.list --dry-run
reSequenceOkLine = re.compile("Inserting (\S+) [(]sequence \S+[)]...")
reEntryCountLine = re.compile("Processed (\d+) CDS entries")
reTotalCountLine = re.compile("[(]out of \d+ CDS entries for this species[)]")
reConnectionRecycling = re.compile(".*sqlalchemy.pool.QueuePool.*exceeded timeout; recycling")
reSkippingAmbiguousSeq = re.compile("Skipping record (\S+), containing non-nucleotide or ambiguous symbols.*")


#Skipping OEU05659.1 (sequence lcl|KV784598.1_cds_OEU05659.1_18098, alternate ids=[])
reSkippingExcludedSeq = re.compile("Skipping (\S+) [(]sequence (\S+), alternate ids=[[][^]]*[]][)]")
#Skipping pseudo-gene entry lcl|KV784402.1_cds_17051

#"Skipping pseudo-gene entry (\S+)"
reSkippingPseudoGene = re.compile("Skipping pseudo-gene entry (\S+)")

reWarningEntriesSkipped = re.compile("Warning: (\d+) entries skipped and (\d+) entries not found")
def loadCDSSequences(cdsSequencesFilename, genesListFilename, args, dryRun=True):
    print("Reading CDS sequences file...")

    arguments = (sys.executable, "store_seqs_from_fasta.py",
                 "--taxid", str(args.taxid),
                 "--input", cdsSequencesFilename,
                 "--variant", args.variant,
                 "--type", "cds",
                 "--output-fasta", "%s.filtered.fna" % cdsSequencesFilename,
                 "--gene-ids-file", genesListFilename)
    if dryRun:
        arguments = arguments + ("--dry-run",)

    if args.ignore_id_check:
        arguments = arguments + ("--ignore-id-check",)

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

        elif reSkippingExcludedSeq.match(line):
            skippedGenes.append( reSkippingExcludedSeq.match(line).group(1) )
            print(line)
            continue

        elif reSkippingPseudoGene.match(line):
            skippedGenes.append( reSkippingPseudoGene.match(line).group(1) )
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
    print(genesListFilename)
    assert(not os.path.exists(genesListFilename))

    # Sanity test 2 -- this species doesn't already exist in the DB
    if r.exists(speciesNameKey % args.taxid):
        if r.exists(speciesCDSList % args.taxid):
            raise Exception("Species with taxid=%d already exists as '%s'" % (args.taxid, getSpeciesName(args.taxid)))
        else:
            pass # Name defined but no CDS list; maybe an earlier run of this script was interrupted?
    

    # Sanity tests passed

    genomefn = None
    cdsfn = None
    gff3fn = None
    
    # Step 1 - get files from Ensembl FTP

    if( args.variant == "Ensembl" and args.fetch_ftp_files ):
    
        ftp = EnsemblFTP(args.local_name, args.remote_name, release=args.release, section=args.section, subsection=args.subsection)
        (genomefn, cdsfn, gff3fn) = ftp.fetchAll()
        ftp.close()

        assert( os.path.exists(genomefn) and os.path.isfile(genomefn) )
        assert( os.path.exists(cdsfn)    and os.path.isfile(cdsfn) )
        assert( os.path.exists(gff3fn)   and os.path.isfile(gff3fn) )
    else:
        gff3fn = args.gff3
        cdsfn = args.cds

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

_knownBoolVals = {"true":True, "false":False}
def parseBool(val):
    _val = val.lower()
    if _val in _knownBoolVals:
        return _knownBoolVals[_val]
    else:
        raise Exception("Unknown bool value '%s'" % _val)
    

def standaloneRun():
    import argparse
    import sys
    
    argsParser = argparse.ArgumentParser()
    #argsParser.add_argument("--verbose", type=int, default=0)
    argsParser.add_argument("--local-name", type=str, required=False)
    argsParser.add_argument("--remote-name", type=str, required=False)
    argsParser.add_argument("--release", type=int, default=37)
    argsParser.add_argument("--section", type=str, default="bacteria")
    argsParser.add_argument("--subsection", type=str, required=False)
    argsParser.add_argument("--nuclear-genetic-code", type=int, required=True)
    argsParser.add_argument("--taxid", type=int, required=True)
    argsParser.add_argument("--full-name", type=str, required=True)
    argsParser.add_argument("--short-name", type=str, required=True)
    argsParser.add_argument("--variant", type=str, default="Ensembl")
    argsParser.add_argument("--fetch-ftp-files", type=parseBool, default=True)
    argsParser.add_argument("--gff3", type=str, required=False)
    argsParser.add_argument("--cds", type=str, required=False)
    argsParser.add_argument("--ignore-id-check", action="store_true", default=False)
    
    args = argsParser.parse_args()

    print(args)

    print(args.fetch_ftp_files)

    if( args.fetch_ftp_files ):
        if( (not args.local_name) or (not args.remote_name) ):
            raise Exception("--remote-name and --local-name are required when --fetch-ftp-files=True")
        
        if( args.gff3 or args.cds ):
            raise Exception("--gff3 and --cds are not supported when --fetch-ftp-files=True")
        
        if args.variant != "Ensembl":
            raise Exception("Only --variant=Ensembl is supported with --fetch-ftp-files=True")
    else:
        if( args.local_name or args.remote_name ):
            raise Exception("--remote-name and --local-name cannot be used when --fetch-ftp-files=False")
        
        if( (not args.gff3) or (not args.cds) ):
            raise Exception("--gff3 and --cds are required when --fetch-ftp-files=False")
        

    sys.exit(ingestGenome(args))
    
    
if __name__=="__main__":
    standaloneRun()
