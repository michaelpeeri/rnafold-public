from ftplib import FTP
import os
import re
from config import ensembl_data_dir


class LinesBuffer(object):
    def __init__(self):
        self._lines = []
        
    def __call__(self, line):
        self._lines.append(line)

    def get(self):
        return self._lines

"""
drwxr-sr-x    2 ftp      ftp          4096 May 18 23:15 cdna
-rw-r--r--    1 ftp      ftp         11320 May 18 21:20 README
"""
reDirListingLine = re.compile("([d-])[r-][w-][xs-][r-][w-][xs-][r-][w-][xs-]\s+\d+\s+\w+\s+\w+\s+\d+\s+\w+\s+\d+\s+\d+[:]\d+([:]\d+)?\s+(\w+.*)")
def parseDirListingLine_returnFilenames(line):
    match = reDirListingLine.match(line)
    if match is None:
        raise Exception("Unknown format for dir listing line: %s" % line)
    isDir = True if match.group(1)=="d" else False
    name = match.group(3)
    
    return (isDir, name)
    

class EnsemblFTP(object):
    def __init__(self, localDir, speciesDirName, release=37, section="bacteria", subsection=None):
        self._release = release
        self._section = section
        self._subsection = subsection
        self._speciesDirName = speciesDirName

        self._localDir = "%s/%s" % (ensembl_data_dir, localDir)
        self.prepareLocalDir()
        
        self._ftp = FTP("ftp.ensemblgenomes.org")
        self._ftp.login()

    def prepareLocalDir(self):
        if os.path.exists(self._localDir):
            if os.path.isdir(self._localDir):
                for fn in os.listdir(self._localDir):
                    raise Exception("Data path %s not empty, refusing to overwrite..." % self._localDir)
                
                return 
            else:
                raise Exception("Data path %s exists, but is not a directory..." % self._localDir)
        else:
            os.mkdir(self._localDir)

    def getLocalFilename(self, remoteName):
        return "%s/%s" % (self._localDir, remoteName)

    def getDirName(self, category="fasta", subdir=None):
        #ftp.dir("/pub/bacteria/release-36/fasta/bacteria_3_collection/wolinella_succinogenes_dsm_1740", x)
        #ftp.dir("/pub/protists/release-36/fasta/protists_stramenopiles1_collection/thalassiosira_oceanica_ccmp1005", x)
        if subdir is None:
            return "/pub/%s/release-%d/%s/%s/%s/"    % (self._section, self._release, category, self._subsection, self._speciesDirName)
        else:
            return "/pub/%s/release-%d/%s/%s/%s/%s/" % (self._section, self._release, category, self._subsection, self._speciesDirName, subdir)
        

    def listSpeciesItems(self, category="fasta", subdir=None):
        names = LinesBuffer()
        dirName = self.getDirName(category, subdir)
        
        self._ftp.dir(dirName, names)
                                                                 
        return list(map(parseDirListingLine_returnFilenames, names.get()))

    def determineCDSPath(self):
        items = self.listSpeciesItems("fasta")
        try:
            if( items.index((True, "cds")) >= 0 ):
                return "cds"
        except ValueError:
            pass
        raise Exception("Failed to determine 'cds' path")

    def determineGenomePath(self):
        items = self.listSpeciesItems("fasta")

        try:
            if( items.index((True, "dna")) >= 0 ):
                return "dna"
        except ValueError:
            pass
        raise Exception("Failed to determine 'dna' path")

    def fetchGenomeFiles(self):
        genomePath = self.determineGenomePath()
        items = self.listSpeciesItems("fasta", genomePath)

        genomeMatches    = list(filter( lambda x:x[1].endswith("dna_rm.toplevel.fa.gz"), items))
        readmeMatches    = list(filter( lambda x:x[1] == "README"                      , items))
        checksumsMatches = list(filter( lambda x:x[1] == "CHECKSUMS"                   , items))

        localGenomeFilename = self.getLocalFilename(genomeMatches[0][1])

        if genomeMatches:
            path = "%s%s" % (self.getDirName("fasta", genomePath), genomeMatches[0][1])
            print("Fetching genome from %s..." % path)
            with open(localGenomeFilename, "wb") as f:
                resp = self._ftp.retrbinary("RETR %s" % path, f.write)
                if not resp.startswith("226 "):
                    print("Warning: transfer failed? (response: %s)" % resp)
        else:
            raise Exception("Couldn't find genome file")

        if readmeMatches:
            path = "%s%s" % (self.getDirName("fasta", genomePath), readmeMatches[0][1])
            print("Fetching readme from %s..." % path)
            with open(self.getLocalFilename(readmeMatches[0][1]) + ".genome", "wb") as f:
                resp = self._ftp.retrbinary("RETR %s" % path, f.write)
                if not resp.startswith("226 "):
                    print("Warning: transfer failed? (response: %s)" % resp)
        
        if checksumsMatches:
            path = "%s%s" % (self.getDirName("fasta", genomePath), checksumsMatches[0][1])
            print("Fetching checksums from %s..." % path)
            with open(self.getLocalFilename(checksumsMatches[0][1] + ".genome"), "wb") as f:
                resp = self._ftp.retrbinary("RETR %s" % path, f.write)
                if not resp.startswith("226 "):
                    print("Warning: transfer failed? (response: %s)" % resp)

        return localGenomeFilename
    

    def fetchCDSFiles(self):
        cdsPath = self.determineCDSPath()
        items = self.listSpeciesItems("fasta", cdsPath)

        cdsMatches       = list(filter( lambda x:x[1].endswith(".cds.all.fa.gz")       , items))
        readmeMatches    = list(filter( lambda x:x[1] == "README"                      , items))
        checksumsMatches = list(filter( lambda x:x[1] == "CHECKSUMS"                   , items))

        localCDSfilename = self.getLocalFilename(cdsMatches[0][1])
        
        if cdsMatches:
            path = "%s%s" % (self.getDirName("fasta", cdsPath), cdsMatches[0][1])
            print("Fetching CDSs from %s..." % path)
            with open(localCDSfilename, "wb") as f:
                resp = self._ftp.retrbinary("RETR %s" % path, f.write)
                if not resp.startswith("226 "):
                    print("Warning: transfer failed? (response: %s)" % resp)
        else:
            raise Exception("Couldn't find CDS file")

        if readmeMatches:
            path = "%s%s" % (self.getDirName("fasta", cdsPath), readmeMatches[0][1])
            print("Fetching readme from %s..." % path)
            with open(self.getLocalFilename(readmeMatches[0][1]) + ".cds", "wb") as f:
                resp = self._ftp.retrbinary("RETR %s" % path, f.write)
                if not resp.startswith("226 "):
                    print("Warning: transfer failed? (response: %s)" % resp)
        
        if checksumsMatches:
            path = "%s%s" % (self.getDirName("fasta", cdsPath), checksumsMatches[0][1])
            print("Fetching checksums from %s..." % path)
            with open(self.getLocalFilename(checksumsMatches[0][1]) + ".cds", "wb") as f:
                resp = self._ftp.retrbinary("RETR %s" % path, f.write)
                if not resp.startswith("226 "):
                    print("Warning: transfer failed? (response: %s)" % resp)

        return localCDSfilename
        

    def fetchGFF3Files(self):
        items = self.listSpeciesItems("gff3")

        gff3Matches      = list(filter( lambda x:x[1].endswith(".%d.gff3.gz" % self._release), items))
        readmeMatches    = list(filter( lambda x:x[1] == "README"                            , items))
        checksumsMatches = list(filter( lambda x:x[1] == "CHECKSUMS"                         , items))

        localGFF3Filename = self.getLocalFilename(gff3Matches[0][1])

        if gff3Matches:
            path = "%s%s" % (self.getDirName("gff3"), gff3Matches[0][1])
            print("Fetching GFF3 from %s..." % path)
            with open(localGFF3Filename, "wb") as f:
                resp = self._ftp.retrbinary("RETR %s" % path, f.write)
                if not resp.startswith("226 "):
                    print("Warning: transfer failed? (response: %s)" % resp)
        else:
            raise Exception("Couldn't find CDS file")

        if readmeMatches:
            path = "%s%s" % (self.getDirName("gff3"), readmeMatches[0][1])
            print("Fetching readme from %s..." % path)
            with open(self.getLocalFilename(readmeMatches[0][1]) + ".gff3", "wb") as f:
                resp = self._ftp.retrbinary("RETR %s" % path, f.write)
                if not resp.startswith("226 "):
                    print("Warning: transfer failed? (response: %s)" % resp)
        
        if checksumsMatches:
            path = "%s%s" % (self.getDirName("gff3"), checksumsMatches[0][1])
            print("Fetching checksums from %s..." % path)
            with open(self.getLocalFilename(checksumsMatches[0][1]) + ".gff3", "wb") as f:
                resp = self._ftp.retrbinary("RETR %s" % path, f.write)
                if not resp.startswith("226 "):
                    print("Warning: transfer failed? (response: %s)" % resp)
                    
        return localGFF3Filename

    def fetchAll(self):
        fn1 = self.fetchGenomeFiles()
        fn2 = self.fetchCDSFiles()
        fn3 = self.fetchGFF3Files()
        return (fn1, fn2, fn3)

    def close(self):
        self._ftp.quit()

def standaloneRun():
    import argparse
    import sys
    
    argsParser = argparse.ArgumentParser()
    #argsParser.add_argument("--verbose", type=int, default=0)
    argsParser.add_argument("--local-name", type=str, required=True)
    argsParser.add_argument("--remote-name", type=str, required=True)
    argsParser.add_argument("--release", type=int, default=37)
    argsParser.add_argument("--section", type=str, default="bacteria")
    argsParser.add_argument("--subsection", type=str)
    args = argsParser.parse_args()


    #ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/fasta/bacteria_3_collection/wolinella_succinogenes_dsm_1740/dna/
    #ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/fasta/bacteria_3_collection/wolinella_succinogenes_dsm_1740
    #ftp://ftp.ensemblgenomes.org/pub/bacteria/release-36/gff3/bacteria_3_collection/wolinella_succinogenes_dsm_1740

    #f = EnsemblFTP("Wsuccinogenes", "wolinella_succinogenes_dsm_1740", release=36, subsection="bacteria_3_collection")
    f = EnsemblFTP(args.local_name, args.remote_name, release=args.release, section=args.section, subsection=args.subsection)
    #f.fetchGenomeFiles()
    #f.fetchCDSFiles()
    #f.fetchGFF3Files()
    print(f.fetchAll())

    sys.exit(0)
    
    
if __name__=="__main__":
    standaloneRun()
