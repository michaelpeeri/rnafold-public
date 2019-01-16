import argparse


def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert

profilesSpecs = []
def parseProfileSpec():
    def convert(value):
        o = value.split(':')
        profilesSpecs.append(o)
        return o
    return convert

argsParser = argparse.ArgumentParser()
argsParser.add_argument("--taxid", type=parseList(int))
argsParser.add_argument("--profile", type=parseProfileSpec())
#argsParser.add_argument("--variant", type=parseOption(set(("yeastgenome", "NCBI")), "variant"))
#argsParser.add_argument("--type", type=parseOption(set(("cds", "shuffle", "fixCDSkey")), "sequence type"))
#argsParser.add_argument("--dry-run", action="store_true", default=False)
#argsParser.add_argument("--output-fasta")
#argsParser.add_argument("--gene-ids-file")
#argsParser.add_argument("--alt-protein-ids", type=parseOption(set(("locus_tag",)), "alt-protein-id"))
#argsParser.add_argument("--headers-from-another-fasta")
args = argsParser.parse_args()


#print(args.taxid)
print(profilesSpecs)
