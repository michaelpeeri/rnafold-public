#from BCBio.GFF import GFFExaminer # pip install bcbio-gff
# See http://biopython.org/wiki/GFF_Parsing
# See https://github.com/chapmanb/bcbb/tree/master/gff
import gffutils



def createGffDb( filename ):
    db = gffutils.create_db(filename, dbfn="%s.db" % filename, force=True, keep_order=True,
                            sort_attribute_values=True,
                            merge_strategy='merge',
                            id_spec=['ID', 'Name'])
    return db



