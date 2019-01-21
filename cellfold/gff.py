# Access GFF data sources
import gffutils  # See: http://daler.github.io/gffutils/contents.html
import os

# Get a GFF FeatureDB object for a given .gff file; Implement policies for this
def createGffDb( filename, variant ):

    if os.path.exists( filename ):
        db = gffutils.FeatureDB("{}.db".format(filename), keep_order=True )
        
    else:
        db = gffutils.create_db(filename, dbfn="{}.db".format(filename), force=True, keep_order=True,
                                sort_attribute_values=True,
                                merge_strategy='merge',
                                id_spec=['ID', 'Name'])
    return db



