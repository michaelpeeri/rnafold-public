# RNAFold - Analyze mRNA folding bias (Local Folding Energy) based on randomizations.
# Copyright (C) 2016-2020 Michael Peeri
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import os.path
import hashlib
import cPickle as pickle
from collections import Sequence


class ParamHash(object):
    def __init__(self):
        self.m = hashlib.md5()
        
    def add(self, var):

        if isinstance(var, Sequence) and not isinstance(var, type("")):
            for elem in var:
                self.m.update(str(elem))
        else:
            self.m.update(str(var))

    def hexdigest(self):
        return self.m.hexdigest()

"""
tag - string identifier to identify cache files
idvars - dictionary specifying which of the parameters define the caching key
e.g.: @pcache("gene-CAI")
"""
def pcache(tag, idvars={}):
    assert(not(idvars)) # Not impl. currently
    
    def wrapper(wrappedFunc):
        
        def cachedCall(*args, **kw):
            
            # calculate lookup key for this call (based on tag + args + kw)
            m = ParamHash()
            # # Old impl. - allow identity vars to be explicitly specified
            # for varpos, vartype in idvars.items():
            #     value = None

            #     # Get this variable's value, accorind to varpos
            #     if type(varpos)==type(0):
            #         # Positional
            #         value = args[varpos]

            #     elif type(varpos)==type(""):
            #         value = kw[varpos]
            #     else:
            #         raise Exception("Don't know how to find variable '{}'".format(varpos))

            #     # Add this variable to the cache
            #     m.add(value)

            # New impl: All arguments (args+kw) are identity variables
            for arg in args:
                #print("-: [x] {}".format(arg))
                m.add(arg)

            for ident in sorted(kw.keys()): # iterate over the keys in well-defined order
                #print("-: [{}] <val> ".format(ident))
                m.add( kw[ident] )

            lookupKey = "pc_{}_{}.pickle".format( tag, m.hexdigest() )
            print("Key: {}".format(lookupKey))

            retVal = None

            if os.path.isfile( lookupKey ):
                # Cached value found; return it
                print("--Cached")
                with open(lookupKey, "rb") as f:
                    retVal = pickle.load( f )
            else:
                # Cached value not found# call the function
                print("--New")
                retVal = wrappedFunc( *args, **kw )
                # store in cache
                with open(lookupKey, "wb") as f:
                    pickle.dump( retVal, f, pickle.HIGHEST_PROTOCOL )
                
            return retVal
            
        return cachedCall
            
    return wrapper

    
def getRandomNucleotideSeq(length=100):
    import random # test only
    return ''.join([random.choice("acgt") for _ in range(length)])

@pcache("gene-CAI")
def slowFunc( allSeqs, highlyExpressedSeqs, geneticCode ):
    from time import sleep
    print("(doing something slow...)")
    sleep(2)
    a = []
    for _ in range(20):
        a.append( getRandomNucleotideSeq(1000) )  # Note: this is random and violates the deterministic assumption behind the caching. But this is just a test, so it doesn't matter...
        
    return (len(allSeqs), a, geneticCode)

def testPcache():
    slowFunc( ["1","2","3"], ["4","5","6"], 909 )
    slowFunc( ["1","2","3"], ["4","5","6"], 909 )
    return 0
    
if __name__=="__main__":
    import sys
    sys.exit( testPcache() )

