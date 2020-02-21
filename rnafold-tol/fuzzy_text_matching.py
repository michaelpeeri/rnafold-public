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
import swalign # https://github.com/mbreese/swalign
# conda install --channel https://conda.anaconda.org/bioconda swalign


"""
Scoring matrix for Needleman-Wunsch edit distance, for fuzzy matching taxon names
"""
class FuzzySpeciesNameScoringMatrix(object):
    def __init__(self):
        self._wildcard =  0
        self._match    =  5
        self._mismatch = -5
        
    def score( self, one, two, wildcard=None):
        #if wildcard and (one in wildcard and two in wildcard):
        #    return self._wildcard

        if one.lower()==two.lower() and one.isalnum() and two.isalnum():  # Allow matching of characters, prevent treating whitespace as a match
            #print(" '%s' '%s'" % (one, two))
            return self._match

        elif (not one.isalnum()) and (not two.isalnum()):  # Whitespace matches don't count
            return 0
        
        elif one.lower()!=two.lower() and (one.isdigit() or two.isdigit()):   # Give a severe penalty for mismatching digits
            return 10*self._mismatch
        
        else:
            assert(one != two)
            return self._mismatch
    

sw = swalign.LocalAlignment(
        FuzzySpeciesNameScoringMatrix(),
        -9, -1,
    verbose=False, globalalign=False, full_query=True)

def test():
    #NCBI    = "Candidatus Curtissbacteria bacterium GW20 11 GWA1 40 16"
    NCBI2    = "Curtissbacteria bacterium GW20 11 GWA1 40 16"
    Nmicros6 = "Curtissbacteria GWA1 OP11 40 13 partial"
    Nmicros7 = "Curtissbacteria GWA1 OP11 40 16 partial"
    Nmicros8 = "Curtissbacteria GWA1 OP11 40 16"
    Nmicros82 = "Curtissbacteria GWA1 OP11 40 13"

    #sw.align( Nmicros6, NCBI2 ).dump()
    #sw.align( NCBI2, Nmicros6 ).dump()
    #sw.align( Nmicros6, Nmicros6 ).dump()
    #sw.align( NCBI2, Nmicros8 ).dump()
    sw.align( Nmicros8, NCBI2 ).dump()
    #sw.align( Nmicros8, Nmicros8 ).dump()

    #sw.align( NCBI2, Nmicros82 ).dump()
    sw.align( Nmicros82, NCBI2 ).dump()

    sw.align("Woesebacteria GWD2 OP11 40 19", "Candidatus Woesebacteria bacterium GW2011 GWD2 40 19").dump()

    return 0


def standalone():
    import sys
    sys.exit(test())
if __name__=="__main__":
    standalone()
