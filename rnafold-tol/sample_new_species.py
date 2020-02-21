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
from data_helpers import allSpeciesSource
import random
import argparse

def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert

argsParser = argparse.ArgumentParser()
argsParser.add_argument("--exclude-species", type=parseList(int), default=())
argsParser.add_argument("--size", type=int, default=10)
argsParser.add_argument("--sample-from-taxon", type=int)



args = argsParser.parse_args()

allSpecies = frozenset(allSpeciesSource())
excludedSpecies = frozenset(args.exclude_species)
candidateSpecies = allSpecies - excludedSpecies

print("All species:\t\t{}".format(len(allSpecies)))
print("Excluded species:\t\t{}".format(len(excludedSpecies)))
print("Candidate species:\t\t{}".format(len(candidateSpecies)))



species = list(candidateSpecies)
random.shuffle(species)

ret = []

if args.sample_from_taxon is None: # unrestricted sampling
    #print(random.sample( candidateSpecies, args.size ))
    ret = species[:args.size]

else: # sample from species restricted by taxon
    from ncbi_taxa import ncbiTaxa

    ret = []

    for taxid in species:
        lineage = ncbiTaxa.get_lineage( taxid )
        if args.sample_from_taxon in lineage:
            ret.append(taxid)

        if len(ret) >= args.size: break

print(",".join(map(str,ret)))
