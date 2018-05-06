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



args = argsParser.parse_args()

allSpecies = frozenset(allSpeciesSource())
excludedSpecies = frozenset(args.exclude_species)
candidateSpecies = allSpecies - excludedSpecies

print("All species:\t\t{}".format(len(allSpecies)))
print("Excluded species:\t\t{}".format(len(excludedSpecies)))
print("Candidate species:\t\t{}".format(len(candidateSpecies)))

#print(random.sample( candidateSpecies, args.size ))
species = list(candidateSpecies)
random.shuffle(species)
print(",".join(map(str,species[:args.size])))
