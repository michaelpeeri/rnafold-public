from __future__ import print_function

test1 = ("Borreliella_bavariensis_pbi.ASM19621v1.37.abinitio.gff3.gz",
         "Borreliella_bavariensis_pbi.ASM19621v1.37.chromosome.Chromosome.gff3.gz",
         "Borreliella_bavariensis_pbi.ASM19621v1.37.chromosome.cp26.gff3.gz",
         "Borreliella_bavariensis_pbi.ASM19621v1.37.chromosome.lp54.gff3.gz",
         "Borreliella_bavariensis_pbi.ASM19621v1.37.gff3.gz",
         "README")


test2 = ("CHECKSUMS",
         "Paramecium_tetraurelia.ASM16542v1.90.abinitio.gff3.gz",
         "Paramecium_tetraurelia.ASM16542v1.90.chr.gff3.gz",
         "Paramecium_tetraurelia.ASM16542v1.90.chromosome.undetermined.gff3.gz",
         "Paramecium_tetraurelia.ASM16542v1.90.gff3.gz",
         "README")

test3 = ("CHECKSUMS",
         "README",
         "Theileria_annulata.ASM322v1.90.abinitio.gff3.gz",
         "Theileria_annulata.ASM322v1.90.chr.gff3.gz",
         "Theileria_annulata.ASM322v1.90.chromosome.2.gff3.gz",
         "Theileria_annulata.ASM322v1.90.chromosome.3.gff3.gz",
         "Theileria_annulata.ASM322v1.90.chromosome.4.gff3.gz",
         "Theileria_annulata.ASM322v1.90.chromosome.Mt.gff3.gz",
         "Theileria_annulata.ASM322v1.90.gff3.gz")


test4 = ("Aspergillus_fumigatus.CADRE.37.abinitio.gff3.gz",
         "Aspergillus_fumigatus.CADRE.37.chromosome.I.gff3.gz",
         "Aspergillus_fumigatus.CADRE.37.chromosome.II.gff3.gz",
         "Aspergillus_fumigatus.CADRE.37.chromosome.III.gff3.gz",
         "Aspergillus_fumigatus.CADRE.37.chromosome.IV.gff3.gz",
         "Aspergillus_fumigatus.CADRE.37.chromosome.MT.gff3.gz",
         "Aspergillus_fumigatus.CADRE.37.chromosome.V.gff3.gz",
         "Aspergillus_fumigatus.CADRE.37.chromosome.VI.gff3.gz",
         "Aspergillus_fumigatus.CADRE.37.chromosome.VII.gff3.gz",
         "Aspergillus_fumigatus.CADRE.37.chromosome.VIII.gff3.gz",
         "Aspergillus_fumigatus.CADRE.37.gff3.gz",
         "CHECKSUMS",
         "README")



test5 = ("CHECKSUMS",
         "README",
         "Ustilago_maydis.Umaydis521_2.0.37.abinitio.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chr.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.1.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.10.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.11.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.12.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.13.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.14.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.15.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.16.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.17.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.18.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.19.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.2.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.20.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.21.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.22.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.23.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.3.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.4.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.5.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.6.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.7.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.8.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.chromosome.9.gff3.gz",
         "Ustilago_maydis.Umaydis521_2.0.37.gff3.gz")


"""
Return the longest common prefix of all strings, or "" if no prefix exists
"""
def findCommonPrefix(strs):
    for i,c in enumerate(strs[0]):
        allMatch = True
        
        for s in strs:
            if s[i]!=c:
                allMatch = False
                break
            
        if not allMatch:
            return strs[0][:i]
        
    return ""

def filterSpecialFiles(options):
    return [x for x in options if x!="README" and x!="CHECKSUMS"]

def selectGff3File(options):
    prefix = findCommonPrefix(filterSpecialFiles(options))
    assert(prefix)
    assert(prefix[-1]==".")

    expectedName = prefix+"gff3.gz"

    assert(expectedName.find("chromosome") == -1)

    if expectedName in options:
        return expectedName
    else:
        return None


print(selectGff3File(test1))
print(selectGff3File(test2))
print(selectGff3File(test3))
print(selectGff3File(test4))
print(selectGff3File(test5))

