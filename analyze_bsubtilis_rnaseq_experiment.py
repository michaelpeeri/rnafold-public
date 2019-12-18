from csv import reader
from tempfile import NamedTemporaryFile
import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, mannwhitneyu
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
import seaborn as sns
from data_helpers import CDSHelper, SpeciesCDSSource, getSpeciesName, getSpeciesTranslationTable
from process_series_data import readSeriesResultsForSpecies, convertResultsToMFEProfiles, sampleProfilesFixedIntervals, profileElements #, profileLength, MeanProfile, calcSampledGCcontent, profileEdgeIndex, readSeriesResultsForSpeciesWithSequence
from codonw import readCodonw, calcCAI
import gff
from mysql_rnafold import Sources
from pcache import pcache
from paxdb import getSpeciesPaxdbData
from create_identifiers_mapping import createMappingForSpeciesProteins

# Configuration
speciesConfig = {224308: dict( rs = ('./data/sra/SRR3466199.trimmed.vs.Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.dna_rm.toplevel.fa.bam.union.htseq',
                                     './data/sra/SRR3466200.trimmed.vs.Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.dna_rm.toplevel.fa.bam.union.htseq',
                                     './data/sra/SRR3466201.trimmed.vs.Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.dna_rm.toplevel.fa.bam.union.htseq'),
                               refvalsfn   = './data/sra/GSE80786_counts_matrix.txt',
                               gff3fn      = './data/sra/Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.33.chromosome.Chromosome.gff3.gz',
                               gff3variant = "Ensembl"),
                 
                 511145: dict(
                     rs = (),
                     gff3fn      = './data/sra/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 208964: dict(  #SRR10259076.trimmed.vs.GCF_000006765.1_ASM676v1_genomic.fna.bam.union.nostrand.htseq
                     rs = ('./data/sra/SRR10259076.trimmed.vs.GCF_000006765.1_ASM676v1_genomic.fna.bam.union.nostrand.htseq',
                           './data/sra/SRR10259076.trimmed.vs.GCF_000006765.1_ASM676v1_genomic.fna.bam.union.nostrand.htseq',
                           './data/sra/SRR10259076.trimmed.vs.GCF_000006765.1_ASM676v1_genomic.fna.bam.union.nostrand.htseq'),
                     gff3fn      = './data/sra/GCF_000006765.1_ASM676v1_genomic.gff.gz',
                     gff3variant = "NCBI"
                     ),
                 100226: dict(
                     rs = ('./data/sra/SRR6782725.trimmed.vs.Streptomyces_coelicolor_a3_2_.ASM20383v1.dna_rm.toplevel.fa.byname.bam.union.htseq',
                           './data/sra/SRR6782725.trimmed.vs.Streptomyces_coelicolor_a3_2_.ASM20383v1.dna_rm.toplevel.fa.byname.bam.union.htseq',
                           './data/sra/SRR6782725.trimmed.vs.Streptomyces_coelicolor_a3_2_.ASM20383v1.dna_rm.toplevel.fa.byname.bam.union.htseq'),
                     gff3fn      = './data/sra/Streptomyces_coelicolor_a3_2_.ASM20383v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 220668: dict(
                     rs = ('./data/sra/SRR7276154.trimmed.vs.Lactobacillus_plantarum_wcfs1.ASM20385v3.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR8349567.trimmed.vs.Lactobacillus_plantarum_wcfs1.ASM20385v3.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR8349567.trimmed.vs.Lactobacillus_plantarum_wcfs1.ASM20385v3.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Lactobacillus_plantarum_wcfs1.ASM20385v3.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 269084: dict(
                     rs = ('./data/sra/SRR10421010.trimmed.vs.Synechococcus_elongatus_pcc_6301.ASM1006v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR10421010.trimmed.vs.Synechococcus_elongatus_pcc_6301.ASM1006v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR10421010.trimmed.vs.Synechococcus_elongatus_pcc_6301.ASM1006v1.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Synechococcus_elongatus_pcc_6301.ASM1006v1.32.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 64091: dict(
                     rs = ('./data/sra/SRR5651597.trimmed.vs.Halobacterium_salinarum_nrc_1.ASM680v1.dna_rm.toplevel.fa.bam.union.htseq',
                           './data/sra/SRR5651597.trimmed.vs.Halobacterium_salinarum_nrc_1.ASM680v1.dna_rm.toplevel.fa.bam.union.htseq',
                           './data/sra/SRR5651597.trimmed.vs.Halobacterium_salinarum_nrc_1.ASM680v1.dna_rm.toplevel.fa.bam.union.htseq'),
                     gff3fn      = './data/sra/Halobacterium_salinarum_nrc_1.ASM680v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 882: dict(
                     rs = ('./data/sra/SRR5871740.trimmed.vs.Desulfovibrio_vulgaris_str_hildenborough.ASM19575v1.dna_rm.toplevel.union.htseq',
                           './data/sra/SRR5871740.trimmed.vs.Desulfovibrio_vulgaris_str_hildenborough.ASM19575v1.dna_rm.toplevel.union.htseq',
                           './data/sra/SRR5871740.trimmed.vs.Desulfovibrio_vulgaris_str_hildenborough.ASM19575v1.dna_rm.toplevel.union.htseq'),
                     gff3fn      = './data/sra/Desulfovibrio_vulgaris_str_hildenborough.ASM19575v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 211586: dict(
                     rs = ('./data/sra/SRR6170086.trimmed.vs.Shewanella_oneidensis_mr_1.ASM14616v2.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR6170086.trimmed.vs.Shewanella_oneidensis_mr_1.ASM14616v2.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR6170086.trimmed.vs.Shewanella_oneidensis_mr_1.ASM14616v2.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Shewanella_oneidensis_mr_1.ASM14616v2.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 226186: dict(
                     rs = ('./data/sra/SRR8874378.trimmed.vs.Bacteroides_thetaiotaomicron_vpi_5482.ASM1106v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR8874378.trimmed.vs.Bacteroides_thetaiotaomicron_vpi_5482.ASM1106v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR8874378.trimmed.vs.Bacteroides_thetaiotaomicron_vpi_5482.ASM1106v1.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Bacteroides_thetaiotaomicron_vpi_5482.ASM1106v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 190304: dict(
                     rs = ('./data/sra/SRR9010011.trimmed.vs.Fusobacterium_nucleatum_subsp_nucleatum_atcc_25586.ASM732v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR9010011.trimmed.vs.Fusobacterium_nucleatum_subsp_nucleatum_atcc_25586.ASM732v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR9010011.trimmed.vs.Fusobacterium_nucleatum_subsp_nucleatum_atcc_25586.ASM732v1.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Fusobacterium_nucleatum_subsp_nucleatum_atcc_25586.ASM732v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 264462: dict(
                     rs = ('./data/sra/SRR3605968.trimmed.vs.Bdellovibrio_bacteriovorus_hd100.ASM19617v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR3605968.trimmed.vs.Bdellovibrio_bacteriovorus_hd100.ASM19617v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR3605968.trimmed.vs.Bdellovibrio_bacteriovorus_hd100.ASM19617v1.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Bdellovibrio_bacteriovorus_hd100.ASM19617v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 93061: dict(
                     rs = ('./data/sra/SRR10355925.trimmed.vs.Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.chromosome.Chromosome.bam.union.htseq',
                           './data/sra/SRR10355925.trimmed.vs.Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.chromosome.Chromosome.bam.union.htseq',
                           './data/sra/SRR10355925.trimmed.vs.Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.chromosome.Chromosome.bam.union.htseq'),
                     gff3fn      = './data/sra/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.32.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 169963: dict(
                     rs = ('./data/sra/SRR9167346.trimmed.vs.Listeria_monocytogenes_egd_e.ASM19603v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR9167346.trimmed.vs.Listeria_monocytogenes_egd_e.ASM19603v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR9167346.trimmed.vs.Listeria_monocytogenes_egd_e.ASM19603v1.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Listeria_monocytogenes_egd_e.ASM19603v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 486041: dict(
                     rs = ('./data/sra/SRR6709977.trimmed.vs.Laccaria_bicolor_s238n_h82.V1.0.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR6709977.trimmed.vs.Laccaria_bicolor_s238n_h82.V1.0.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR6709977.trimmed.vs.Laccaria_bicolor_s238n_h82.V1.0.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Laccaria_bicolor_s238n_h82.V1.0.37.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 6669: dict(
                     rs = ('./data/sra/SRR3356112.trimmed.vs.Daphnia_pulex.GCA_000187875.1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR3356112.trimmed.vs.Daphnia_pulex.GCA_000187875.1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR3356112.trimmed.vs.Daphnia_pulex.GCA_000187875.1.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Daphnia_pulex.GCA_000187875.1.33.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 400682: dict(
                     rs = ('./data/sra/SRR8951896.trimmed.vs.Amphimedon_queenslandica.Aqu1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR8951896.trimmed.vs.Amphimedon_queenslandica.Aqu1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR8951896.trimmed.vs.Amphimedon_queenslandica.Aqu1.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Amphimedon_queenslandica.Aqu1.33.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 45351: dict(
                     rs = ('./data/sra/SRR5839397.trimmed.vs.Nematostella_vectensis.ASM20922v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR5839397.trimmed.vs.Nematostella_vectensis.ASM20922v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR5839397.trimmed.vs.Nematostella_vectensis.ASM20922v1.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Nematostella_vectensis.ASM20922v1.37.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 1397361: dict(
                     rs = ('./data/sra/SRR9602168.trimmed.vs.Sporothrix_schenckii_1099_18.S_schenckii_v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR9602168.trimmed.vs.Sporothrix_schenckii_1099_18.S_schenckii_v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR9602168.trimmed.vs.Sporothrix_schenckii_1099_18.S_schenckii_v1.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Sporothrix_schenckii_1099_18.S_schenckii_v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 5061: dict(
                     rs = ('./data/sra/SRR7772837.trimmed.vs.Aspergillus_niger.ASM285v2.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR7772837.trimmed.vs.Aspergillus_niger.ASM285v2.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR7772837.trimmed.vs.Aspergillus_niger.ASM285v2.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Aspergillus_niger.ASM285v2.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 176299: dict(
                     rs = ('./data/sra/SRR6227211.trimmed.vs.Agrobacterium_fabrum_str_c58.ASM9202v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR6227211.trimmed.vs.Agrobacterium_fabrum_str_c58.ASM9202v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR6227211.trimmed.vs.Agrobacterium_fabrum_str_c58.ASM9202v1.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Agrobacterium_fabrum_str_c58.ASM9202v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 314225: dict(
                     rs = ('./data/sra/SRR8571520.trimmed.vs.Erythrobacter_litoralis_htcc2594.ASM1300v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR8571520.trimmed.vs.Erythrobacter_litoralis_htcc2594.ASM1300v1.dna_rm.toplevel.bam.union.htseq',
                           './data/sra/SRR8571520.trimmed.vs.Erythrobacter_litoralis_htcc2594.ASM1300v1.dna_rm.toplevel.bam.union.htseq'),
                     gff3fn      = './data/sra/Erythrobacter_litoralis_htcc2594.ASM1300v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                     ),
                 1148: dict( #Synechocystis sp. PCC 6803
                     rs = (),
                     gff3fn      = './data/sra/Synechocystis_sp_pcc_6803.ASM972v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 158878: dict( #Staphylococcus aureus subsp. aureus Mu50
                     rs = (),
                     gff3fn      = './data/sra/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.32.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 160490: dict( #Streptococcus pyogenes M1 GAS
                     rs = (),
                     gff3fn      = './data/sra/Streptococcus_pyogenes_m1_gas.ASM678v2.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 192222: dict( #Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819
                     rs = (),
                     gff3fn      = './data/sra/Campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819.ASM908v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 243159: dict( #Acidithiobacillus ferrooxidans ATCC 23270
                     rs = (),
                     gff3fn      = './data/sra/Acidithiobacillus_ferrooxidans_atcc_23270.ASM2148v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 267671: dict( #Leptospira interrogans serovar Copenhageni str. Fiocruz L1-130
                     rs = (),
                     gff3fn      = './data/sra/Leptospira_interrogans_serovar_copenhageni_str_fiocruz_l1_130.ASM768v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 272623: dict( #Lactococcus lactis subsp. lactis Il1403
                     rs = (),
                     gff3fn      = './data/sra/Lactococcus_lactis_subsp_lactis_il1403.ASM686v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 283166: dict( #Bartonella henselae str. Houston-1
                     rs = (),
                     gff3fn      = './data/sra/Bartonella_henselae_str_houston_1.ASM4670v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 449447: dict( #Microcystis aeruginosa NIES-843
                     rs = (),
                     gff3fn      = './data/sra/Microcystis_aeruginosa_nies_843.ASM1062v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 # 546414: dict( #Deinococcus deserti VCD115
                 #     rs = (),
                 #     gff3fn      = './data/sra/', # ???????
                 #     gff3variant = "Ensembl"
                 # ),
                 593117: dict( #Thermococcus gammatolerans EJ3
                     rs = (),
                     gff3fn      = './data/sra/Thermococcus_gammatolerans_ej3.ASM2236v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 722438: dict( #Mycoplasma pneumoniae FH
                     rs = (),
                     gff3fn      = './data/sra/Mycoplasma_pneumoniae_m129.ASM2734v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 83332: dict( #Mycobacterium tuberculosis H37Rv
                     rs = (),
                     gff3fn      = './data/sra/Mycobacterium_tuberculosis_h37rv.ASM19595v2.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 85962: dict( #Helicobacter pylori 26695
                     rs = (),
                     gff3fn      = './data/sra/Helicobacter_pylori_26695_gca_000307795.ASM30779v1.36.gff3.gz',
                     gff3variant = "Ensembl"
                 ),
                 99287: dict( #Salmonella enterica subsp. enterica serovar Typhimurium str. LT2
                     rs = (),
                     gff3fn      = './data/sra/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.36.gff3.gz',
                     gff3variant = "Ensembl"
                 )                 
                 #
                 #
                 }
config = None

#taxId       = 224308 # = B. subtilis
#taxId       = 511145 # = B. subtilis
# Experimental data for WT (control) reps, processed by me
#rs          = ('./data/sra/SRR3466199.trimmed.vs.Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.dna_rm.toplevel.fa.bam.union.htseq',
#               './data/sra/SRR3466200.trimmed.vs.Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.dna_rm.toplevel.fa.bam.union.htseq',
#               './data/sra/SRR3466201.trimmed.vs.Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.dna_rm.toplevel.fa.bam.union.htseq')
# Reference processed values (analyzed by the original authors)
#refvalsfn   = './data/sra/GSE80786_counts_matrix.txt'
# GFF3 for genome
#gff3fn      = './data/sra/Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.33.chromosome.Chromosome.gff3.gz'
#gff3fn      = './data/sra/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz'
#gff3variant = "Ensembl"
# Misc
debugMode = False
# DLFE data
# mfe_v2_40nt_num_genes_vs_position_310_10_begin_t11_Bsubtilis.pdf
# mfe_v2_40nt_num_genes_vs_position_310_10_end_t11_Bsubtilis.pdf
#dLFEComputationTag = Sources.RNAfoldEnergy_SlidingWindow40_v2
#dLFEShuffleType    = Sources.ShuffleCDSv2_python 
numShuffledGroups  = 20

def readMyResults(fn):
    ret = {}
    with open(fn) as csvfile:
        for row in reader( csvfile, delimiter='\t' ):
            ret[ row[0] ] = int(row[1])
    return ret

def readRefResults(fn):
    ret = {}
    with open(fn) as csvfile:
        for row in reader( csvfile, delimiter=' ' ):
            if row[0]=='WT_rep1': continue
            ret[ row[0] ] = ( int(row[1]), int(row[2]), int(row[3]) )  # get data for WT (control) reps
    return ret

# def parseList(conversion=str):
#     def convert(values):
#         return map(conversion, values.split(","))
#     return convert

def parseProfileSpec():
    def convert(value):
        o = value.split(':')
        assert(len(o) >= 3 and len(o) <= 4)
        
        o[0] = int(o[0])
        assert(o[0]>0)
        
        o[1] = int(o[1])
        assert(o[1]>0)
        
        assert(o[2]=="begin" or o[2]=="end")

        if( len(o) == 4 ):
            o[3] = int(o[3])
        else:
            o.append(0)
        
        return (o[0], o[1], o[2], o[3])
    return convert


def getRegions(db, variant):
    out = []
    if( variant=="NCBI"):
        for region in db.features_of_type("region"):
            out.append(region[0])
            
    elif( variant=="Ensembl"):
        out=[None]
        
    else:
        assert(False)
    return out

def CDS_and_mRNA_source():
    print("Loading gff3 file...")
    db = gff.createGffDb(config['gff3fn'], config['gff3variant'])
    print("Done.")

    for area in getRegions( db, config['gff3variant'] ):
        if( config['gff3variant']=="Ensembl"):
            for cds in db.features_of_type('CDS'):
                #print("CDS.id: %s" % cds.id)
                mRNA = None
                gene = None
                for p in db.parents(cds.id):
                    #print("mRNA.featuretype: %s" % (mRNA.featuretype))
                    if p.featuretype   == 'gene':
                        gene = p
                    elif p.featuretype in ('mRNA', 'miRNA'): # miRNA actually applies to protein-coding genes...
                        mRNA = p
                    elif p.featuretype == 'transcript':
                        mRNA = p

                yield cds, mRNA, gene
                
        elif( config['gff3variant']=="NCBI"):
            for cds in db.features_of_type('CDS', limit=(area, 0, 100000000), completely_within=False):
                #for mRNA in db.parents(cds.id):
                for gene in db.parents(cds.id):

                    yield cds, None, gene

gene2prot = None
prot2gene = None
rePutativeHighlyExpressedGeneFunctions = (re.compile(".*ribosom.*"),
                                          re.compile(".*chaperon.*"),
                                          re.compile(".*heat[ -]shock.*"),
                                          re.compile(".*elongation.*"),
                                          re.compile(".*dna[ -]binding.*"),
                                          re.compile(".*trna.*"),
                                          re.compile(".*rna[ -]polymerase.*")          )
rePutativeLowlyExpressedGeneFunctions = (re.compile(".*pseudo.*"),
                                         re.compile(".*transcription[- ]factor.*")  )
putativeHighlyExpressedGenes = None
putativeLowlyExpressedGenes  = None
geneDescriptions = None

def dedup(xs):
    ret = []
    e = set()

    for x in xs:
        if not x in e:
            ret.append(x)
            e.add(x)
            
    return ret

def resetMappings():
    global gene2prot, prot2gene, putativeHighlyExpressedGenes, putativeLowlyExpressedGenes, geneDescriptions
    gene2prot = {}
    prot2gene = {}
    putativeHighlyExpressedGenes = set()
    putativeLowlyExpressedGenes  = set()
    geneDescriptions = {}

def getGeneNameMapping():
    global gene2prot, prot2gene, putativeHighlyExpressedGenes, putativeLowlyExpressedGenes, geneDescriptions
    
    if gene2prot:
        return gene2prot, prot2gene

    # Load the data
    for cds, mRNA, gene in CDS_and_mRNA_source( ):
        if mRNA:
            biotype = mRNA.attributes['biotype'][0]
        else:
            biotype = gene.attributes['gene_biotype'][0]
            
        if( biotype != 'protein_coding' ):
            continue
        #print(cds.attributes)
        protId = cds.attributes['protein_id'][0]

        #geneId = mRNA.attributes['transcript_id'][0]
        if mRNA and (not 'Name' in mRNA.attributes) and (not 'gene_id' in gene.attributes):
            continue

        geneIds = []
        if mRNA and 'Name' in mRNA.attributes:
            geneIds.append( mRNA.attributes['Name'][0] )
            
        if 'Name' in gene.attributes:
            geneIds.append( gene.attributes['Name'][0] )

        if 'gene_id' in gene.attributes:
            geneIds.append( gene.attributes['gene_id'][0] )
        assert(geneIds) # at least one gene-id is required

        for i in range(len(geneIds)):
            if len(geneIds[i])>3 and geneIds[i][-2]=='-': # remove trailing '-1'
                geneIds[i] = geneIds[i][:-2]

        desc = None
        if 'description' in gene.attributes:
            desc = gene.attributes['description'][0].lower()
        elif 'product' in cds.attributes:
            desc = cds.attributes['product'][0].lower()
        else:
            print("Warning: {}".format(cds.attributes))

        if desc:
            matches = [r.match(desc) for r in rePutativeHighlyExpressedGeneFunctions]
            if any(matches):
                print("*** {}".format(desc))
                putativeHighlyExpressedGenes.add( geneIds[0] )
            else:
                print("    {}".format(desc))

            matches = [r.match(desc) for r in rePutativeLowlyExpressedGeneFunctions]
            if any(matches):
                print("xxx {}".format(desc))
                putativeLowlyExpressedGenes.add( geneIds[0] )

            geneDescriptions[geneIds[0]] = desc # store the description
            

        prot2gene[protId] = tuple(dedup(geneIds))
        gene2prot[geneIds[0]] = protId
        
    return gene2prot, prot2gene

def calcCodonwMeasures( seqs ):
    assert( bool(seqs) )
    print("Calculating codonw measures...")
    fout = NamedTemporaryFile( mode="w", delete=(not debugMode) )
    SeqIO.write( seqs, fout.name, "fasta")  # write the full sequences into the file
    return readCodonw( fout.name )          # run codonw and get the gene-level results
    #u'T3s', u'C3s', u'A3s', u'G3s', u'CAI', u'CBI', u'Fop', u'Nc', u'GC3s',  u'GC', u'L_sym', u'L_aa', u'Gravy', u'Aromo', u'Unnamed: 15'


@pcache("gene-CAI")
def calcCAIusingReferenceGenes( allSeqs, highlyExpressedSeqs, geneticCode ):
    print("Calculating CAI for {} sequences (ref={}; code={})".format( len(allSeqs), len(highlyExpressedSeqs), geneticCode))
    fall   = NamedTemporaryFile( mode="w", delete=(False) ) #not debugMode) )
    SeqIO.write( allSeqs,             fall.name,   "fasta")  # write the full sequences into the file
    
    fhiexp = NamedTemporaryFile( mode="w", delete=(False) ) #not debugMode) )
    SeqIO.write( highlyExpressedSeqs, fhiexp.name, "fasta")  # write the highly-expressed sequences into the file

    return calcCAI( fall.name, fhiexp.name, geneticCode=geneticCode )          # call CAI and get the gene-level results


@pcache("series-results")
def readSeriesResultsForSpecies_cached( seriesSourceNumber, species, minShuffledGroups, maxShuffledGroups, shuffleType, cdsFilter=None, returnCDS=True ):
    return list(readSeriesResultsForSpecies(seriesSourceNumber, species, minShuffledGroups, maxShuffledGroups, shuffleType, cdsFilter=cdsFilter, returnCDS=returnCDS ))
    


def getGeneDLFEProfiles(taxId, args, protIds):

    geneDLFEs = pd.DataFrame( columns=profileElements(args.profile), index=protIds, dtype=float )
    geneNames2ProtName, protName2GeneName = getGeneNameMapping()

    # ------------------------------------
    # Process all CDS for this species
    # ------------------------------------
    for result in sampleProfilesFixedIntervals(
            convertResultsToMFEProfiles(
                readSeriesResultsForSpecies_cached( (args.computation_tag,), taxId, args.num_shuffles, args.num_shuffles, shuffleType=args.shuffle_type )
                , args.num_shuffles)
            , args.profile[3], args.profile[0], args.profile[1], args.profile[2]):

        #fullCDS = result["cds-seq"]
        #seq = fullCDS[args.profile[3]:args.profile[0]]

        #if not seq:
        #    continue

        protId = result["cds"].getProtId()
        geneId = protName2GeneName.get( protId, None )[0]
        if geneId is None: continue

        profileData = result["profile-data"]
        assert(profileData.shape[0] >= args.num_shuffles )

        #if profileData.shape[1] < len(profileElements(args.profile)):
        #    print("Skipping {}".format( protId ))
        #    continue

        # Get dLFE profile for this gene
        nativeLFE   = profileData[0,None].flatten()
        shuffledLFE = profileData[1:].mean( axis=0 )
        dLFE        = nativeLFE - shuffledLFE

        geneDLFEs.loc[ geneId, tuple(range(0,profileData.shape[1]*10,10)) ] = dLFE.astype(float)
        
    return geneDLFEs


def getAllCDSSeqs(taxId, protIds, args):
    geneNames2ProtName, protName2GeneName = getGeneNameMapping()

    print("Getting CDS seqs...")
    CDSseqs = []
    for protId in protIds: 
        #protId = geneNames2ProtName.get(geneId, None)
        #if protId is None: continue
        geneIds = protName2GeneName[protId]

        cds = CDSHelper(taxId, protId)

        if not cds.exists():
            continue

        if( cds.length()%3 != 0 ):
            continue

        seq = cds.sequence()

        record = SeqRecord( Seq(seq, NucleotideAlphabet), id=protId, description=" ".join(geneIds) )
        CDSseqs.append( record )
    return CDSseqs


def safeCorr( v1, v2 ):
    notnans = np.logical_not( np.logical_or( np.isnan(v1), np.isnan(v2) ) )
    return( pearsonr( v1[notnans], v2[notnans] ) )

def annotatePvaluesHighVsLow( df, args, var, highLowPercentile ):
    #assert( ranges[1] >= ranges[0] )
    assert( highLowPercentile > 1e-3 and highLowPercentile <= 0.5)
    # calculate percentiles

    loThres, hiThres = df.loc[:, var].quantile( (highLowPercentile, 1-highLowPercentile) )
    dfLo = df[ df.loc[:, var] <= loThres ]
    dfHi = df[ df.loc[:, var] >  hiThres ]
    
    pvalsLo = []
    pvalsHi = []
    for pos in profileElements( args.profile ):
        xs = dfHi.loc[:,pos].values
        ys = dfLo.loc[:,pos].values
        pvalsLo.append( mannwhitneyu( xs, ys, alternative='less'    ).pvalue )
        pvalsHi.append( mannwhitneyu( xs, ys, alternative='greater' ).pvalue )
    
    return (pvalsLo, pvalsHi)


def plotStuff( var, xpos, df, pvalsLess, pvalsGreater, args, highLowPercentile, taxId ):
    fig, ax1 = plt.subplots()

    #nonzero = df[ df.loc[:,'Exp_rep1'] > 0 ]

    
    pcs = df.loc[:, var].quantile( (highLowPercentile, 1-highLowPercentile) )
    ranges = ( pcs.iloc[0], pcs.iloc[1] )
    assert(ranges[1] >= ranges[0])


    def plotGroup(group, color="blue", annotateR=False):
        xs = group.loc[:,xpos].values
        ys = group.loc[:, var].values.argsort().argsort() # convert values to ranks
        plt.scatter(  xs, ys, alpha=0.2, color=color )
        if annotateR:
            r = safeCorr( xs, ys )
            ax1.annotate( s="r={:.3}".format(r[0]), xy=(5, 3500) )

    nonSpecial = df[ df.loc[:,'Putative'] == False ]
    special    = df[ df.loc[:,'Putative'] == True  ]
    plotGroup(nonSpecial, color="blue" )
    plotGroup(special,    color="red" )
    

    #ax1.set_yscale( "log" )
    #plt.ylim( (0, 1e6) )
    ax1.axvline( x=0, c="black" )

    plt.savefig("sra_{}_scatter_{}.pdf".format(taxId, var))
    plt.savefig("sra_{}_scatter_{}.svg".format(taxId, var))
    plt.close(fig)

    
    fig, (ax1,ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [2, 1]})

    dfHi = df[ df.loc[:,var] >  ranges[1] ]
    print(dfHi.shape)
    dfLo = df[ df.loc[:,var].between( 0, ranges[0]) ]
    print(dfLo.shape)

    xvals = profileElements( args.profile )
    ax1.plot( xvals, df.loc[  :,xvals].mean(axis=0), c='#5c70b8',  linewidth=4, alpha=0.7, label="All genes (N={})".format( df.shape[0] ))
    ax1.plot( xvals, dfHi.loc[:,xvals].mean(axis=0), c='#344277',  linewidth=4, alpha=0.7, label="High {} (N={})".format(var, dfHi.shape[0]) )
    ax1.plot( xvals, dfLo.loc[:,xvals].mean(axis=0), c='#b9c2e1',  linewidth=4, alpha=0.7, label="Low {} (N={})".format(var, dfLo.shape[0])  )

    logpval1 = np.log10(pvalsLess)
    logpval2 = np.log10(pvalsGreater)
    ax2.plot( xvals, logpval1, c="#8d6f3e", linewidth=4, alpha=0.7, label="Less"    )
    ax2.plot( xvals, logpval2, c="#3e8e8e", linewidth=4, alpha=0.7, label="Greater" )

    actualPvalRange = (max( np.max(logpval1), np.max(logpval2) ),
                       min( np.min(logpval1), np.min(logpval2) ) )
    minYforLogPvals = min( np.max(logpval1)*1.5, -5.0 )

    ax1.axhline( y=0, c="black" )
    ax1.set_ylim( (-1.5, 0.5) )
    ax1.set_ylabel(u"\u0394LFE")
    ax1.legend(fontsize=12, loc=(0,0), ncol=2)
    ax2.axhline( y=-1.3, c="black" )
    ax2.set_ylim( (minYforLogPvals, 0.0) )
    ax2.set_ylabel('log10(p-value)')
    ax2.legend(fontsize=12, loc=(0,0), ncol=3)
    plt.title( getSpeciesName(taxId) )
    
    plt.savefig("sra_{}_dlfe_high_low_{}_{}.pdf".format(taxId, var, highLowPercentile))
    plt.savefig("sra_{}_dlfe_high_low_{}_{}.svg".format(taxId, var, highLowPercentile))
    plt.close(fig)


    # -------------------------------

    fig, ax1 = plt.subplots()
    p = list(range(100,310,10))
    dlfes1 = df[ df.loc[:,var] >  ranges[1]  ].loc[:, p ].values
    dlfes1 = dlfes1[ np.logical_not( np.isnan( dlfes1 ) ) ]
    dlfes2 = df[ df.loc[:,var] <=  ranges[0]  ].loc[:, p ].values
    dlfes2 = dlfes2[ np.logical_not( np.isnan( dlfes2 ) ) ]
    data = ( dlfes1, dlfes2 )
    print("////////////")
    print(len(data[0]), len(data[1]))
    plt.boxplot( data, labels=("> {}".format(ranges[1]),"<= {}".format(ranges[0])) )
    
    #ax1.set_yscale( "log" )
    #plt.ylim( (0, 1e6) )
    plt.savefig("sra_{}_reads_vs_dlfe{}_{}.pdf".format(taxId, xpos, var))
    plt.savefig("sra_{}_reads_vs_dlfe{}_{}.svg".format(taxId, xpos, var))
    plt.close(fig)

    # -------------------------------

    # fig, ax1 = plt.subplots()
    # data = ( df[ df.loc[:,xpos] >  -0.10  ][var],
    #          df[ df.loc[:,xpos] <= -0.40  ][var] )
    # print("////////////")
    # print(len(data[0]), len(data[1]))
    # plt.boxplot( data, labels=("> -0.10","<=-0.40") )
    
    # ax1.set_yscale( "log" )
    # #plt.ylim( (0, 1e6) )
    # plt.savefig("sra_dlfe{}_vs_reads_{}.pdf".format(xpos, var))
    # plt.savefig("sra_dlfe{}_vs_reads_{}.svg".format(xpos, var))
    # plt.close(fig)
    
    # -------------------------------

    # fig, ax1 = plt.subplots()
    # data = ( df[                 df.loc[:,'Putative']  ][var ],
    #          df[ np.logical_not( df.loc[:,'Putative']) ][var ] )
    # print("////////////")
    # print(len(data[0]), len(data[1]))
    # plt.boxplot( data, labels=("Yes","No") )
    
    # ax1.set_yscale( "log" )
    # #plt.ylim( (0, 1e6) )
    # plt.savefig("sra_reads_{}.pdf".format(var))
    # plt.savefig("sra_reads_{}.svg".format(var))
    # plt.close(fig)

def plotLengthDiagnostics( df, taxId ):
    fig, ax1 = plt.subplots()
    #df.plot.scatter( "Length", "Exp_rep1_norm", ax=ax1 )
    sns.jointplot( "Length", "Exp_rep1_norm", data=df, ax=ax1, kind="reg" )
    plt.title( getSpeciesName(taxId) )
    plt.savefig("sra_{}_length1.pdf".format(taxId))
    plt.savefig("sra_{}_length1.svg".format(taxId))
    plt.close(fig)

    fig, ax1 = plt.subplots()
    #df.plot.scatter( "Exp_rep1", "Exp_rep1_norm", ax=ax1 )
    sns.jointplot( "Exp_rep1", "Exp_rep1_norm", data=df, ax=ax1, kind="reg" )
    plt.title( getSpeciesName(taxId) )
    plt.savefig("sra_{}_length2.pdf".format(taxId))
    plt.savefig("sra_{}_length2.svg".format(taxId))
    plt.close(fig)
    

def processGenome(args, taxId):

    geneticCode = getSpeciesTranslationTable(taxId)
    resetMappings()

    exps = [readMyResults(r) for r in config['rs']]
    if 'refvalsfn' in config:
        ref  = readRefResults(config['refvalsfn'])
    else:
        ref = {}
    geneNames2ProtName, protName2GeneName = getGeneNameMapping()

    allGenes = geneNames2ProtName.keys()
        

    # Populate combined df with all values
    df = pd.DataFrame( {'Exp_rep1': pd.Series( dtype='float' ),
                        'Exp_rep2': pd.Series( dtype='float' ),
                        'Exp_rep3': pd.Series( dtype='float' ),
                        'Ref_rep1': pd.Series( dtype='float' ),
                        'Ref_rep2': pd.Series( dtype='float' ),
                        'Ref_rep3': pd.Series( dtype='float' ),
                        'ProtId':   pd.Series( dtype='str'   ),
                        'Nc':       pd.Series( dtype='float' ),
                        'CAI':      pd.Series( dtype='float' ),
                        'CAI0':     pd.Series( dtype='float' ),
                        'CBI':      pd.Series( dtype='float' ),
                        'Fop':      pd.Series( dtype='float' ),
                        'Putative': pd.Series( dtype='bool'  ),
                        'Length':   pd.Series( dtype='int'   ),
                        'Exp_rep1_norm':
                                    pd.Series( dtype='float' ),
                        'PA':       pd.Series( dtype='float' )},
                        index=allGenes )

    pa = getSpeciesPaxdbData( taxId )
    additionalIdentifiers = {}
    if pa and not args.PA_simple_mapping:
        additionalIdentifiers = createMappingForSpeciesProteins( pa.keys(), taxId )

    for geneId in allGenes:
        protId = geneNames2ProtName.get(geneId, None)

        # fill in prot-id
        df.loc[geneId, 'ProtId'] = protId

        # fill in experimental values (processed by me)
        if not protId is None:
            for n, exp_n in enumerate( ('Exp_rep1','Exp_rep2','Exp_rep3') ):
                if n<len(exps) and protId in exps[n]:
                    df.loc[geneId, exp_n] = exps[n][protId]

        for k in (protId,) + protName2GeneName[protId]: # + additionalIdentifiers.get( 
            if k in pa:
                df.loc[geneId, "PA"] = pa[k]
                #break          

        # fill in reference values
        if geneId in ref:
            for n, ref_n in enumerate( ('Ref_rep1','Ref_rep2','Ref_rep3') ):
                df.loc[geneId, ref_n] = ref[geneId][n]


    CDSseqs = getAllCDSSeqs( taxId, tuple(df.ProtId), args )
    assert(bool(CDSseqs))
    codonwDf = calcCodonwMeasures( CDSseqs )

    highlyExpressedGenes = [seq for seq in CDSseqs if seq.id in frozenset([geneNames2ProtName.get(x,None) for x in putativeHighlyExpressedGenes])]
    maxCDSforCAI = len(CDSseqs)
    if args.limit_CAI:
        maxCDSforCAI = args.limit_CAI
    caiValues = calcCAIusingReferenceGenes( CDSseqs[:maxCDSforCAI], highlyExpressedGenes, geneticCode=geneticCode )

    CDSlengths = dict([(seq.id, len(seq)) for seq in CDSseqs])
        
    #assert(len(caiValues)==len(df.index.values))

    for geneId in allGenes:
        protId = geneNames2ProtName.get(geneId, None)

        if protId in codonwDf.index:
            df.loc[geneId, 'Nc' ]  = codonwDf.loc[protId, 'Nc' ]
            df.loc[geneId, 'CAI0'] = codonwDf.loc[protId, 'CAI']
            df.loc[geneId, 'CBI']  = codonwDf.loc[protId, 'CBI']
            df.loc[geneId, 'Fop']  = codonwDf.loc[protId, 'Fop']
            
        if protId in caiValues:
            df.loc[geneId, 'CAI']  = caiValues[protId]  # codonwDf.loc[protId, 'CAI']
            

    for geneId in allGenes:
        protId = geneNames2ProtName.get(geneId, None)

        # Add 'putative' annotation (i.e., membership in references-set determined by annotation)
        df.loc[geneId, 'Putative'] = geneId in putativeHighlyExpressedGenes
        
        if protId in CDSlengths:
            df.loc[geneId, 'Length']   = CDSlengths[protId]

    df = df.assign( Exp_rep1_norm = lambda x: x.Exp_rep1 / x.Length )

    # Print correlations between all variables
    print(df.corr())

    dLFEdf = getGeneDLFEProfiles(taxId, args, protIds=allGenes)
    #dLFEdf.loc[:,290]
    x = pd.merge(df, dLFEdf, left_index=True, right_index=True)
    print( safeCorr( x.loc[:,0].values, x.loc[:,10].values ))
    print( safeCorr( x.loc[:,0].values, x.loc[:,'Exp_rep1'].values ))
    print( safeCorr( x.loc[:,'Exp_rep1'].values, x.loc[:,'Exp_rep2'].values  ))
    print( safeCorr( x.loc[:,'Exp_rep1'].values, x.loc[:,'Exp_rep3'].values  ))
    print( safeCorr( x.loc[:,'Exp_rep2'].values, x.loc[:,'Exp_rep3'].values  ))

    cutoff = 0.2
    if exps:
        print("====================")
        (pvalsLess, pvalsGreater) = annotatePvaluesHighVsLow( x, args, var='Exp_rep1', highLowPercentile=cutoff )
        plotStuff('Exp_rep1', 150, x, pvalsLess, pvalsGreater, args,   highLowPercentile=cutoff, taxId=taxId )

        print("====================")
        (pvalsLess, pvalsGreater) = annotatePvaluesHighVsLow( x, args, var='Exp_rep1', highLowPercentile=0.5 )
        plotStuff('Exp_rep1', 150, x, pvalsLess, pvalsGreater, args,   highLowPercentile=0.5, taxId=taxId )

        print("====================")
        (pvalsLess, pvalsGreater) = annotatePvaluesHighVsLow( x, args, var='Exp_rep1_norm', highLowPercentile=cutoff )
        plotStuff('Exp_rep1_norm', 150, x, pvalsLess, pvalsGreater, args,   highLowPercentile=cutoff, taxId=taxId )

        print("====================")
        (pvalsLess, pvalsGreater) = annotatePvaluesHighVsLow( x, args, var='Exp_rep1_norm', highLowPercentile=0.5 )
        plotStuff('Exp_rep1_norm', 150, x, pvalsLess, pvalsGreater, args,   highLowPercentile=0.5, taxId=taxId )
    
    print("====================")
    (pvalsLess, pvalsGreater) = annotatePvaluesHighVsLow( x, args, var='CAI',  highLowPercentile=0.5 )
    plotStuff('CAI',     150,  x, pvalsLess, pvalsGreater, args,   highLowPercentile=0.5, taxId=taxId )

    print("====================")
    if np.sum(~x.PA.isna().values) > 50:
        (pvalsLess, pvalsGreater) = annotatePvaluesHighVsLow( x, args, var='PA',  highLowPercentile=0.5 )
        plotStuff('PA',     150,  x, pvalsLess, pvalsGreater, args,   highLowPercentile=0.5, taxId=taxId )
    
    print("====================")
    (pvalsLess, pvalsGreater) = annotatePvaluesHighVsLow( x, args, var='Nc',   highLowPercentile=cutoff )
    plotStuff('Nc',     150,  x, pvalsLess, pvalsGreater, args,    highLowPercentile=cutoff, taxId=taxId )
    
    print("====================")
    (pvalsLess, pvalsGreater) = annotatePvaluesHighVsLow( x, args, var='CAI0', highLowPercentile=cutoff )
    plotStuff('CAI0',     150,  x, pvalsLess, pvalsGreater, args,  highLowPercentile=cutoff, taxId=taxId )
    
    print("====================")
    (pvalsLess, pvalsGreater) = annotatePvaluesHighVsLow( x, args, var='CBI', highLowPercentile=cutoff )
    plotStuff('CBI',     150,  x, pvalsLess, pvalsGreater, args,   highLowPercentile=cutoff, taxId=taxId )
    
    print("====================")
    (pvalsLess, pvalsGreater) = annotatePvaluesHighVsLow( x, args, var='Nc', highLowPercentile=cutoff )
    plotStuff('Nc',     150,  x, pvalsLess, pvalsGreater, args,    highLowPercentile=cutoff, taxId=taxId )


    if exps:
        plotLengthDiagnostics( df, taxId=taxId )

        # topGenes = x[ x.loc[:,'Exp_rep1'] > 5000 ].index.values
        # for geneId in topGenes:

        #     if geneId in putativeHighlyExpressedGenes:
        #         c = "***"
        #     else:
        #         c = "   "

        #     if geneId in geneDescriptions:
        #         print( "{} {} {}".format( geneId, c, geneDescriptions[geneId] ) )
        #     else:
        #         print( "{} {} NA".format( geneId, c ) )

    print("================= All done with taxid={} =================".format(taxId))

def parseList(conversion=str):
    def convert(values):
        return map(conversion, values.split(","))
    return convert


if __name__=="__main__":
    import sys
    import argparse
    
    argsParser = argparse.ArgumentParser()
    argsParser.add_argument( "--taxid",           type=parseList(int) )
    argsParser.add_argument( "--profile",         type=parseProfileSpec(), default=parseProfileSpec()('310:10:begin:0') )
    argsParser.add_argument( "--computation-tag", type=int,                default=Sources.RNAfoldEnergy_SlidingWindow40_v2 )
    argsParser.add_argument( "--shuffle-type",    type=int,                default=Sources.ShuffleCDSv2_python )
    argsParser.add_argument( "--num-shuffles",    type=int,                default=20 )
    argsParser.add_argument( "--limit-CAI",       type=int,                default=0 )
    argsParser.add_argument( "--PA-simple-mapping",  default=False, action="store_true" )
    args = argsParser.parse_args()

    for taxId in args.taxid:
        assert(int(taxId))
        assert(taxId in speciesConfig)
        config = speciesConfig[ taxId ]
        processGenome(args, taxId)
        
    sys.exit()
