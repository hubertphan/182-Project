import pandas as pd
from collections import defaultdict
from pathlib import Path
from pybedtools import BedTool


#define constants
#fgfr1 used for distance calculations
FGFR1_CHR = "chr8"

#these intervals are teh Reference build: GRCh37
FGFR1_START = 38268661
FGFR1_END = 38326153

#these intervals are the Reference build: GRCh38
# FGFR1_START = 38411143
# FGFR1_END = 38468635

#interval bed file, used in initial check to determine if ecDNA file contains FGFR1
FGFR1_BED = "fgfr1.bed.txt"
#FGFR1_BED = "FGFR1_pedi.bed.txt"

#ref genome, used in the second search
#finding all genes 
GENES_BED = "refseq_genes_GRCh37.bed" #for the GRCh37
#GENES_BED = "refseq_genes_GRCh38.bed" #for the GRCh38

#directory with the downloaded bedfiles
#ECDNA_DIR = "TRUE_ecDNA/"
#ECDNA_DIR = "PediEC/"

#testing
ECDNA_DIR = "ecdna_fixed/"


fgfr1_positive_samples = 0
#create dictionary to track specific gene count
def make_entry():
    #we need to take into account gene count and distance
    #store distances as a list, could be multiple
    return {"count":0, "distances": []}
geneDict = defaultdict(make_entry)

def run_bedTools_initial(ecDNA_bed, FGFR1_Bed):
    #this checks if FGFR1 is in ecDNA file
    #return a bool
    ec = BedTool(str(ecDNA_bed))
    return bool(ec.intersect(str(FGFR1_Bed), u=True))


def run_bedTools_further(ecDNA_bed, refgenome):
    result = []

    ref = BedTool(str(refgenome))
    overlap = ref.intersect(str(ecDNA_bed), wa=True)
    #we need to now the chrom where the found gene came from
    #start 
    #end 
    #name
    for i in overlap:
        chrom = i.chrom
        start = int(i.start)
        end = int(i.end)
        name = i.name
        result.append((chrom, start, end, name))
    return result

for ecDNA_bed in Path(ECDNA_DIR).glob("*.bed"):
    #first check if ecDNA contains FGFR1
    #we technically do nott need this IF we know we donwloaded the correct ecDNA files (contains FGFR1)
    # if not run_bedTools_initial(ecDNA_bed, FGFR1_BED):
    #     continue
    # #we now now the ecDNA does have FGFR1
    fgfr1_positive_samples += 1

    #now find all the genes 
    overlapping_genes = run_bedTools_further(ecDNA_bed, GENES_BED)
    
    for chrom, start, end, gene in overlapping_genes:
        #will only have a distance IF on the same chrom as human FGFR1
        #hard to tell distances for genes on different chromosomes
        if chrom == FGFR1_CHR:
            dist = min(abs(start - FGFR1_START), abs(end - FGFR1_END))
        else:
            dist = float("inf")

        geneDict[gene]["count"] += 1
        geneDict[gene]["distances"].append(dist)

#output to a csv
rows = []
for gene, stats in geneDict.items():
    min_dist = min(stats["distances"]) #i think this is ok, it will help us know which gene it is later on anwyays
    count = stats["count"]
    fraction = count / fgfr1_positive_samples
    rows.append([gene, min_dist, count, fraction])
    #rows.append([gene, min_dist, count, fraction_allEC])

df = pd.DataFrame(rows, columns=["Gene", "Distance", "Co-occurrence Count","Fraction in FGFR1 positive samples"])
df.to_csv("grch37EVERYTHING_fgfr1_cooccurrence_table.csv", index=False)

