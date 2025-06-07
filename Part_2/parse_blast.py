#!/usr/bin/env python3

from Bio.Blast import NCBIXML

# open the BLAST XML
with open("FGFR1_vs_ecDNA.xml") as handle:
    blast_qresult = NCBIXML.read(handle)

covered_positions = set()
for alignment in blast_qresult.alignments:
    for hsp in alignment.hsps:
        covered_positions.update(range(hsp.query_start, hsp.query_end+1))

total_len = blast_qresult.query_length
covered = len(covered_positions)
print(f"Protein length: {total_len}")
print(f"Residues covered by ecDNA hits: {covered}")
print(f"Coverage: {covered/total_len:.1%}")

all_pos   = set(range(1, 821+1))
missing   = sorted(all_pos - covered_positions)
print("Missing positions:", missing)

runs = []
for p in missing:
    if not runs or p != runs[-1][-1] + 1:
        runs.append([p])
    else:
        runs[-1].append(p)
print("Contiguous missing runs:", runs)
