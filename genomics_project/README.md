# Phylogeny of Prevotella

Genetic variation and evolutionary reconstruction of selected *Prevotella* species, implemented as a single-script Python bioinformatics pipeline.

## Overview

*Prevotella* is a genus of anaerobic bacteria within the family *Prevotellaceae*, playing a central role in carbohydrate and protein fermentation in mammalian microbiomes. Given ongoing taxonomic revisions in the genus and the known influence of cattle contact on the human microbiome, this project analyzes eight human-associated and two cattle-associated *Prevotella* taxa to assess evolutionary relationships across host niches.

## Repository Structure

- `tree_pipeline.py` — main pipeline script (gene prediction, BLAST search, alignment, tree construction)
- `functions.py` — helper functions used by the main pipeline
- `check_filter_settings.py` — script used to iteratively test and determine BLAST filtering thresholds
- `filter.log` — log output from filter threshold testing

## Pipeline

The full workflow was implemented within a single Python script:

1. **Data retrieval** — Ten *Prevotella* genomes downloaded from NCBI GenBank.
2. **Gene prediction** — Structural gene prediction via **Prodigal**, producing annotated protein-coding sequences (`.fna` nucleotide / `.faa` amino acid).
3. **Reference selection** — One genome selected as reference (highest number of predicted genes; used consistently for both nucleotide and protein analysis).
4. **Ortholog search** — All protein-coding genes from the reference genome used as queries in **BLASTp** (protein) and **BLASTn** (nucleotide) against the remaining nine genomes.
5. **Filtering** — BLAST hits filtered on E-value, sequence identity, and query coverage to distinguish true orthologs from paralogs/artifacts. Final thresholds (Proteins: 100% coverage, 70% identity; Nucleotides: 70% coverage, 65% identity) were determined through iterative testing (`check_filter_settings.py`, `filter.log`).
6. **Homolog set selection** — Five protein and five nucleotide homologous sets randomly selected (to avoid gene-family/genome bias). Genomes missing from a set (9/10 present) were gap-filled with `N` (nucleotide) or `X` (amino acid) placeholders to preserve alignment structure.
7. **Alignment** — All homolog sets aligned using **MUSCLE**.
8. **Phylogenetic reconstruction** — Individual gene trees built with **IQ-TREE**; nucleotide and protein alignments each concatenated into a "supermatrix" and re-analyzed with IQ-TREE for a robust, genome-scale phylogeny.
9. **Visualization** — Resulting Newick tree files visualized using **BioPython's Phylo module** (tip labels, branch lengths, bootstrap values).
10. **Functional verification** — Orthology of extracted markers confirmed via BLASTp (protein) / BLASTx (nucleotide) searches.

## Tools & Software

- Python, BioPython (Phylo module)
- Prodigal (gene prediction)
- BLAST (BLASTp, BLASTn, BLASTx)
- MUSCLE (multiple sequence alignment)
- IQ-TREE (phylogenetic inference)
- BioEdit Sequence Aligner (alignment visualization)

## Key Results

- Functional annotation confirmed all extracted markers as core genes (except Set 1), suitable for species-level phylogenetics. Nucleotide-level sets recovered exclusively ribosomal genes — consistent with these being the "hard core" genes stable enough for reliable phylogenetic signal at that level.
- Strict filtering thresholds successfully excluded paralogs and artificial matches; resulting phylogenies show distinct habitat-associated lineages with moderate-to-high bootstrap support.

## Identified Marker Genes

**Protein sets:**
| Set | Protein | Conservation |
|---|---|---|
| 1 | MlaE family ABC transporter permease | medium–high |
| 2 | UDP-glucose 4-epimerase (GalE) | high |
| 3 | Asparagine–tRNA ligase | very high |
| 4 | Do family serine endopeptidase (HtrA) | high |
| 5 | 50S ribosomal protein L25 | very high |

**Nucleotide sets:**
| Set | Protein | Conservation |
|---|---|---|
| 1 | 30S ribosomal protein S6 | very high |
| 2 | 50S ribosomal protein L11 | very high |
| 3 | 50S ribosomal protein L2 | very high |
| 4 | 30S ribosomal protein S12 | very high |
| 5 | 50S ribosomal protein L3 | very high |


## References

1. Betancur-Murillo, C. L., Aguilar-Marín, S. B., & Jovel, J. (2023). Prevotella: A Key Player in Ruminal Metabolism. *Microorganisms*, 11(1), 1.
2. Cuperus, T. et al. (2025). Prevotella as the main driver for the association between dairy farming and human gut microbiome composition. *Frontiers in Microbiomes*, 4.
3. Hitch, T. C. A. et al. (2022). A taxonomic note on the genus Prevotella. *Systematic and Applied Microbiology*, 45(6), 126354.
4. Tett, A. et al. (2021). Prevotella diversity, niches and interactions with the human host. *Nature Reviews Microbiology*, 19(9), 585–599.
5. NCBI BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi
