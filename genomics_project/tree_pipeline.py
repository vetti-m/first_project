#!/usr/bin/env python

import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pandas as pd
import functions


# -------- 1. Define variables --------

# Define the number of threads/tasks that can be done simultaneously
THREADS = 4

# Define your input and output directories
INPUT_DIR = Path("dataflow/01-nucl")
PRODIGAL_OUTPUT_DIR = Path("dataflow/02-prot")
BLAST_DB_DIR_PROT = Path("dataflow/03-blast-db-prot")
BLAST_DB_DIR_NUCL = Path("dataflow/03-blast-db-nucl")
BLAST_OUT = Path("dataflow/04-blast-out")
SETS_DIR = Path("dataflow/05-sets")
SET_PREP = Path("dataflow/05-sets_complete")
ALIGNMENT_DIR = Path("dataflow/06-alignment")
TREE_DIR = Path("dataflow/07-iqtree")
FINAL_TREE_DIR = Path("dataflow/08-final_tree")

# Define the extension of the input files
NUCL_EXT = ".fasta"
PROT_EXT = ".faa"
NUCL_OUT_EXT = ".fna"
ALN_EXT = ".afa"

# Cutoffs
CUT_OFFS = {
    "blastp": {"evalue": 1e-10, "identity": 70, "coverage": 100},  # in %  # in %
    "blastn": {"evalue": 1e-5, "identity": 65, "coverage": 70},  # in %  # in %
}

# Logging:
logging.basicConfig(
    filename="tree_pipeline.log",
    filemode="w",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


# -------- 2. Make directories --------

# Ensure the input direcotry exists
logging.info("--- Checking/Creating folders. ---")
INPUT_DIR.mkdir(parents=True, exist_ok=True)

# Ensure the output directories exists
PRODIGAL_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
BLAST_DB_DIR_PROT.mkdir(parents=True, exist_ok=True)
BLAST_DB_DIR_NUCL.mkdir(parents=True, exist_ok=True)
BLAST_OUT.mkdir(parents=True, exist_ok=True)
SETS_DIR.mkdir(parents=True, exist_ok=True)
SET_PREP.mkdir(parents=True, exist_ok=True)
ALIGNMENT_DIR.mkdir(parents=True, exist_ok=True)
TREE_DIR.mkdir(parents=True, exist_ok=True)
FINAL_TREE_DIR.mkdir(parents=True, exist_ok=True)

# Prepare a list of input files with .fasta extension
list_input_files = list(INPUT_DIR.glob(f"*{NUCL_EXT}"))

# print(f"Input Files: {list_input_files}")
logging.info("Inputfiles found: %d", len(list_input_files))

# -------- 3. Define parallel Threads ---------

# Run prodigal in parallel using ThreadPoolExecutor with max. 4 threads
logging.info("--- Starting Prodigal Gene Prediction. ---")
if PRODIGAL_OUTPUT_DIR.exists() and any(PRODIGAL_OUTPUT_DIR.iterdir()):
    logging.info("Annotation files already exists.")
else:
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        executor.map(
            functions.run_prodigal,
            list_input_files,
            [PRODIGAL_OUTPUT_DIR] * len(list_input_files),
        )


# -------- 4. Make Blast Databases --------

list_prot_files = list(PRODIGAL_OUTPUT_DIR.glob(f"*{PROT_EXT}"))
list_nucl_files = list(PRODIGAL_OUTPUT_DIR.glob(f"*{NUCL_OUT_EXT}"))

# --- Create protein BLAST databases ---
with ThreadPoolExecutor(max_workers=THREADS) as executor:
    executor.map(
        functions.make_blast_db,
        list_prot_files,
        [BLAST_DB_DIR_NUCL] * len(list_prot_files),
        [BLAST_DB_DIR_PROT] * len(list_prot_files),
        ["prot"] * len(list_prot_files),
    )

# --- Create nucleotide BLAST databases ---
with ThreadPoolExecutor(max_workers=THREADS) as executor:
    executor.map(
        functions.make_blast_db,
        list_nucl_files,
        [BLAST_DB_DIR_NUCL] * len(list_nucl_files),
        [BLAST_DB_DIR_PROT] * len(list_nucl_files),
        ["nucl"] * len(list_nucl_files),
    )


# -------- 5. Decide what's the best Reference genome - most genes--------

logging.info("--- Counting genes per Genome. ---")
gene_stats = []
list_faa_files = list(PRODIGAL_OUTPUT_DIR.glob(f"*{PROT_EXT}"))
for f in list_faa_files:
    with open(f, "r", encoding="utf-8") as file:
        # oneliner for gene_counter (Genename starts with ">") in file
        count = sum(1 for line in file if line.startswith(">"))
        gene_stats.append((count, f))  # tuple with number of genes + filename
        logging.info("Found %d in %s", count, f.stem)

reference_path_prot = max(gene_stats)[1]  # with most genes

# nucleotide form (same genes)
reference_path_nucl = reference_path_prot.with_suffix(NUCL_OUT_EXT)

logging.info("Reference Genome is: %s", reference_path_prot.stem)


# -------- 6. BLAST --------

logging.info("--- Starting BLAST ---")
all_db_files_prot = [
    filename.stem
    for filename in BLAST_DB_DIR_PROT.glob("*.pin")
    if filename.stem != reference_path_prot.stem
]

all_db_files_nucl = [
    filename.stem
    for filename in BLAST_DB_DIR_NUCL.glob("*.nin")
    if filename.stem != reference_path_prot.stem
]

logging.info("Found %d databases to search against (AA).", len(all_db_files_prot))
logging.info("Found %d databases to search against (NUCL.).", len(all_db_files_prot))

# BLASTP
if BLAST_OUT.exists() and any(BLAST_OUT.glob("*.blastp.txt")):
    logging.info("BLASTP results already exist. Skipping BLAST.")
else:
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        executor.map(
            functions.run_blast,
            all_db_files_prot,
            [reference_path_prot] * len(all_db_files_prot),
            ["blastp"] * len(all_db_files_prot),
            [BLAST_DB_DIR_PROT] * len(all_db_files_prot),
            [BLAST_OUT] * len(all_db_files_prot),
        )

# BLASTN
if BLAST_OUT.exists() and any(BLAST_OUT.glob("*.blastn.txt")):
    logging.info("BLASTN results already exist. Skipping BLAST.")
else:
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        executor.map(
            functions.run_blast,
            all_db_files_nucl,
            [reference_path_nucl] * len(all_db_files_nucl),
            ["blastn"] * len(all_db_files_nucl),
            [BLAST_DB_DIR_NUCL] * len(all_db_files_nucl),
            [BLAST_OUT] * len(all_db_files_nucl),
        )


# -------- 7. Filter and Parse Blast results --------

logging.info("--- Starting Filtering ---")

all_proteins_filtered = []
all_nucleotides_filtered = []

for blast_file in BLAST_OUT.glob("*.txt"):
    if ".blastp." in blast_file.name:
        tool = "blastp"
    elif ".blastn." in blast_file.name:
        tool = "blastn"
    else:
        logging.info("Right tool not found.")

    # import Cutoffs
    rule = CUT_OFFS[tool]

    # Call parse and filter function
    data_frame_filtered = functions.parse_blast_results(
        blast_file,
        id_cutoff=rule["identity"],  # removes distant paralogs
        ev_cutoff=rule["evalue"],  # sorts out random matches and artifacts
        cov_cutoff=rule["coverage"],  # most important orthology filter ("full length")
    )

    genome_name = blast_file.name.split(".")[0]
    data_frame_filtered["genome"] = genome_name

    # Store output
    if not data_frame_filtered.empty:
        if tool == "blastp":
            all_proteins_filtered.append(data_frame_filtered)
        else:
            all_nucleotides_filtered.append(data_frame_filtered)
    else:
        logging.error("No results for %s", blast_file.name)


# -------- 8. Selecting 5 random sequences --------

# Combine all individual species results into two master tables
try:
    master_protein_table = pd.concat(all_proteins_filtered)
except ValueError:
    logging.error("No protein BLAST results to concatenate.")
    master_protein_table = pd.DataFrame()

try:
    master_nucl_table = pd.concat(all_nucleotides_filtered)
except ValueError:
    logging.error("No protein BLAST results to concatenate.")
    master_nucl_table = pd.DataFrame()

logging.info(
    "FILTER SUCCESS: %d protein rows and %d nucleotide rows stored in master tables.",
    len(master_protein_table),
    len(master_nucl_table),
)
logging.info("--- Creating sets. ---")

# Select random homologs for 5 sets - CAVE: returns homolog IDs not sequences
protein_sets = functions.build_sets(master_protein_table)
nucleotide_sets = functions.build_sets(master_nucl_table)

logging.info(
    "Final selection: %d protein sets and %d nucleotide sets.",
    len(protein_sets),
    len(nucleotide_sets),
)

# -------- 9. Concetenate the final_selection in one file --------

if SETS_DIR.exists() and any(SETS_DIR.iterdir()):
    logging.info("Sets already exists.")
else:
    for i, set_items in enumerate(protein_sets, start=1):
        out_prot_file = SETS_DIR / f"protein_set_{i}.fasta"
        functions.write_sets_to_fasta(
            set_items,
            PRODIGAL_OUTPUT_DIR,
            PROT_EXT,
            out_prot_file,
            reference_path_prot.stem,
        )

    for i, set_items in enumerate(nucleotide_sets, start=1):
        out_nucl_file = SETS_DIR / f"nucl_set_{i}.fasta"
        functions.write_sets_to_fasta(
            set_items,
            PRODIGAL_OUTPUT_DIR,
            NUCL_OUT_EXT,
            out_nucl_file,
            reference_path_nucl.stem,
        )


# -------- 10. Fill in Gaps for missing genomes --------
logging.info("--- Filling in gaps. ---")

genomes_prot = (
    set(master_protein_table["genome"].unique())
    if not master_protein_table.empty
    else set()
)
genomes_nucl = (
    set(master_nucl_table["genome"].unique()) if not master_nucl_table.empty else set()
)

all_genomes = (
    genomes_prot
    | genomes_nucl
    | {str(reference_path_prot.stem), str(reference_path_nucl.stem)}
)

logging.info("Total genomes to account for: %d", len(all_genomes))

if SET_PREP.exists() and any(SET_PREP.iterdir()):
    logging.info("Sets (with gaps) already exists.")
else:
    for fasta in SETS_DIR.glob(f"*{NUCL_EXT}"):
        output_fasta = SET_PREP / fasta.name
        functions.complete_fasta(fasta, output_fasta, all_genomes)
        logging.info("Checked: %s", fasta)


# -------- 11. Starting MUSCLE Alignment --------

logging.info("--- Starting MUSCLE Alignment ---")
unaligned_files = list(SET_PREP.glob(f"*{NUCL_EXT}"))

if ALIGNMENT_DIR.exists() and any(ALIGNMENT_DIR.iterdir()):
    logging.info("Sets (with gaps) already exists.")
else:
    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        executor.map(
            functions.run_muscle,
            unaligned_files,
            [ALIGNMENT_DIR] * len(unaligned_files),
        )
    logging.info("Alignment finished!")


# -------- 12. Start Building a tree for each Set --------

logging.info("--- Starting Building Trees with IQTREE for each set! ---")

aligned_files = list(ALIGNMENT_DIR.glob(f"*{ALN_EXT}"))

with ThreadPoolExecutor(max_workers=THREADS) as executor:
    executor.map(
        functions.run_iqtree_filtered,
        aligned_files,
        [TREE_DIR] * len(aligned_files),
        ["1"] * len(aligned_files),
    )

# --------- 13. Building a supermatrix --------

prot_alignments = [f for f in aligned_files if "protein_set" in f.name]
nucl_alignments = [f for f in aligned_files if "nucl_set" in f.name]

if prot_alignments:
    prot_matrix_file = TREE_DIR / "supermatrix_proteins.fasta"
    functions.build_supermatrix(prot_alignments, prot_matrix_file)
    logging.info(
        "Building supermatrix (Protein) from %d alignments.", len(prot_alignments)
    )
if nucl_alignments:
    nucl_matrix_file = TREE_DIR / "supermatrix_nucleotides.fasta"
    functions.build_supermatrix(nucl_alignments, nucl_matrix_file)
    logging.info(
        "Building supermatrix (Nucleotides) from %d alignments.", len(nucl_alignments)
    )

# -------- 14. Building a tree of the supermatrix --------

logging.info("Building trees of the supermatrix")

supermatrix_files = [
    TREE_DIR / "supermatrix_proteins.fasta",
    TREE_DIR / "supermatrix_nucleotides.fasta",
]

for sm_file in supermatrix_files:
    tree_result = sm_file.with_suffix(".fasta.treefile")

    if tree_result.exists() and tree_result.stat().st_size > 0:
        logging.info("Final tree for %s already exists.", sm_file.name)
    else:
        logging.info("Building final tree for supermatrix: %s", sm_file.name)
        functions.run_iqtree(sm_file, TREE_DIR)


# ------ 15. Visualize the tree ------

logging.info("--- Creating Tree as Image (PNG) ---")

trees = list(TREE_DIR.glob("*.treefile"))

for treefile in trees:
    logging.info("Creating Tree for %s", treefile)
    functions.plot_tree(treefile, FINAL_TREE_DIR)
