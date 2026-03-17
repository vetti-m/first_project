from collections import Counter
import logging
import functions
from pathlib import Path
import pandas as pd

INPUT_DIR = Path("dataflow/01-nucl")
PRODIGAL_OUTPUT_DIR = Path("dataflow/02-prot")
BLAST_DB_DIR_PROT = Path("dataflow/03-blast-db-prot")
BLAST_DB_DIR_NUCL = Path("dataflow/03-blast-db-nucl")
BLAST_OUT = Path("dataflow/04-blast-out")
SETS_DIR = Path("dataflow/05-sets")
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
    "blastp": {
        "evalue": 1e-10,
        "identity": 70,  # in %
        "coverage": 100,  # in %
    },
    "blastn": {
        "evalue": 1e-5,
        "identity": 65,  # in %
        "coverage": 70,  # in %
    },
}

logging.basicConfig(
    filename="filter.log",
    filemode="w",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

all_proteins_filtered = []
all_nucleotides_filtered = []

for blast_file in BLAST_OUT.glob("*.txt"):
    if ".blastp." in blast_file.name:
        tool = "blastp"
    elif ".blastn." in blast_file.name:
        tool = "blastn"
    else:
        continue  # Skip files that are not correctly saved

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

try:
    master_protein_table = pd.concat(all_proteins_filtered)
except ValueError:
    master_protein_table = pd.DataFrame()

try:
    master_nucl_table = pd.concat(all_nucleotides_filtered)
except ValueError:
    master_nucl_table = pd.DataFrame()


def log_genome_distribution(table, label):
    counts = Counter()
    for qseqid, group in table.groupby("qseqid"):
        n = len(set(group["genome"]))
        counts[n] += 1

    logging.info("=== %s genome distribution ===", label)
    for n in sorted(counts.keys()):
        logging.info("  %d genomes: %d genes", n, counts[n])


log_genome_distribution(master_protein_table, "Protein")
log_genome_distribution(master_nucl_table, "Nucleotide")
