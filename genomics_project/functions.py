"""
Core utility functions for the phylogenomic tree-building pipeline.
"""

import subprocess
import logging
from pathlib import Path
import random
import pandas
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo._io import read
from Bio.Phylo._utils import draw
import matplotlib.pyplot as plot


def run_prodigal(input_file, prodigal_output_dir):
    """
    Predicts genes using Prodigal.
    """
    # 1. Define Output paths
    output_path_prot = prodigal_output_dir / input_file.with_suffix(".faa").name
    output_path_nucl = prodigal_output_dir / input_file.with_suffix(".fna").name

    # 2. Run Prodigal
    if not output_path_prot.exists():
        logging.info("--- Predicting genes for: %s ---", input_file.stem)
        subprocess.run(
            [
                "prodigal",
                "-i",
                str(input_file),
                "-a",
                str(output_path_prot),  # Proteins
                "-d",
                str(output_path_nucl),  # Nucleotides
                "-p",
                "single",
            ],
            check=True,
        )

    return output_path_prot, output_path_nucl


def make_blast_db(input_file, db_dir_nucl, db_dir_prot, db_type):
    """
    Creates a BLAST database.
    """

    if db_type == "prot":
        output_db = db_dir_prot / input_file.stem
        index_suffix = ".pin"

    elif db_type == "nucl":
        output_db = db_dir_nucl / input_file.stem
        index_suffix = ".nin"

    else:
        logging.error("Unknown db_type: %s!", db_type)

    # Check if DB already exists
    if (output_db.with_suffix(index_suffix)).exists():
        logging.info("Database %s (%s) already exists, skipping.", input_file, db_type)
        return

    # Create DB
    logging.info("Creating BLAST DB for %s (%s)", input_file.stem, db_type)

    subprocess.run(
        [
            "makeblastdb",
            "-in",
            str(input_file),
            "-dbtype",
            db_type,
            "-out",
            str(output_db),
        ],
        check=True,
    )


def run_blast(db_name, reference_path, blast_program, blast_db_dir, blast_out):
    """
    Run a single BLAST job (either BLASTP or BLASTN).
    """
    logging.info("Running %s against %s", blast_program.upper(), db_name)

    suffix = ".blastp.txt" if blast_program == "blastp" else ".blastn.txt"
    output_file = blast_out / f"{db_name}{suffix}"

    try:
        subprocess.run(
            [
                blast_program,
                "-query",
                str(reference_path),
                "-db",
                str(blast_db_dir / db_name),
                "-max_target_seqs",
                "1",
                "-evalue",
                "1e-5",
                "-outfmt",
                "6 qseqid sseqid pident evalue sstart send length qcovs",
                "-out",
                str(output_file),
                "-num_threads",
                "1",
            ],
            check=True,
        )
    except subprocess.CalledProcessError as e:
        logging.error("%s failed for %s: %s", blast_program.upper(), db_name, e)


# Filter and Parse BLAST Text Files
def parse_blast_results(file_path, id_cutoff, ev_cutoff, cov_cutoff):
    """
    Reads a BLAST-Table and filters quality.
    """

    # Define colums
    columns = [
        "qseqid",  # Query Gene ID
        "sseqid",  # Subject Gene ID
        "pident",  # Identity %
        "evalue",  # E-Value
        "sstart",  # Alignment Start
        "send",  # Alignment End
        "length",  # Alignment Length
        "qcovs",  # Query Coverage
    ]

    # 1. Read file
    try:
        data_frame = pandas.read_csv(file_path, sep="\t", names=columns)
    except Exception as e:
        logging.error(
            "Could not read %s: %s - Empty, damaged or not found", file_path, e
        )
        return pandas.DataFrame()  # returns empty table

    # 2. Filter
    filtered_file = data_frame[
        (data_frame["pident"] >= id_cutoff)
        & (data_frame["evalue"] <= ev_cutoff)
        & (data_frame["qcovs"] >= cov_cutoff)
    ]
    return filtered_file


def build_sets(master_table, n_sets=5):
    """
    Build sets (randomly by selecting n_sets unique homolog IDs.
    Returns List of ID-Tuples ([("Ref, Homolog), ...]).
    """
    if master_table.empty:
        logging.warning("Master table is empty.")
        return []

    # All genomes in BLAST result
    all_genomes = set(master_table["genome"].unique())

    # Group to reference gene (qseqid)
    groups = master_table.groupby("qseqid")
    valid_groups = []

    min_required = len(all_genomes) - 1

    for qseqid, group in groups:
        genomes_present = set(group["genome"])

        if len(genomes_present) >= min_required:
            valid_groups.append(group)

    if len(valid_groups) == 0:
        logging.error("No complete ortholog groups found.")
        return []

    # Randomly pick n_sets orthologous groups
    random.seed(11)  # to be reproducable
    selected = random.sample(valid_groups, min(n_sets, len(valid_groups)))

    final_sets = []
    for group in selected:
        seq_pairs = list(zip(group["qseqid"], group["sseqid"], group["genome"]))
        final_sets.append(seq_pairs)

    return final_sets


def write_sets_to_fasta(set_items, fasta_dir, fasta_ext, output_file, ref_genome_name):
    """
    Concetenates ortholog sequence with reference sequence in one file.
    Allow missing genome of max. 1
    """

    # Load sequences from FASTA file
    ortholog_database = {}
    for fasta in Path(fasta_dir).glob(f"*{fasta_ext}"):
        ortholog_database.update(SeqIO.to_dict(SeqIO.parse(fasta, "fasta")))

    records_to_write = []

    # Add Reference sequence
    if set_items:
        # Index 0 = qseqid
        ref_id = set_items[0][0]

        if ref_id in ortholog_database:
            ref_record = ortholog_database[ref_id]
            ref_record.id = f"{ref_genome_name} | {ref_id}"
            ref_record.description = ""
            records_to_write.append(ref_record)
        else:
            logging.warning("Reference sequence %s not found in database!", ref_id)

    # Add Homologes
    for qseqid, sseqid, genome in set_items:
        if sseqid not in ortholog_database:
            logging.warning("Reference ID %s not found", sseqid)
            continue
        seq_record = ortholog_database[sseqid]
        seq_record.id = f"{genome} | {sseqid}"
        seq_record.description = ""
        records_to_write.append(seq_record)

        if not records_to_write:
            logging.warning("No sequences found for this set %s!", genome)

    SeqIO.write(records_to_write, output_file, "fasta")
    logging.info(
        "Generated fasta file: %s with %d sequences",
        output_file,
        len(records_to_write),
    )


def complete_fasta(input_fasta, output_fasta, all_genomes):
    """
    Set all Fasta Files to 10 genomes -> if 1 missing: Fill in gaps.
    """

    records = {
        rec.id.split("|")[0].strip(): rec for rec in SeqIO.parse(input_fasta, "fasta")
    }

    target_len = (
        len(next(iter(records.values())).seq) if records else 1
    )  # lenght of placeholder
    placeholder = "X" if "protein" in str(input_fasta).lower() else "N"

    final_records = []
    for genome in sorted(all_genomes):
        if genome in records:
            final_records.append(records[genome])
        else:
            # Hier ist dein Gap-Filling mit korrekter Länge
            final_records.append(
                SeqRecord(
                    Seq(placeholder * target_len), id=f"{genome} | gap", description=""
                )
            )
            logging.info(
                "Missing genome %s filled with gap in %s", genome, input_fasta.name
            )

    SeqIO.write(final_records, output_fasta, "fasta")


def run_muscle(input_file, output_dir):
    """
    Runs MUSCLE alignment for a single file.
    """
    input_path = Path(input_file)
    output_path = Path(output_dir)
    output_afa = output_path / f"{input_path.stem}.afa"

    # Skip if already done
    if output_afa.exists():
        logging.info("Alignment exists for %s, skipping.", input_path.name)
        return

    try:
        subprocess.run(
            [
                "muscle",
                "-align",
                str(input_path),
                "-output",
                str(output_afa),
            ],
            check=True,
        )
        logging.info("Aligned: %s", input_path.name)

    except subprocess.CalledProcessError as e:
        logging.error("MUSCLE failed for %s: %s", input_path.name, e)


def get_filtered_alignment(input_path):
    """Filters out sequences consisting only of N or X using SeqIO."""
    # Wir lesen die Datei als Liste von Records
    records = list(SeqIO.parse(str(input_path), "fasta"))

    # Der Filter bleibt derselbe: Zähle alles außer N und X
    filtered = [
        rec
        for rec in records
        if len(str(rec.seq).upper().replace("N", "").replace("X", "")) > 10
    ]
    return filtered


def run_iqtree(aln_file, out_dir, threads="1"):
    """
    Making tree files with IQTree.
    """
    output_file = out_dir / f"{aln_file.stem}.treefile"
    if output_file.exists() and output_file.stat().st_size > 0:
        return logging.info(f"Tree for {aln_file} already exists.")

    try:
        subprocess.run(
            [
                "iqtree",
                "-s",
                str(aln_file),
                "-m",
                "MFP",
                "-bb",
                "1000",
                "-T",
                str(threads),
                "-pre",
                str(out_dir / aln_file.stem),
            ],
            check=True,
        )
        logging.info("Creating tree for %s", aln_file)
    except subprocess.CalledProcessError as e:
        logging.error("Tree building failed for %s: %s", aln_file, e)


def run_iqtree_filtered(aln_file, out_dir, threads="1"):
    """
    Helperfunction: Filters the alignment first, then runs IQ-TREE.
    """
    # 1. Filtern
    filtered_recs = get_filtered_alignment(aln_file)

    if len(filtered_recs) < 3:
        return logging.warning("Skipping %s: Not enough data.", aln_file.name)

    # 2. temporary file
    temp_aln = out_dir / f"{aln_file.stem}_temp.afa"

    try:
        SeqIO.write(filtered_recs, str(temp_aln), "fasta")

        # call defined function
        run_iqtree(temp_aln, out_dir, threads)

    finally:
        if temp_aln.exists():
            temp_aln.unlink()


def build_supermatrix(alignment_files, output_file):
    """
    Concetenates alignments into a supermatrix.
    """
    # Read all alignments into a list of dicts: {taxon: sequence}
    logging.info("Building a supermatrix for %s", output_file.name)
    if output_file.exists() and output_file.stat().st_size > 0:
        return logging.info("Tree for %s already exists.", output_file)

    alignments = []
    taxa_order = None

    for aln in alignment_files:
        records = list(SeqIO.parse(aln, "fasta"))
        aln_dict = {}

        lengths = {len(rec.seq) for rec in records}
        if len(lengths) > 1:
            logging.error("Inconsistent alignment lengths in file %s: %d", aln, lengths)

        for rec in records:
            genome_id = (
                rec.id.split("|")[0].strip() if "|" in rec.id else rec.id.strip()
            )
            if genome_id in aln_dict:
                logging.error("Duplicate taxon %d in %s", rec.id, aln)
            aln_dict[rec.id] = str(rec.seq)

        # Establish taxon order from the first alignment
        if taxa_order is None:
            taxa_order = list(aln_dict.keys())
            logging.info("Taxa order establishe with %d for %s", len(taxa_order), aln)

        alignments.append(aln_dict)

    logging.info("All alignment loaded, building concatenated sequences.")
    concatenated_records = []
    for taxon in taxa_order:
        try:
            concat_seq = "".join(aln[taxon] for aln in alignments)
            rec = SeqRecord(Seq(concat_seq), id=taxon, description="")
            concatenated_records.append(rec)
        except:
            logging.error("Taxon %s missing in one of the alignments!", taxon)

    logging.info("Writing supermatrix to %s", output_file)
    SeqIO.write(concatenated_records, output_file, "fasta")
    logging.info("Supermatrix construction completed successfully!")


def plot_tree(tree_file_path, output_dir_path):
    tree = read(tree_file_path, "newick")

    fig = plot.figure(figsize=(10, 20))
    ax = fig.add_subplot(1, 1, 1)

    draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=True,
    )

    output_png = output_dir_path / f"{tree_file_path.stem}.png"
    plot.savefig(output_png, dpi=300, bbox_inches="tight")
    logging.info("%s tree saved.", output_png)
    plot.close()
