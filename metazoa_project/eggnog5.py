"""
EggNOG Evolutionary Analysis Library
-------------------------------------

This module provides functions for processing eggNOG data, including identification
of homologous genes between different species, orthologous groups, and functional annotation.


Dependencies:
 - pandas: For time-efficient data manipulation.
 - pathlib: For cross-platform file system operations.
 - logging: For tracking process execution and errors.

"""

import logging
from pathlib import Path
import pandas as pd


def retrieve_taxid(taxid_document: Path) -> dict:
    """
    Retrieves the TaxID from a given species out of the given file from the eggNOG DB.

    Args:
        taxid_document (Path): Path to the taxid information file, e.g. "data/e5.taxid_info.tsv"

    Returns:
        dict: A dictionary mapping scientific names to their Tax IDs.
    """
    # Store File in a Dataframe
    try:
        taxid = pd.read_csv(
            str(
                taxid_document
            ),  # so even if a path would be passed it is a str afterwards
            sep="\t",  # Columns seperated by Tab
            header=0,  # first row (Index 0)
            usecols=[0, 1, 4],  # only get needed columns
            names=[
                "TaxID",  # from species
                "Scientific Name",
                "TaxID Lineage",  # from kingdom
            ],
        )
        logging.info("The tsv-file %s could be read successfully!", taxid_document)
    # raise not necessary, since the next code would not work without
    except FileNotFoundError:
        logging.error("Could not read the file!")
    except Exception as e:
        logging.error("An Error occured: %s", e)

    # get the rows from Metazoa (ID 33208)
    metazoa_df = taxid[taxid["TaxID Lineage"].str.contains("\\b33208\\b", na=False)]

    # safe Scientific Name and TaxID in a dictionary
    taxid_dict = metazoa_df.set_index("Scientific Name")["TaxID"].to_dict()
    logging.info(
        "Successfully created dictionary with %d TaxIDs belonging to Metazoa, "
        "mapped to their scientific name",
        len(taxid_dict),
    )
    return taxid_dict


def parse_members_file(filename: Path) -> pd.DataFrame:
    """
    Parses eggNOG 'members' file and creates an exploded DataFrame (one row per Protein ID).

    Args:
        filename: Path for the eggNOG members file (e.g., '33208_members.tsv').

    Returns:
        A DataFrame with columns: ['OG_id', 'protein_ids', 'species_taxids'].
    """

    # Read the file and create a pandas dataframe
    logging.info("Reading and parsing file %s", filename)

    # Define columms we need to keep. Do not include other columns to save memory
    keep_cols = ["OG_id", "protein_ids"]

    try:
        members_df = pd.read_csv(
            Path(filename),
            sep="\t",
            header=None,
            names=[
                "tax_level",
                "OG_id",
                "n_proteins",
                "n_species",
                "protein_ids",
                "species_taxids",
            ],
            usecols=keep_cols,
        )

        # Split the protein_id column into separate rows
        # This can be done using the pandas explode function, but first we need to split the string of protein_ids into a list
        members_df["protein_ids"] = members_df["protein_ids"].str.split(",")
        # Now we can explode the dataframe to get one row per protein_id
        members_df = members_df.explode("protein_ids").reset_index(
            drop=True
        )  # reset the index after exploding
        # split the protein_id at the first dot to get the species_taxids and save it in a new column
        members_df["species_taxids"] = (
            members_df["protein_ids"].str.split(".", n=1).str[0]
        )
        logging.info(
            "Successfully parsed %s file with %d protein_id entries.",
            filename,
            len(members_df),
        )
        return members_df

    except FileNotFoundError:
        logging.error("%s file not found. Check file path", filename)
    except Exception as e:
        logging.error("An error occurred while reading the file %s: %s", filename, e)
        # return empty dataframe in case of failure
    return pd.DataFrame()


def get_ogs(taxid: str, df: pd.DataFrame) -> set[str]:
    """
    Returns a set of OGs that belong to given taxon ID.

    Args:
        taxid (str): The taxon ID for which to retrieve OGs.
        df (pd.DataFrame): The dataframe containing OG information, with columns "OG_id", "protein_ids", and "species_taxids".

    Returns:
        set[str]: A set of OGs belonging to the specified taxon ID.
    """

    logging.info("Getting OGs for taxid %s ...", taxid)
    # filter the dataframe for the given taxid and get the unique OGs
    # logic:
    #   df["species_taxids"] == taxid gives us a boolean series of the species_taxids column,
    #   which we can use to filter the dataframe.
    #   The outer filter [] filters only the rows where the condition is true.
    #   Then we can get the "OG_id" column and get the unique values as a set.

    try:
        og_set = set(df[df["species_taxids"] == str(taxid)]["OG_id"])
        logging.info("Finished. Found %s OGs belonging to taxid.", len(og_set))
    except KeyError:
        logging.error("The dataframe does not contain the expected columns.")
    except Exception as e:
        logging.error("An error occurred while getting OGs: %s", e)

    return og_set  # returns a set of strings containing the OGs belonging to taxid


# new get_sequences function for exploded dataset (BL)
def get_sequences(taxid: str, df: pd.DataFrame, ogs: set[str] | str = None) -> set[str]:
    """
    Returns sequences that belong to given taxon. If ogs is given,
    then returns only sequences that belong to these OGs.

    Args:
        taxid (str): The taxon ID for which to retrieve sequences.
        df (pd.DataFrame): The dataframe containing sequence information.
        ogs (set[str] | str, optional): A single OG string or a set of OGs to filter the sequences by.

    Returns:
        set[str]: A set of sequences belonging to the specified taxon.
    """
    logging.info("Getting sequences for taxid %s ...", taxid)

    try:
        is_taxid = df["species_taxids"] == str(taxid)
        taxid_data = df[is_taxid]
        if ogs is not None:
            if isinstance(
                ogs, str
            ):  # if only one OG is given as a string, convert it to a set with one element
                ogs = {ogs}
            is_og = taxid_data["OG_id"].isin(
                ogs
            )  # filter df for the given ogs=['...', '...']
            taxid_data = taxid_data[
                is_og
            ]  # filter the dataframe to only include rows that belong to the specified OG
        taxid_sequences_set = set(taxid_data["protein_ids"])
        logging.info(
            "Finished. Found %s sequences belonging to taxid %s.",
            len(taxid_sequences_set),
            taxid,
        )
    except KeyError:
        logging.error("The dataframe does not contain the expected columns.")
    except Exception as e:
        logging.error("An error occurred while getting sequences: %s", e)

    # returns a set of strings containing the sequences belonging to taxid (and if specified only of the given OG)
    return taxid_sequences_set


def _functional_category_dict(
    function_category_file: Path,
) -> dict:
    """
    Creates a mapping of functional category letters to their full function descriptions.

    Example:
        'J' -> 'Translation, ribosomal structure and biogenesis'

    Args:
        function_category_file: Path to the functional categories definition file.

    Returns:
        A dictionary with single-letters as keys and corresponding category names as values.
    """
    categories = {}
    # headers = None  # to keep track of current header

    with open(function_category_file, "r", encoding="utf-8") as file:
        lines = [line.strip() for line in file]
        for l in lines:
            if l.isupper():
                continue

            elif l.startswith("["):
                letter = l[1]  # between square brackets
                categories[letter] = l[4:]  # starts at 4th index

    logging.info(
        "The dictionary of the functional categories has been created: %s", categories
    )

    return categories


def annotate_dataframe(
    members: pd.DataFrame,
    protein_ids: set[str],
    annotationfile: Path,
    func_cat_file: Path,
) -> pd.DataFrame:
    """
    Annotates a subset of proteins ids with functional categories using eggNOG annotations.

    Args:
        members: The exploded DataFrame containing ['OG_id', 'protein_ids'].
        protein_ids: A set of protein ID strings to be annotated.
        annotationfile: Path to the eggNOG annotations (e.g., 'data/33208_annotations.tsv').

    Returns:
        A DataFrame containing merged protein IDs, OG IDs, and descriptive categories.
    """

    logging.info("Reading annotation file %s", annotationfile)
    try:
        annotation_df = pd.read_csv(
            str(annotationfile),
            sep="\t",
            usecols=[1, 2],  # just getting those columns we need
            names=[
                "OG_id",
                "single_letter_code",
            ],
        )
    except FileNotFoundError:
        logging.error("%s file not found", annotationfile)
    except Exception as e:
        logging.error("An error occurred while reading the file: %s", e)

    # merge annotation file with protein sequences in specific taxa to get just the relevant file
    members_filtered = members[members["protein_ids"].isin(protein_ids)]

    logging.info(
        "Filtered %s dimension out of %d proteins",
        members_filtered.shape,
        len(protein_ids),
    )

    annotation_df_merge = pd.merge(
        annotation_df,
        members_filtered[["OG_id", "protein_ids"]],
        on="OG_id",
        how="inner",
    )

    logging.info("Dimension of Filtered Annotation file: %s", annotation_df_merge.shape)

    categories = _functional_category_dict(func_cat_file)

    annotation_df_merge["functional_category"] = annotation_df_merge[
        "single_letter_code"
    ].map(lambda x: ", ".join(categories[char] for char in str(x)))
    logging.info("Finished merging.")

    return annotation_df_merge


def write_to_file(output_file_path: Path, data: set[str], description: str = None):
    """
    Write data (iterable) to a text file with specified description.

    Args:
        output_file_path: Path to the output file, e.g.: "results/2_orthologs_all.txt"
        data: iterable object that is written to the file
        description: description of the file, default = None
    """

    try:
        with open(output_file_path, "w", encoding="utf-8") as f:
            f.write(f"{description}\n")
            for item in data:
                f.write(f"{item}\n")
        logging.info("Written file for: %s", description)
    except Exception as e:
        logging.error("Error writing to file %s: %s", output_file_path, e)
