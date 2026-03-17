#!/usr/bin/env python
import logging
from pathlib import Path
import pandas as pd
import eggnog5 as e5


# logging.basicConfig?

logging.basicConfig(
    filename="doc/metazoa_pipeline.log",
    filemode="w",  # We want to overwrite each time!
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",  # Realpython
    force=True,
)


# Setting variables -> relative path, maybe need to adapt, depends from where the runallscript is running!
METAZOA_GENES = Path("data/33208_members.tsv")
TAX_ID = Path("data/e5.taxid_info.tsv")
ANNOTATION = Path("data/33208_annotations.tsv")
FUNC_CATEGORY = Path("data/eggnog4.functional_categories.txt")

# Get necessary TaxID
SPECIES_TAXID = e5.retrieve_taxid(TAX_ID)
PRIMATES = set()
with open("tmp/primates.txt", "r", encoding="utf-8") as file:
    lines = [line.strip() for line in file]
    for line in lines:
        PRIMATES.add(line)


# Parse the members.tsv file
MEMBERS_DF = e5.parse_members_file(METAZOA_GENES)


def main():
    logging.info("Starting the main function.")

    # ----------------------------------------------------------------------------------------------------------
    # 1.A:
    # get OGs of human, chimp and mouse in separate sets
    human_ogs = e5.get_ogs(SPECIES_TAXID["Homo sapiens"], MEMBERS_DF)
    chimp_ogs = e5.get_ogs(SPECIES_TAXID["Pan troglodytes"], MEMBERS_DF)
    mouse_ogs = e5.get_ogs(SPECIES_TAXID["Mus musculus"], MEMBERS_DF)

    # get OGs in human and chimps but not in mice and store in a set
    human_chimp_not_mouse = human_ogs & chimp_ogs - mouse_ogs

    # store human sequences of the previously filtered OG ids
    human_seq = e5.get_sequences(
        SPECIES_TAXID["Homo sapiens"], MEMBERS_DF, human_chimp_not_mouse
    )

    result_A = len(
        human_seq
    )  # this is the number of human genes that are homolog to chimps but not mice
    logging.info(
        "Found %d human genes with chimp homologs but no mouse homologs.", result_A
    )

    # ----------------------------------------------------------------------------------------------------------
    # 1.B:
    # What are their protein IDs? Store them in a result file
    output_path_1A = Path("results/1b_human_genes_chimp_not_mouse.txt")

    logging.info(
        "Starting to write %d protein ids to %s", len(human_seq), output_path_1A
    )

    e5.write_to_file(
        output_path_1A,
        human_seq,
        "Genes that are found in humans and chimps, but not in mice",
    )

    # ----------------------------------------------------------------------------------------------------------
    # 1.C:
    annotation_df = e5.annotate_dataframe(
        MEMBERS_DF, set(human_seq), ANNOTATION, FUNC_CATEGORY
    )
    output_path_1C = Path("results/gene_counts.csv")
    counter = annotation_df["functional_category"].value_counts()
    counter.to_csv("results/1c_gene_counts_functional_categories.csv", sep="\t")
    logging.info(
        "Created a file with all Functional groups and the gene count of them for Genes in Human and Chimp but not Mouse: %s",
        output_path_1C,
    )

    # ---------------------------------------------------------------------
    # 1.D:
    # Are these genes only found in human and chimp, or do some (or all?)
    # of them have homologs in other metazoan species?

    # make a new df containing only our filtered OG ids
    subset_MEMBERS_DF = MEMBERS_DF[MEMBERS_DF["OG_id"].isin(human_chimp_not_mouse)]

    # count the number of unique species taxids for each OG_id in the subset_MEMBERS_DF
    og_species_counts = subset_MEMBERS_DF.groupby("OG_id")["species_taxids"].nunique()

    # if the number of species is bigger than 2, then we keep the OG_id, and put the species in a list
    # example: OG_id: [species1, species2, species3]
    ogs_over2species = og_species_counts[og_species_counts > 2].index.tolist()

    logging.info(
        "Found %d OGs with more than 2 species (not only homologs in human and chimp)",
        len(ogs_over2species),
    )  # 1052
    logging.info("Originally we had %d OGs", len(human_chimp_not_mouse))  # 1073

    # get the OGs that only have human and chimp proteins, but no other species
    og_only_human_chimp = set(human_chimp_not_mouse) - set(ogs_over2species)
    logging.info("OG IDs only in Human and Chimpanzees: %s", og_only_human_chimp)

    # human proteins that only have chimp homologs but no other species homologs
    e5.get_sequences(SPECIES_TAXID["Homo sapiens"], MEMBERS_DF, og_only_human_chimp)

    logging.info(
        "Out of %d human genes with chimp homologs, only %d are unique to these two species.",
        len(human_seq),
        len(
            e5.get_sequences(
                SPECIES_TAXID["Homo sapiens"], MEMBERS_DF, og_only_human_chimp
            )
        ),
    )

    # --------------------------------------------------------------------------------
    # 1.E:
    # How many of these genes are primate-specific, i.e. found only in primates?
    # Get set of primate taxIDs

    primate_taxids = set()
    for species in PRIMATES:
        primate_taxids.add(int(SPECIES_TAXID[species]))

    # get set of non-primate taxIDs after subtracting from total number of taxids
    species_taxid_total = set(MEMBERS_DF["species_taxids"].astype(int))  # 161
    not_primates = species_taxid_total - primate_taxids
    logging.info(
        "Of %d Total species, we have %d non-primate species.",
        len(species_taxid_total),
        len(not_primates),
    )  # 146

    # look for rows in df that have non-primates
    rows_not_primates = MEMBERS_DF["species_taxids"].astype(int).isin(not_primates)

    # get the OG group of the non-primates using .loc and make a set out of it
    og_not_primates = set(MEMBERS_DF.loc[rows_not_primates, "OG_id"])

    # get the OG groups of all and subtract the non-primate OG groups
    all_og_ids = set(MEMBERS_DF["OG_id"])
    og_only_primates = all_og_ids - og_not_primates

    logging.info(
        "There are %d OGs primate-specific.", len(og_only_primates)
    )  # 1545 OG ids

    # get human genes that only have homologs in primates:
    human_genes_only_primates = e5.get_sequences(
        SPECIES_TAXID["Homo sapiens"], MEMBERS_DF, og_only_primates
    )
    logging.info(
        "There are %d human genes primate-specific.", len(human_genes_only_primates)
    )  # 332 human genes

    # ----------------------------------------------------------------------------------------------------------
    # 2:
    # Identify the orthologs that are present in primates, chicken and fish
    # get OGs for the fish species, human and chimp OGs from 1A
    danio_ogs = e5.get_ogs(SPECIES_TAXID["Danio rerio"], MEMBERS_DF)
    takifugu_ogs = e5.get_ogs(SPECIES_TAXID["Takifugu rubripes"], MEMBERS_DF)

    # Define groups
    primates = human_ogs | chimp_ogs  # union
    chicken = e5.get_ogs(SPECIES_TAXID["Gallus gallus"], MEMBERS_DF)
    fish = danio_ogs | takifugu_ogs

    # get set of OGs common in all groups
    orthologs = primates & chicken & fish  # intersection

    # Are there OGs that were lost in mouse and rat or only in mouse / only in rat

    # get OGs for rat
    rat_ogs = e5.get_ogs(SPECIES_TAXID["Rattus norvegicus"], MEMBERS_DF)

    lost_mouse = orthologs - mouse_ogs  # OGs lost in mouse
    lost_rat = orthologs - rat_ogs  # OGs lost in rat
    lost_mouse_rat = lost_mouse & lost_rat  # OGs that are lost in both
    lost_mouse_only = lost_mouse - lost_rat  # OGs that are lost in mouse only
    lost_rat_only = lost_rat - lost_mouse  # OGs that are lost in rat only

    # sanity check
    # print(len(lost_mouse_only) + len(lost_mouse_rat) == len(lost_mouse)) #True
    # print(len(lost_rat_only) + len(lost_mouse_rat) == len(lost_rat)) #True

    logging.info("Writing results to output files...")

    output_path_2_1 = Path("results/2_orthologs_all.txt")
    output_path_2_2 = Path("results/2_orthologs_lost_mouse_rat.txt")
    output_path_2_3 = Path("results/2_orthologs_lost_mouse_only.txt")
    output_path_2_4 = Path("results/2_orthologs_lost_rat_only.txt")

    e5.write_to_file(output_path_2_1, orthologs, "Orthologs found in all groups.")
    e5.write_to_file(
        output_path_2_2, lost_mouse_rat, "Orthologs lost in both mouse and rat."
    )
    e5.write_to_file(output_path_2_3, lost_mouse_only, "Orthologs lost in mouse only.")
    e5.write_to_file(output_path_2_4, lost_rat_only, "Orthologs lost in rat only.")

    # ----------------------------------------------------------------------------------------------------------
    # 3:
    # Get 99% of species -> threshold
    import math  # for rounding up

    nspecies_99 = math.ceil(len(SPECIES_TAXID) * 0.99)

    # just take unique species/og combinations
    uniques = MEMBERS_DF[["OG_id", "species_taxids"]].drop_duplicates()

    # count how many uniques are and safe it in a data frame (so next steps will be easier)
    count_ogs = uniques["OG_id"].value_counts().to_frame(name="count")

    # combine output with threshold and safe in file
    metazoa_ogs = count_ogs[count_ogs["count"] >= nspecies_99]
    metazoa_ogs.to_csv("results/3_metazoa_universal_ogs.csv", sep="\t")

    logging.info(
        "There have been %d values clipped, so %d ortholog groups are in >= 99 percent of metazoa with max of %d of %d metazoa species.",
        len(count_ogs) - len(metazoa_ogs),
        len(metazoa_ogs),
        count_ogs["count"].max(),
        len(SPECIES_TAXID),
    )


if __name__ == "__main__":
    main()
