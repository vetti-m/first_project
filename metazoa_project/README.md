Metazoa_pipeline 
==============
ABOUT
======
This project aims to study metazoan genes in different lineages to answer questions about gene conservation and gene loss. It uses the eggnog5 database as a resource to retrieve files containing information about orthologous groups in metazoa with corresponding proteinIDs and species TaxIDs and the annotation of the function of genes belonging to a specific orthologous group. 

PREREQUESITES
==============
-	Linux operating systems
-	Python3
-	Pandas module

INSTALLATION
============
Download metazoa_pipeline.tar.gz
This archive contains all essential scripts to run this program:
runall.sh, prepare_pipeline.sh, eggnog5.py, metazoa_pipeline.py
Execute the bash script runall.sh. This executes helper scripts to prepare folder structures and download necessary files from the eggnog5 database, activates the base conda environment, and executes the metazoa_pipeline.py script.
$ bash runall.sh

FEATURES
=========
Files from eggnog5 used for this project: 
--------------------------------------------------
e5.taxid_info.tsv (Tax ID | Scientific name | Rank | Named lineage | TaxID lineage)
33208_annotations.tsv (Taxon level ID | Orthologous group ID | COG functional category | Functional description)
33208_members.tsv (Taxon level ID | Orthologous group ID | Protein Count | Species Count | Protein ID | Species TaxID)
eggnog4.funtional_categories.txt (COG functional category | Description)
For more detailed description see http://eggnog5.embl.de/
- To extract the TaxID of a given species out of the e5.taxid_info.tsv, a function (retrieve_taxid) was written to create a dictionary mapping the scientific name to its corresponding TaxID 
- To parse the members.tsv file to efficiently extract needed information, a function (parse_ members_file) was written to create a pandas dataframe containing information on the OG-ID, protein-ID and TaxID for each Protein-ID in the file (one row per proteinID)
Additionally, functions were created to extract information out of this parsed dataframe: 
-	“get_OGs”: function to extract OG-IDs belonging to a given taxon
-	“get_sequences”: function to extract protein-IDs belonging to a given taxon, filtered on specified OG-IDs were written.
Questions answered in this project:
---------------------------------------------
1.	To study what makes primates like humans and chimps special compared to other mammals like mice, identify all human (Homo sapiens) genes that have at least one homolog in chimp (Pan troglodytes), but don't have a homolog in mouse (Mus musculus).
-	First sets of orthologous group IDs (OG IDs) for human (Homo sapiens), chimp (Pan troglodytes) and mouse (Mus musculus) were created using the previously created function get_ogs. Then, using set operations (intersection and difference) a set was created containing only OG-IDs that contain both human and chimp genes but no mouse genes. Then sequences (ProteinIDs) belonging to the human TaxID were retrieved, additionally filtering on the previously created set of OGs (human_chimp_not_mouse). The resulting ProteinIDs were saved in a new set.
A.	How many such genes are there?
-	The length of the created set gives us the number of genes. We identified 1316 human genes that have homologs in chimpanzees but no homologs in mice. 
B.	What are their protein IDs? Store them in a result file.
-	The set of Protein IDs were written to the file “1b_human_genes_chimp_not_mouse.txt”
C.	It would be interesting to get an overview of their functions. If a functional description from eggNOG is available: Which functional categories do they have, how many genes are in each category? Store the result table (function -> number of genes) in a result file, sorted by the most common categories.
-	A dictionary with a single letter description as keys and the matching functional category as values gives the base. Next, a filtered data frame with only the protein IDs for the given species was created. After merging those with the original created data frame (created from the eggNOG annotation file), it was possible to map the functional categories by using the single letter as a common key. 
-	It is important to mention that some orthologous groups have multiple single letter codes and therefore more than 1 functional group is merged.
-	From 1613 genes, 666 have an unknown function, 190 have a transcriptional function, and 124 would be responsible for a signal transduction mechanism. There are more functional groups which are all saved (including the number of corresponding genes) in the 1c_gene_counts_functional_categories.csv-file. Many functional groups have only one corresponding gene. 
D.	Are these genes only found in human and chimp, or do some (or all?) of them have homologs in other metazoan species?
-	To answer the question if the genes found in 1B also have homologs in other metazoan species, a data frame that contained only the previously filtered OGs was created. Then, the number of unique species TaxIDs in this data frame were counted for each OG. To find out which orthologous groups contain more than two TaxIDs (these are human and chimp because we filtered on them), the OG IDs were filtered based on having a count >2 and saved as a list. To obtain the number of OGs that only contain homologs from human and chimps, the list was converted to a set and the containing OGs were removed from the human_chimp_not_mouse set. 
We identified 1052 OGs that also contain homologs of other metazoan species (more than 2 species in the group) out of the 1073 originally identified OGs. Out of 1316 human genes with chimp homologs, only 21 are unique to these two species.
E.	How many of these genes are primate-specific, i.e. found only in primates? (You can use the list of primate species below.)
-	To isolate primate-specific genes, first a set of non-primates was created by subtracting the set of primate species from the set containing all species. Then every OG that contained at least one non-primate species was excluded, so that only primate-specific OGs were left. Finally, these OGs were used to find the human sequences associated with the orthologous groups.
-	Out of the 57839 metazoan orthologous groups, only 1545 groups were defined as primate-specific. 332 human genes were found in these OGs.


2.	Identify the orthologs that are present in primates (in human, chimp, or both), chicken (Gallus gallus) and fish (Danio rerio, Takifugu rubripes, or both). Were any of these orthologous groups lost both in mouse and rat? Were some lost only in mouse, and some lost only in rat? (And if so, which were the groups lost in each species?)
-	Here we looked at orthologs that are widely conserved across primates, birds and fish, but lost in specific rodent species (mouse and rat). Sets containing the OGs of the considered species (human OR chimp, chicken & fish) were constructed. Then the intersection of these groups was stored as a set. Afterwards, this set was compared against the OGs present in mice and rats and orthologs lost in both species, lost only in rat and lost only in mouse were identified.
-	Out of 9861 orthologs that are present in primates, chicken and fish --> 108 orthologs were lost in both mouse and rat; 63 orthologs were lost only in mouse; 138 orthologs were lost only in rat. All identified orthologous groups were exported to the results/ directory and can be analyzed further.

3.	Can you identify orthologous groups which are universal to all animals, e.g. that occur in 99% or more of all animal species?
-	To identify these universal orthologous groups, a species number threshold was calculated by taking 99% of the total count of metazoa species of our members' dataset and rounding up to the nearest integer. Duplicates of OGs were removed and then we counted how many different species were associated with each OG group. If a group contained more than 99% of all used metazoa species, it was saved and counted together afterwards.
-	Out of the 57839 total orthologous groups, only 694 OGs were universal to all animals and therefore present in at least 160 out of 161 species. For a detailed view of these 694 universal orthologous groups, the file in /results/3_metazoa_universal_ogs.csv can be looked at.


OUTPUT FILES
============
saved in "data" Folder:
files downloaded from eggnog database:
33208_annotations.tsv
33208_members.tsv
e5.taxid_info.tsv
eggnog4.functional_categories.txt


Saved in “results” folder:
1b_human_genes_chimp_not_mouse.txt
1c_gene_counts_functional_categories.csv
2_orthologs_all.txt
2_orthologs_lost_mouse_rat.txt
2_orthologs_lost_mouse_only.txt
2_orthologs_lost_rat_only.txt
3_metazoa_universal_ogs.csv

Saved in “doc” folder:
metazoa_pipeline.log


CONTRIBUTORS
===============
Martina Gruber: Flowchart; basic configuration logging, runall.sh, prepare_pipeline.sh, function development: retrieve_taxid(), _functional_category_dict(), annotate_dataframe(); exercises 1A, 1C, 1E (getting primate taxIDs), 3, merging code snipplets and journal snipplets together.
Katharina Karner: function development: parse_members_file(), get_ogs(), write_to_file(); exercises 1A, 1C, 1E, 2
Bianca Lang: runall.sh, prepare_pipeline.sh; function development: get_ogs(), get_sequences(), write_to_file(); exercises 1A, 1B, 1D, 1E, 3

(ALL: General plan for Coding, Debugging, Journaling, README file,...)
