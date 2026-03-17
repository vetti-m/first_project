#!/bin/bash
mkdir -p results
mkdir -p tmp
mkdir -p data


# Download the data files from eggNOG and place them in the data directory (gz to be faster)
cd data

if [ -e ./33208_members.tsv ]; then
    echo "File exists."
else
    wget http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/33208/33208_members.tsv.gz && gunzip 33208_members.tsv.gz
fi

if [ -e ./33208_annotations.tsv ]; then
    echo "File exists."
else
    wget http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/33208/33208_annotations.tsv.gz && gunzip 33208_annotations.tsv.gz
fi

if [ -e ./e5.taxid_info.tsv ]; then
    echo "File exists."
else
    wget http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv
fi

if [ -e ./eggnog4.functional_categories.txt ]; then
    echo "File exists."
else
    wget http://eggnog5.embl.de/download/eggnog_4.5/eggnog4.functional_categories.txt
fi

# Get the Primate Scientific Names of eggNOG
cd ..
grep Primates data/e5.taxid_info.tsv | cut -f 2 > tmp/primates.txt