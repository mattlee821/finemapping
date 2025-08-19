#!/bin/bash
VAR_study="ferkingstad_2021_PMID34857953"
FILE_IN="/data/GWAS_data/work/${VAR_study}/window/"
FILE_OUT="/data/IARC_Biostat/work/finemapping/src/${VAR_study}/"

mkdir -p "$FILE_OUT"
mkdir -p "/data/IARC_Biostat/work/finemapping/analysis/${VAR_study}/"

find "$FILE_IN" -type f > "${FILE_OUT}files.txt"
wc -l < "${FILE_OUT}files.txt"
