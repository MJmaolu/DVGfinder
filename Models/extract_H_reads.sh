#!/bin/bash
# extract_H_reads.sh

# Extrae del fichero sam de alineamiento las reads que contienen H en su CIGAR

complete_sam=$1
reads_H=$2

awk '$1 ~ /^@/ || $6 ~/H/ {print $0}' $complete_sam > $reads_H