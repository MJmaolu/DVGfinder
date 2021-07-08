#!/bin/bash

# extract_recombination_events_ditector.sh

sample_name=$1
recombinations_file=$2

mkdir dir_temp_ditector

# Eliminamos las separaciones
sed $'s/ //g' $recombinations_file | sed $'s/\'//g' |
awk  -v id=$sample_name '$1~/^3|^5|^Del|^Ins/ {print $1, $3, $4, $7}' >> dir_temp_ditector/${sample_name}_temp_events.txt

# Reducimos la terminología de los DVG_type:
#   Deletion(Fwd.strand) --> Deletion_forward
#   Deletion(Rev.strand) --> Deletion_reverse
#   Insertion(Fwd.strand) --> Insertion_forward
#   Insertion(Rev.strand) --> Insertion_reverse

sed 's/DVG//' dir_temp_ditector/${sample_name}_temp_events.txt | 
sed $'s/Deletion(Fwd.strand)/Deletion_forward/g' |
sed $'s/Deletion(Rev.strand)/Deletion_reverse/g' |
sed $'s/Insertion(Fwd.strand)/Insertion_forward/g' |
sed $'s/Insertion(Rev.strand)/Insertion_reverse/g' > dir_temp_ditector/${sample_name}_temp2_events.txt

# Format to tsv file and remove "DVG" from the text
sed $'s/ /\t/g' dir_temp_ditector/${sample_name}_temp2_events.txt > Outputs/${sample_name}_from_raw_ditector.tsv

# Borramos carpeta de archivos temporales
rm -r dir_temp_ditector