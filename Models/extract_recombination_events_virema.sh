#!/bin/bash

# extract_recombination_events_virema.sh

# This script extract the recombination events from the virema output file
# '{sample_name}_Virus_Recombinations_Results.txt'. Split the info in lines
# with the info: sense, BP, RI, read_counts_virema

sample_name=$1
ref=$2
recomb_file=$3

# create directory
mkdir dir_temp

# senses
for_for="@NewLibrary: ${ref}_to_${ref}"
rev_rev="@NewLibrary: ${ref}_RevStrand_to_${ref}_RevStrand"
for_rev="@NewLibrary: ${ref}_to_${ref}_RevStrand"
rev_for="@NewLibrary: ${ref}_RevStrand_to_${ref}"

# Extract each event in a different file based on its sense

cat $recomb_file | grep -w -A1 -- "$for_for" | sed $'s/\t/\\\n/g' | sed $'s/_to_/\t/g' | sed $'s/_#_/\t/g' > dir_temp/${sample_name}_temp_events_++.tsv
cat $recomb_file | grep -w -A1 -- "$rev_rev" | sed $'s/\t/\\\n/g' | sed $'s/_to_/\t/g' | sed $'s/_#_/\t/g' > dir_temp/${sample_name}_temp_events_--.tsv
cat $recomb_file | grep -w -A1 -- "$for_rev" | sed $'s/\t/\\\n/g' | sed $'s/_to_/\t/g' | sed $'s/_#_/\t/g' > dir_temp/${sample_name}_temp_events_+-.tsv
cat $recomb_file | grep -w -A1 -- "$rev_for" | sed $'s/\t/\\\n/g' | sed $'s/_to_/\t/g' | sed $'s/_#_/\t/g' > dir_temp/${sample_name}_temp_events_-+.tsv

# For events file: add a new column with the sense of the event and clean last line

awk '$1!~/^@/ {print $0 "\t++"}' dir_temp/${sample_name}_temp_events_++.tsv > dir_temp/${sample_name}_temp_sense_events_++.tsv 
awk '$0!~/^\t/ {print $0}' dir_temp/${sample_name}_temp_sense_events_++.tsv > dir_temp/${sample_name}_temp_sense2_events_++.tsv 

awk '$1!~/^@/ {print $0 "\t--"}' dir_temp/${sample_name}_temp_events_--.tsv > dir_temp/${sample_name}_temp_sense_events_--.tsv 
awk '$0!~/^\t/ {print $0}' dir_temp/${sample_name}_temp_sense_events_--.tsv > dir_temp/${sample_name}_temp_sense2_events_--.tsv 

awk '$1!~/^@/ {print $0 "\t+-"}' dir_temp/${sample_name}_temp_events_+-.tsv > dir_temp/${sample_name}_temp_sense_events_+-.tsv 
awk '$0!~/^\t/ {print $0}' dir_temp/${sample_name}_temp_sense_events_+-.tsv > dir_temp/${sample_name}_temp_sense2_events_+-.tsv 

awk '$1!~/^@/ {print $0 "\t-+"}' dir_temp/${sample_name}_temp_events_-+.tsv > dir_temp/${sample_name}_temp_sense_events_-+.tsv 
awk '$0!~/^\t/ {print $0}' dir_temp/${sample_name}_temp_sense_events_-+.tsv > dir_temp/${sample_name}_temp_sense2_events_-+.tsv 


# Join in all the recombination events in a single file
cat dir_temp/${sample_name}_temp_sense2_events_++.tsv dir_temp/${sample_name}_temp_sense2_events_--.tsv dir_temp/${sample_name}_temp_sense2_events_+-.tsv dir_temp/${sample_name}_temp_sense2_events_-+.tsv > Outputs/${sample_name}_from_raw_virema.tsv

# Remove temporary files
rm -r dir_temp/	
