#!/usr/bin/env bash

dataset=$(basename $1)
minimap=
rala=
rast=
gepard=
matrix=

# preproces
$minimap -t12 -L100 -Sw5 -m0 $1 $1 > "$dataset"_overlaps.paf

#$minimap -t12 -x ava-pb $1 $1 > "$dataset"_overlaps.paf

$rala -t12 -d "$dataset"_debug $1 "$dataset"_overlaps.paf > "$dataset"_layout.fasta


# assembly
$minimap -t12 -L100 -w5 -m0 -f0 "$dataset"_debug.fasta $1 > "$dataset"_all_overlaps.paf

#$minimap -t12 -x ava-pb --dual=yes -f 0.0 "$dataset"_debug.fasta $1 > "$dataset"_all_overlaps.paf

$rala -t12 -d "$dataset"_debug_final -r "$dataset"_all_overlaps.paf $1 "$dataset"_overlaps.paf > "$dataset"_layout_final.fasta


# statistics
$rast -r $2 "$dataset"_layout_final.fasta 2>> rast_statistics.txt

java -cp $gepard org.gepard.client.cmdline.CommandLine -seq1 $2 -seq2 "$dataset"_layout_final.fasta -matrix $matrix -outfile "$dataset"_debug.png -lower 0 -upper 100

# cleanup
rm "$dataset"_overlaps.paf
rm "$dataset"_all_overlaps.paf
rm "$dataset"_layout.fasta
