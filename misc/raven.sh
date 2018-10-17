#!/usr/bin/env bash

dataset=$(basename $1)
minimap=
rala=
rast=
gepard=
matrix=

$minimap -t12 -L100 -Sw5 -m0 $1 $1 > "$dataset"_overlaps.paf

$rala -t12 -d "$dataset"_debug $1 "$dataset"_overlaps.paf > "$dataset"_layout.fasta

$minimap -t12 -L100 -w5 -m0 "$dataset"_debug_knots.fasta $1 > "$dataset"_all.paf

$rala -t12 -d "$dataset"_debug_after -r "$dataset"_all.paf $1 "$dataset"_overlaps.paf > "$dataset"_layout_final.fasta

$rast -r $2 "$dataset"_layout_final.fasta 2>> rast_statistics.txt

java -cp $gepard org.gepard.client.cmdline.CommandLine -seq1 $2 -seq2 "$dataset"_layout_final.fasta -matrix $matrix -outfile "$dataset"_gepard.png -lower 0 -upper 100

rm "$dataset"_overlaps.paf
rm "$dataset"_all.paf
rm "$dataset"_layout.fasta
rm "$dataset"_layout_final.fasta
