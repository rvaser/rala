#!/usr/bin/env bash

dataset=$(basename $1)
minimap=
rala=

# preproces
$minimap -t12 -x ava-ont $1 $1 > "$dataset"_overlaps.paf

$rala -t12 -p $1 "$dataset"_overlaps.paf > "$dataset"_preconstructed.fasta


# assembly
$minimap -t12 -x ava-ont --dual=yes -f 0.0001 "$dataset"_preconstructed.fasta $1 > "$dataset"_sensitive_overlaps.paf

$rala -t12 -s "$dataset"_sensitive_overlaps.paf $1 "$dataset"_overlaps.paf > "$dataset"_layout_final.fasta

# cleanup
rm "$dataset"_overlaps.paf
rm "$dataset"_sensitive_overlaps.paf
rm "$dataset"_preconstructed.fasta
rm "$dataset"_layout_final.fasta
