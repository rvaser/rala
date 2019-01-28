#!/usr/bin/env bash

dataset=$(basename $1)

minimap=
#minimap2=
rala=

# preproces
$minimap -t12 -L100 -Sw5 -m0 $1 $1 > "$dataset"_overlaps.paf
#$minimap2 -t12 -x ava-pb $1 $1 > "$dataset"_overlaps.paf

$rala -t12 $1 "$dataset"_overlaps.paf -p > "$dataset"_preconstructed.fasta

# assembly
$minimap -t12 -L100 -w5 -m0 -f0.00001 "$dataset"_preconstructed.fasta $1 > "$dataset"_sensitive_overlaps.paf
#$minimap2 -t12 -x ava-pb --dual=yes -f 0.00001 "$dataset"_preconstructed.fasta $1 > "$dataset"_sensitive_overlaps.paf

$rala -t12 $1 "$dataset"_overlaps.paf -s "$dataset"_sensitive_overlaps.paf > "$dataset"_layout.fasta

# cleanup
rm "$dataset"_overlaps.paf
rm "$dataset"_sensitive_overlaps.paf
rm "$dataset"_preconstructed.fasta
