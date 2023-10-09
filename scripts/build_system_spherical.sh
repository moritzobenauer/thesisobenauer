#!/bin/bash
# Moritz Lennart Obenauer, JGU Mainz, August 2023
#
# Input Variables
# Example. ./build_system.sh 10 15 single.gro
mol=$1
box=$2
filename=$3

for i in $(seq 1 $mol)

# This translational operation leads only to spherical systems!

do
   translate=$((5*$i))
   rotate1=$((0 + 60*$i))
   rotate2=$((0 - 40*$i))
   rotate3=$((0 - 30*$i))
   gmx editconf -f $filename -o single$i.pdb -rotate $rotate1 $rotate2 $rotate3 -translate 5 5 5
   tail -n +5 single$i.pdb > temp_$i.pdb
   head -n -2 temp_$i.pdb > temp_2_$i.pdb
   echo $j
done

# PDB File Adjustments

head -n -2 single0.pdb > system0.pdb

for i in $(seq 0 $(($mol - 1)))
do
   j=$(($i + 1))
   cat system$i.pdb temp_2_$j.pdb > system$j.pdb
   echo $i
   echo $j
   end=$j
done


# Add TER ENDMDL

printf "TER" >> system$j.pdb
printf "ENDMDL" >> system$j.pdb

# Output GRO Conversion

gmx editconf -f system$j.pdb -o system.gro -box $box $box $box -bt cubic

rm *.pdb
