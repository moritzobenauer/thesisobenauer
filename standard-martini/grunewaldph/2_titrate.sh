#!/bin/bash 

#pHrange=(3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5 5.75 6.0 6.25 6.5 6.75 7.0 7.25 7.5 7.75 8.0)
pHrange=(3.5 4.5 5.5 6.5 7.5)
#pHrange=(3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75 7.25 7.75)
for pH in ${pHrange[*]}; do

  rm -rf ${pH}
  mkdir ${pH}
  cd ${pH}
  
  sed "s/<value>/${pH}/g" ../system.top > system.top 
  
  mkdir min
  mkdir eq
  mkdir NVT
  
  cd min 
  
  gmx grompp -f ../../../mdp_files/min.mdp -c ../../start.gro -p ../system.top -maxwarn 55 >> grompp.out 2>&1
  gmx mdrun -v
  
  cd ../eq
  
  gmx grompp -f ../../../mdp_files/eq.mdp -c ../min/confout.gro -p ../system.top -maxwarn 55 >> grompp.out 2>&1
  gmx mdrun -nsteps 100000 -v 
  
  cd ../NVT;
  gmx grompp -f ../../../mdp_files/NVT.mdp -c ../eq/confout.gro -p ../system.top -maxwarn 55;
  gmx mdrun -nsteps 10000000 -v -pin on
  
  cd ../..

done
