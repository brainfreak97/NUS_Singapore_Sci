#!/bin/bash

#Assuming proper installation of gromacs through sudo apt-get install gromacs
#Assuming all codes are in the Desktop domain
#Required PDB files: APOBEC-XTEN-Cas13b.pdb
#Other available PDB files: APOBEC-Linker1-Cas13b.pdb, APOBEC-Linker2-Cas13b.pdb, APOBEC-Linker3-Cas13b.pdb
#Required mdp files for simulation: ions.mdp, minim.mdp, nvt.mdp, npt.mdp, md.mdp
#Code should generate MD_FP.tpr at the end of the run

#Generating topology files
grep -v HOH APOBEC-XTEN-Cas13b.pdb > Test1.pdb
cat <<EOF1 > inp1
8 #selection for CHARMM27 force field
EOF1
gmx pdb2gmx -f Test1.pdb -o Test1.gro -water spce <inp1

#Creating a solvation box
gmx editconf -f Test1.gro -o Test1_boxed.gro -c -d 1.0 -bt cubic 
#Change from cubix to dodecagon if required
gmx solvate -cp Test1_boxed.gro -cs spc216.gro -o Test1_solvated.gro -p topol.top
#Delete below if no debug is required
gmx editconf -f Test1_solvated.gro -o solvated.pdb

#Adding ions to the box
gmx grompp -f ions.mdp -c Test1_solvated.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o Test1_solvated_ion.gro -p topol.top -pname NA -nname CL -neutral 

cat<<EOF2 > inp2
13 #Selecting for SOL option here. Change to the SOL option if the num assigned is different
EOF2

#Delete below if no debug is required
gmx editconf -f Test1_solvated_ion.gro -o solvated_ionized.pdb < inp2 #next select for the regions to replace the mocules with ions


#Energy minimization
gmx grompp -f minim.mdp -c Test1_solvated_ion.gro -p topol.top -o test1_em.tpr
gmx mdrun -s test1_em.tpr -v -deffnm em

#To check for minimization output
cat<<EOF6 > inp6
12 #Selecting for potential option here.
EOF6
gmx energy -f em.edr -o potential.xvg <inp6
xmgrace potential.xvg

#Cancel below if not interested to get the PDB files
gmx editconf -f Test1_solvated_ion.gro -o solvated_ionized_minid.pdb

#NVT Equilibriation
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

cat<<EOF3 > inp3
15 #selecting temp, choose 37 degree Celcius to mimic within cell conditions
EOF3

#To check for NVT equilibration output
gmx energy -f nvt.edr -o nvt_img.xvg <inp3
xmgrace nvt_img.xvg

#NPT Equilbration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

cat<<EOF4 > inp4
18 #selecting pressure
EOF4

#To check for NPT equilibration output
gmx energy -f npt.edr -o npt_img.xvg<inp4
xmgrace npt_img.xvg

#Checking density

cat<<EOF5 > inp5
24
EOF5

gmx energy -f npt.edr -o density_img.xvg<inp5
xmgrace density_img.xvg


#If all else passed
#Finally, MD simulations
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o MD_FP.tpr
# Check estimate for relative computational load of PME mesh particle
#Value should be around 0.25
#Selection of core as an automatic process
gmx mdrun -v -deffnm MD_FP



