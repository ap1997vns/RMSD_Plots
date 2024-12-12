#!/bin/bash
# get all RMSD plots together and using awk to enter extra columns so that i can plot from R 
rm -f ALLDomainRMSD.txt
mutants=("WT" "A97V" "R102Q" "S171L" "T180I" "M188I" "E152A" "R71C")
base_dir="/home/user/WORK/DgoR/APO"
temp_dir="temp_rmsd_files"
mkdir -p "$temp_dir"

for mut in "${mutants[@]}"
do
    dirname="$base_dir/$mut/charmm-gui-7169740005/gromacs"
    awk -v val=APO$mut '{print val",Protein,"$1*0.002","$2}' "$dirname/$mut-APOproteinBBrmsd.txt" >> "$temp_dir/$mut-rmsd.txt"
    awk -v val=APO$mut '{print val",Cterm,"$1*0.002","$2}' "$dirname/$mut-APOCterm-proteinBBrmsd.txt" >> "$temp_dir/$mut-rmsd.txt"
    awk -v val=APO$mut '{print val",Nterm,"$1*0.002","$2}' "$dirname/$mut-APONterm-proteinBBrmsd.txt" >> "$temp_dir/$mut-rmsd.txt"
    awk -v val=APO$mut '{print val",Nucleic,"$1*0.002","$2}' "$dirname/$mut-APOnucleicBBrmsd.txt" >> "$temp_dir/$mut-rmsd.txt"
done

cat "$temp_dir"/*.txt > ALLDomainRMSD.txt
rm -rf "$temp_dir"
echo "RMSD values from all mutants have been consolidated into ALLDomainRMSD.txt"
