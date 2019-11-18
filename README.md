# PepGen
Randomised peptide conformation generator from known docked positions

## Requirements
1. PepGen only works on Linux OS due to dependency on the MAFFT alignment software
2. Requires R v3.6.1

## How to run PepGen
1. Create folder with mhc name in /Inputfiles directory
2. Place prob.txt (from I-Tasser/COACH) and your MHC .fasta file into the folder
3. Run the following code from the command line within the /PepGen directory
    Rscript --vanilla Extended_Pocket_Quantification.R <MHC/folder_name> <Pocket B/F>




