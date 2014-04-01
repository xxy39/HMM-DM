HMM-DM
======
The HMM-DM program identifies differentially methylated (DM) CG sites and regions from whole genome and targeted bisulfite sequencing (BS) data. This approach first uses a Hidden Markov Model to identify differentially methylated CG sites accounting for spatial correlation across CGs and variation across samples, and then summarizes identified DM CG sites into regions based on their status and distance. This program takes aligned BS data in multiple samples and outputs identified DM CG sites and regions.

HMM-DM requires a Linux/Unix system, with R installed. There are one document and two folders:
_____________________________________________________________________________________________________________
HMM.DM.user.manual.pdf:	A copy of the user manual
HMM.DM.code: A folder containing all R source code files used for HMM-DM.
example.data: A folder containing all example input data as mentioned in this document, an example.script.txt for running               HMM-DM, and the output files generated from the example.script.txt.
_____________________________________________________________________________________________________________
