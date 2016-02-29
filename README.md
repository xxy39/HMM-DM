HMM-DM
======
The HMM-DM program identifies differentially methylated (DM) CG sites and regions from whole genome and targeted bisulfite sequencing (BS) data. This approach first uses a Hidden Markov Model to identify differentially methylated CG sites accounting for spatial correlation across CGs and variation across samples, and then summarizes identified DM CG sites into regions based on their status and distance. This program takes aligned BS data in multiple samples and outputs identified DM CG sites and regions.

HMM-DM requires R installed. Ideally it is run in a Linux/Unix system. This program includes the following documents and folders:
_____________________________________________________________________________________________________________
 
HMM.DM.user.manual.pdf:	A copy of the user manual
 
HMM.DM.code: A folder containing all R source code files used for HMM-DM.

example.data: A folder containing all example input data, an example.script.txt for running HMM-DM, and the output files generated from the example.script.txt.

Simulation: A folder containing all R scource code used for data simulating.
_____________________________________________________________________________________________________________

**How to cite us**

Yu, X. & Sun, S. (2016). HMM-DM: identifying differentially methylated regions using a hidden Markov model. Statistical Applications in Genetics and Molecular Biology, Doi:10.1515/sagmb-2015-0077

**Access to the published manuscripts**

HMM-DM: http://www.ncbi.nlm.nih.gov/pubmed/26887041

HMM-Fisher: http://www.ncbi.nlm.nih.gov/pubmed/26854292

Comparing five DM methods: http://www.ncbi.nlm.nih.gov/pubmed/26910753
