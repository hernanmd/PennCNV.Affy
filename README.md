# Description

Workflow for running PennCNV with Affymetrix platform files

# Introduction

The workflow requires to open and edit the script files to customize the variables for your environment and species. 

The first script runs the Affymetrix Power Tools command-line programs to generate genotyping calls from CEL files, and extract signal intensity values for PM probes in all the CEL files, the values are then quantile normalized, median polished, generating A & B signal intensity values for each SNP.

The second script runs the PennCNV calculating LRR and BAF.

# Software Requirements

  - Git
    - Linux/MacOS: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
	- Windows: GitBash https://git-scm.com/download/win
  - Perl
    - Linux/MacOS: Already installed.
    - Windows: Download and install ActivePerl: https://www.activestate.com/activeperl/downloads (required by PennCNV)
  - PennCNV: `git clone https://github.com/WGLab/PennCNV.git`
  - (Affymetrix Power Tools (APT) software)[https://www.thermofisher.com/ar/es/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html]
  - UCSC utilities:
    - Linux: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
    - MacOS: http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/
	- Alternative download: `git clone https://github.com/ENCODE-DCC/kentUtils`
  - Clone this script: `git clone https://github.com/hernanmd/PennCNV.Affy.git`

# Library Requirements

  - Download Axiom Library files for your species (CDFs https://www.thermofisher.com/ar/es/home/life-science/microarray-analysis/microarray-data-analysis/genechip-array-library-files.html
    - Linux/MacOS : AAS Library files should be installed in /usr/local/src/AxiomAnalysisSuite/Library/
    - Windows: AAS Library files should be installed in c:\Users\Public\Documents\AxiomAnalysisSuite\Library\
      - Library files example names:
      - 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_GW_Bos_SNP_1.r3'
      - 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_Buffalo_Analysis-r2'
      - 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_GW_GT_Chicken_Analysis-r2'
      - 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_OviCap_Analysis-caprine-r1'
      - 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_MouseHD_Analysis.r1'  
	
# Script configuration
	
  - Put your .CEL files in a subdirectory of this script directory 
    - Edit the ax_cel_dir variable with the subdirectory name in the script run_AxiomGT1_part1.sh.
  - Edit or check the following variables: 
    - prj_prefix
    - penncnv_root
    - ax_library_root 
    - ax_library_dir 
    - ax_params_file 
    - ax_sketch_file 
    - ax_probeset_file 
    - ax_annot_db 
    - ax_penncnv_map 
