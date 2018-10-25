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
  - [Affymetrix Power Tools (APT) software](https://www.thermofisher.com/ar/es/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html)
  - UCSC utilities:
    - Linux: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
    - MacOS: http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/
	- Alternative download: `git clone https://github.com/ENCODE-DCC/kentUtils`
	- Copy or move the downloaded commands to a PATH location:
	  - `echo $PATH | tr \: \\n`

# Library Requirements

  - Download [Axiom Library files](https://www.thermofisher.com/ar/es/home/life-science/microarray-analysis/microarray-data-analysis/genechip-array-library-files.html) for your species (CDFs) 
    - Linux/MacOS : AAS Library files should be installed in /usr/local/src/AxiomAnalysisSuite/Library/
    - Windows: AAS Library files should be installed in c:\Users\Public\Documents\AxiomAnalysisSuite\Library\ or wherever the ax_library_root variable point to.
      - Library files example names:
      - 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_GW_Bos_SNP_1.r3'
      - 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_Buffalo_Analysis-r2'
      - 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_GW_GT_Chicken_Analysis-r2'
      - 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_OviCap_Analysis-caprine-r1'
      - 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_MouseHD_Analysis.r1'  
	
# Script configuration
	
```bash
git clone https://github.com/hernanmd/PennCNV.Affy.git
cd PennCNV.Affy
```

  - Put your .CEL files in a subdirectory of this script directory 
    - Edit the ax_cel_dir variable with the subdirectory name in the script run_AxiomGT1_part1.sh.
  - Edit or check the following variables in Part 1 script: 
    - prj_prefix = Prefix of output files
    - penncnv_root = Directory with PennCNV binaries
    - ax_library_root = Parent directory of Axiom Libraries
    - ax_library_dir = Directory of Axiom Library files for your species (see above Library Requirements)
    - ax_params_file = XML file parameters, ex: Axiom_GW_Bos_SNP_1_96orMore_Step1.r3.apt-probeset-genotype.AxiomGT1.apt2.xml
    - ax_sketch_file = Axiom sketch file: Axiom_GW_Bos_SNP_1.r3.AxiomGT1.sketch
    - ax_probeset_file = Axiom CDF file
    - ax_annot_db = SQLite annotation database file (.annot.db)
  - Edit in Part 2 script:
    - ax_penncnv_map = Axiom probe mappings file (.annot.db-probemappings.txt)
	- ucsc_gen = true if your GC file for cal_gc_snp.pl was generated by UCSC, false if is not present in UCSC or is a pointer to a bigWig file (.bw)
	  - If false, you should run the Part 2 script under Linux or MacOS since UCSC are not available on Windows.

# Usage

```bash
./run_AxiomGT1_part1.sh
```

  - The first script generates a subdirectory ax_output/ (or whatever name in the ax_results_dir).
  - As next step you should run the Axiom CNV Tools which will generate a set of penncnv.cnv.txt files.
  - Move the set of penncnv.cnv.txt files to the ax_output directory.
  

```bash
./run_AxiomGT1_part2.sh
```
