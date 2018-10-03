#!/bin/sh
# Author: Hernán Morales Durand <hernan.morales@gmail.com>	
#
# Software & Library requirements
#
#	ActivePerl: https://www.activestate.com/activeperl/downloads (required by PennCNV for Windows)
# 	PennCNV software & accessible through PATH
#		Windows: cmd.exe -> rundll32 sysdm.cpl,EditEnvironmentVariables
# 	Affymetrix Power Tools (APT) software & accessible through PATH
# 	Axiom BOS1 Library files
# 		Other library files (CDFs) at https://www.thermofisher.com/ar/es/home/life-science/microarray-analysis/microarray-data-analysis/genechip-array-library-files.html


###########################################
# Input parameters
###########################################

set -x

# Default directory for AAS: 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_GW_Bos_SNP_1.r3'
cel_files_dir="004"
pfb_file=".pfb"
results_dir="output/"

library_dir="/c/Users/Public/Documents/AxiomAnalysisSuite/Library/Axiom_GW_Bos_SNP_1.r3/"
xml_file=${library_dir}'Axiom_GW_Bos_SNP_1_96orMore_Step1.r3.apt-probeset-genotype.AxiomGT1.xml'
sketch_file=${library_dir}'Axiom_GW_Bos_SNP_1.r3.AxiomGT1.sketch'
probeset_file=${library_dir}"Axiom_GW_Bos_SNP_1.r3.cdf"
annot_db==${library_dir}'Axiom_GW_Bos_SNP_1.na35.annot.db'
axiom_report=${results_dir}'AxiomGT1.report.txt'
computed_gender_file=${results_dir}gender_computed.txt
list_file="cel_files_list.txt"

# PennCNV params for ActivePerl
perlexec="/c/Perl64/bin/perl"
gen_cluster_exec=$(which generate_affy_geno_cluster.pl)
norm_cluster_exec=$(which normalize_affy_geno_cluster.pl)

###########################################
# Begin
###########################################

[ -d $library_dir ] ||  { echo "ERROR: Library directory not found"; exit 1; }

#
# 1. Make a file containing a list of all CEL files + location
#

echo "cel_files" > $list_file
ls -d $cel_files_dir/*.CEL >> $list_file

#
# 2. Generate genotyping calls from CEL files
#

# apt-probeset-genotype -v 2 --log-file "apt-probeset-genotype.log" --analysis-files-path ${library_dir} --xml-file ${xml_file} --out-dir ${results_dir} --cel-files ${list_file}
apt-genotype-axiom --log-file "apt-genotype-axiom.log" --analysis-files-path ${library_dir} --arg-file ${xml_file} --out-dir ${results_dir} --cel-files ${list_file}

#
# 3. Allele-specific signal extraction from CEL files
#
# The below will extract signal intensity values for PM probes in all the CEL files specified in the listfile. The values are then quantile normalized,  median polished, generating A & B signal intensity values for each SNP. The filehapmap.quant-norm.normalization-target.txt is is the reference quantile distribution.

apt-probeset-summarize --analysis "quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true" --target-sketch $sketch_file --out-dir $results_dir --cel-files $list_file -v 2 --cdf-file $probeset_file

#
# 4. create sex file from report
#
# The file_sex file is a two-column file that annotates the sex information for each CEL file. The birdseed.report.txt file that was generated from from the previous step can bes used to construct the sexfile based on computed gender information.

fgrep male ${axiom_report} | cut -f 1,2 > ${computed_gender_file}

#
# 5. Generate canonical genotype clustering file
#
# The following step will generate canonical genotype clusters.
# file_sex = two-column file that annotates the sex information for each CEL file
# 10918.CEL       male
# 10924.CEL       male
# 11321_2.CEL     female

qnpmpes='quant-norm.pm-only.med-polish.expr.summary.txt'
axcalls='AxiomGT1.calls.txt'
axconfs='AxiomGT1.confidences.txt'
genoclust=${results_dir}.genocluster

$perlexec $gen_cluster_exec ${axcalls} ${axconfs} ${qnpmpes} -locfile ${pfb_file} -sexfile ${computed_gender_file} -out ${genoclust}

#
# 6. Calculate LRR and BAF
#
# Finally, we can extract allele-specific signal intensity measures to calculate the Log R Ratio (LRR) values and the B Allele Frequency (BAF) values for each marker in each individual. The cluster file gw6.genocluster is a file which specifies the chromosome position of each SNP or CN probe.

$perlexec $norm_cluster_exec ${genoclust} ${qnpmpes} -locfile ${pfb_file} -out ${results_dir}.lrr_baf.txt
