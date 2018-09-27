#!/bin/sh

# Software & Library requirements
#
# 	PennCNV software
# 	Affymetrix Power Tools (APT) software
# 	Axiom BOS1 Library files
# 		Other library files (CDFs) at https://www.thermofisher.com/ar/es/home/life-science/microarray-analysis/microarray-data-analysis/genechip-array-library-files.html

###########################################
# Input parameters
###########################################

# Default directory for AAS: 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_GW_Bos_SNP_1.r3'
library_dir='/c/Users/Public/Documents/AxiomAnalysisSuite/Library/'
xml_file=$library_dir'Axiom_GW_Bos_SNP_1_96orMore_Step1.r3.apt-probeset-genotype.AxiomGT1.xml'
list_file='cel_files_list.txt'
pfb_file='.pfb'
results_dir='output'

#
# 1. Make a file containing a list of all CEL files + location
#

ls -d $PWD/*.CEL > $list_file
sed -i '1s/^/cel_files/' $list_file

#
# 2. Generate genotyping calls from CEL files
#

apt-probeset-genotype \
	--analysis-files-path ${library_dir} \
	--xml-file ${xml_file} \
	--out-dir ${results_dir} \
	--cel-files ${$list_file}	
	
	#
# 3. Allele-specific signal extraction from CEL files
#
# The below will extract signal intensity values for PM probes in all the CEL files specified in the listfile. The values are then quantile normalized,  median polished, generating A & B signal intensity values for each SNP. The filehapmap.quant-norm.normalization-target.txt is is the reference quantile distribution from the HapMap3 project.

# Genome-wide 6.0 array
apt-probeset-summarize \
	--cdf-file $probeset_file \
	--analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true \
	--target-sketch $sketch_file \
	--out-dir $results_dir \
	--cel-files $list_file
	
#
# 4. create sex file from birdseed report
#
# The file_sex file is a two-column file that annotates the sex information for each CEL file. The birdseed.report.txt file that was generated from from the previous step can bes used to construct the sexfile based on computed gender information.

fgrep male ${results_dir}birdseed.report.txt | cut -f 1,2 > ${results_dir}gender_computed.txt


#
# 5. Generate canonical genotype clustering file
#
# The following step will generate canonical genotype clusters. 
# file_sex = two-column file that annotates the sex information for each CEL file
# 10918.CEL       male
# 10924.CEL       male
# 11321_2.CEL     female
qnpmpes='quant-norm.pm-only.med-polish.expr.summary.txt'
brcalls='birdseed.calls.txt'
brconfs='birdseed.confidences.txt'
genoclust=${results_dir}.genocluster

generate_affy_geno_cluster.pl \
	${brcalls} \
	${brconfs} \
	${qnpmpes} \
	-locfile ${pfb_file} \
	-sexfile ${results_dir}gender_computed.txt \
	-out ${genoclust}
	
#
# 6. Calculate LRR and BAF
#
# Finally, we can extract allele-specific signal intensity measures to calculate the Log R Ratio (LRR) values and the B Allele Frequency (BAF) values for each marker in each individual. The cluster file gw6.genocluster is a file which specifies the chromosome position of each SNP or CN probe.

normalize_affy_geno_cluster.pl \
	${genoclust} ${qnpmpes} \
	-locfile ${pfb_file}
	-out ${results_dir}.lrr_baf.txt

