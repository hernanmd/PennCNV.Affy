#!/bin/bash
# Author: Hernán Morales Durand <hernan.morales@gmail.com>

# Uncomment to debug this script
# set -x

apt_geno=$(type -p apt-genotype-axiom)
apt_summarize=$(type -p apt-probeset-summarize)

###########################################
# Sanity checks
###########################################

echo -n "Sanity checks..."
[ -d ${ax_library_dir} ] ||  { echo "ERROR: Library directory not found"; exit 1; }
[ -d ${ax_cel_dir} ] || { echo "ERROR: CEL files directory not found"; exit 1; }
[ ! -z "${apt_geno}" ] || { echo "ERROR: APT Power Tools not found in PATH"; exit 1; }
[ ! -z "${apt_summarize}" ] || { echo "ERROR: APT Power Tools not found in PATH"; exit 1; }
echo "done"

###########################################
# Begin processing
###########################################

#
# 1. Make a file containing a list of all CEL files + location
#
echo -n "Building CEL list file..."
echo "cel_files" > $ax_cel_file
ls -d $ax_cel_dir/*.CEL >> $ax_cel_file
echo "done"

#
# 2. Generate genotyping calls from CEL files
# 
echo -n "Generate genotyping calls from CEL files..."
"$apt_geno" --cel-files ${ax_cel_file} --analysis-files-path ${ax_library_dir} --arg-file ${ax_params_file} --summaries --dual-channel-normalization true --write-models --batch-folder ${ax_suite_dir} --log-file ${ax_apt_geno_log}  --out-dir ${ax_results_dir}
echo "done"

#
# 3. Allele-specific signal extraction from CEL files
#
# The below will extract signal intensity values for PM probes in all the CEL files specified in the listfile. The values are then quantile normalized,  median polished, generating A & B signal intensity values for each SNP. The filehapmap.quant-norm.normalization-target.txt is is the reference quantile distribution.

"$apt_summarize" --cel-files ${ax_cel_file} --analysis "quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true" --target-sketch ${ax_sketch_file} --out-dir $ax_results_dir -v 2 --cdf-file $ax_probeset_file

#
# 4. create sex file from report
#
# The file_sex file is a two-column file that annotates the sex information for each CEL file. The birdseed.report.txt file that was generated from from the previous step can bes used to construct the sexfile based on computed gender information.

# %affymetrix-algorithm-param-apt-opt-igender-female-threshold=1.100000
# %affymetrix-algorithm-param-apt-opt-igender-male-threshold=1.200000
fgrep male ${ax_report_file} | cut -f 1,2 > ${computed_gender_file}

######################################################################################
######################################################################################
# Axiom CNV Tool processing (currently not available on command-line tool?)
######################################################################################
######################################################################################


