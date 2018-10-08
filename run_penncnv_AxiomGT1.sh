#!/bin/sh
# Author: Hernán Morales Durand <hernan.morales@gmail.com>	
#
# Software & Library requirements
#
#	Put your .CEL files in a subdirectory of this script directory and edit the cel_files_dir variable below.
#	If you are on Windows download and install ActivePerl: https://www.activestate.com/activeperl/downloads (required by PennCNV for Windows)
# 	Install PennCNV software
#		git clone https://github.com/WGLab/PennCNV.git
# 	Install Affymetrix Power Tools (APT) software 
#		https://www.thermofisher.com/ar/es/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html & 
# 	Download Axiom Library files (CDFs) 
#		https://www.thermofisher.com/ar/es/home/life-science/microarray-analysis/microarray-data-analysis/genechip-array-library-files.html

###########################################
# Input parameters
###########################################

# Uncomment to debug this script
#set -x

prj_prefix=$1
penncnv_dir='/c/git_projects/PennCNV/affy/bin/'
PATH=$PATH:$penncnv_dir:/c/Program\ Files/Thermo\ Fisher\ Scientific/Analysis\ Power\ Tools/APT-2.10.0/bin/

# Default directory for AAS: 'c:\Users\Public\Documents\AxiomAnalysisSuite\Library\Axiom_GW_Bos_SNP_1.r3'
cel_files_dir='004/'
cel_list_file=$(pwd)'cel_files_list.txt'
results_dir="output/"
pfb_file="AxiomGT1.pfb"

ax_library_dir='/c/Users/Public/Documents/AxiomAnalysisSuite/Library/Axiom_GW_Bos_SNP_1.r3/'
ax_params_file=${ax_library_dir}'Axiom_GW_Bos_SNP_1_96orMore_Step1.r3.apt-probeset-genotype.AxiomGT1.apt2.xml'
ax_sketch_file=${ax_library_dir}'Axiom_GW_Bos_SNP_1.r3.AxiomGT1.sketch'
ax_probeset_file=${ax_library_dir}"Axiom_GW_Bos_SNP_1.r3.cdf"
ax_annot_db==${ax_library_dir}'Axiom_GW_Bos_SNP_1.na35.annot.db'
ax_report_file=${results_dir}'AxiomGT1.report.txt'
ax_suite_dir=${results_dir}'suitefiles'
ax_qnpmpes_file=${results_dir}'quant-norm.pm-only.med-polish.expr.summary.txt'
ax_calls_file=${results_dir}'AxiomGT1.calls.txt'
ax_confs_file=${results_dir}'AxiomGT1.confidences.txt'
genoclust=${results_dir}'AxiomGT1.genocluster'

computed_gender_file=${results_dir}'gender_computed.txt'
apt_geno=$(type -p apt-genotype-axiom)
apt_summarize=$(type -p apt-probeset-summarize)

##############################################
# PennCNV general parameters
##############################################

if [ -f "/c/Perl64/bin/perl" ]; then
	perlexec="/c/Perl64/bin/perl"
else
	perlexec=$(which perl)
fi
gen_cluster_exec=$(type -p generate_affy_geno_cluster.pl)
norm_cluster_exec=$(type -p normalize_affy_geno_cluster.pl)


##############################################
# PennCNV: Detect CNV parameters
##############################################

# Call CNVs containing more than or equal to 3 SNPs
minsnp=3
# Last chromosome nr
lastchr=31

############ For main HMM file
hmm1="/usr/local/src/PennCNV-1.0.4/lib/hhall.hmm"
log_file1=hmm1.minsnp_"$minsnp".log
raw_file1=hmm1.minsnp_"$minsnp".rawcnv
qc_log_file1=hmm1_minsnp_"$minsnp".log
qclrrsd1=0.3
qc_passout1=hmm1_minsnp_"$minsnp".qcpass
qc_sumout1=hmm1_minsnp_"$minsnp".qcsum
qc_goodcnv1=hmm1_minsnp_"$minsnp".goodcnv

##############################################################
# Output files:
##############################################################

pfb_file=${prj_prefix}".pfb"
gc_file_model=${prj_prefix}".gcmodel"
signal_outdir="PennCNV_run1/"
signal_outsuffix=".txt"
signal_file_list="signal_fltr_files.txt"

###########################################
# Sanity checks
###########################################

[ -d ${ax_library_dir} ] ||  { echo "ERROR: Library directory not found"; exit 1; }
[ -d ${cel_files_dir} ] || { echo "ERROR: CEL files directory not found"; exit 1; }
[ ! -z ${gen_cluster_exec} ] || { echo "ERROR: PennCNV scripts not found in PATH"; exit 1; }
[ ! -z ${norm_cluster_exec} ] || { echo "ERROR: PennCNV scripts not found in PATH"; exit 1; }
[ ! -z "${apt_geno}" ] || { echo "ERROR: APT Power Tools not found in PATH"; exit 1; }
[ ! -z "${apt_summarize}" ] || { echo "ERROR: APT Power Tools not found in PATH"; exit 1; }

###########################################
# Begin processing
###########################################

#
# 1. Make a file containing a list of all CEL files + location
#
echo "cel_files" > $cel_list_file
ls -d $cel_files_dir/*.CEL >> $cel_list_file

#
# 2. Generate genotyping calls from CEL files
#  --summaries True  --console-add-select debug

"$apt_geno" --cel-files ${cel_list_file} --analysis-files-path ${ax_library_dir} --arg-file ${ax_params_file} --summaries --dual-channel-normalization true --write-models --batch-folder ${ax_suite_dir} --log-file "apt-genotype-axiom.log"  --out-dir ${results_dir}

#
# 3. Allele-specific signal extraction from CEL files
#
# The below will extract signal intensity values for PM probes in all the CEL files specified in the listfile. The values are then quantile normalized,  median polished, generating A & B signal intensity values for each SNP. The filehapmap.quant-norm.normalization-target.txt is is the reference quantile distribution.

"$apt_summarize" --cel-files ${cel_list_file} --analysis "quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true" --target-sketch ${ax_sketch_file} --out-dir $results_dir -v 2 --cdf-file $ax_probeset_file

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



##############################################
# Compile PFB
##############################################

# No .pfb files? Use compile_pfb.pl
# ls='ls -lkF --color=auto'
echo "Creating signal file list"
#unalias ls
signal_files=$(ls -1 $signal_output_dir)
echo "Created signal file list: $signal_files"

# The output is a new file with directory/signal_file_name in each line
for f in $signal_files; do
	echo $signal_output_dir/$f >> $signal_file_list
done

# SNP_Map.txt from SNP_Map.zip in the Illumina raw files
# The output is a .pfb file
echo "About compiling PFB..."
compile_pfb.pl \
	--listfile $signal_file_list \
	--snpposfile $snp_map \
	--output $pfb_file
echo "done compile PFB"


#
# 5. Generate canonical genotype clustering file
#
# The following step will generate canonical genotype clusters.
# file_sex = two-column file that annotates the sex information for each CEL file
# 10918.CEL       male
# 11321_2.CEL     female

# $perlexec $gen_cluster_exec ${ax_calls_file} ${ax_confs_file} ${ax_qnpmpes_file} -locfile ${pfb_file} -sexfile ${computed_gender_file} -out ${genoclust}

#
# 6. Calculate LRR and BAF
#
# Finally, we can extract allele-specific signal intensity measures to calculate the Log R Ratio (LRR) values and the B Allele Frequency (BAF) values for each marker in each individual. The cluster file gw6.genocluster is a file which specifies the chromosome position of each SNP or CN probe.

# $perlexec $norm_cluster_exec ${genoclust} ${ax_qnpmpes_file} -locfile ${pfb_file} -out ${results_dir}.lrr_baf.txt
