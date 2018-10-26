#!/bin/sh
# Author: Hernán Morales Durand <hernan.morales@gmail.com>	
#

# Uncomment to debug this script
#set -x

###########################################
# Input parameters
###########################################

case "$OSTYPE" in
	linux*|darwin*)
		penncnv_root="/usr/local/src/PennCNV/"
		penncnv_axdir=${penncnv_root}"affy/bin/"
		# APT commands should be already available in PATH as indicated in INSTALL file
		PATH=$PATH:$penncnv_root:$penncnv_axdir
		;;
	msys*)
		penncnv_root="/c/PennCNV/"	
		penncnv_axdir=${penncnv_root}"affy/bin/"
		PATH=$PATH:$penncnv_root:$penncnv_axdir:/c/Program\ Files/Thermo\ Fisher\ Scientific/Analysis\ Power\ Tools/APT-2.10.0/bin/		
		;;
	*)
		echo "unknown: $OSTYPE"
		exit 1
		;;			
esac

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

###########################################
# Sanity checks
###########################################

echo -n "Sanity checks..."
[ ! -z ${gen_cluster_exec} ] || { echo "ERROR: PennCNV scripts not found in PATH"; exit 1; }
[ ! -z ${norm_cluster_exec} ] || { echo "ERROR: PennCNV scripts not found in PATH"; exit 1; }
echo "done"

ucsc_gen=1
ax_penncnv_map=${ax_results_dir}"Axiom_GW_Bos_SNP_1.na35.annot.db-probemappings.txt"
ax_qnpmpes_file=${ax_results_dir}'quant-norm.pm-only.med-polish.expr.summary.txt'
ax_calls_file=${ax_results_dir}'AxiomGT1.calls.txt'
ax_confs_file=${ax_results_dir}'AxiomGT1.confidences.txt'
gc_file_prefix="gc5Base"
gc_file=$gc_file_prefix".txt"
gc_file_sorted=$gc_file_prefix".txt"
gc_file_gz=$gc_file".gz"

# If your species has already generated a sorted GC file from UCSC: edit the path to match your species
# If the GC file contains only a pointer to the bigWig file (like /gbdb/mm10/bbi/gc5Base.bw) read below
gc_file_url="http://hgdownload.cse.ucsc.edu/goldenPath/equCab2/database/"$gc_file_gz
gc_file_url="http://hgdownload.cse.ucsc.edu/goldenPath/loxAfr3/database/"$gc_file_gz

# If your species has not a sorted GC file from UCSC, it should be available as bigWig file: http://hgdownload.cse.ucsc.edu/gbdb/
gc_bigWig=${gc_file_prefix}".bw"
gc_wig=${gc_file_prefix}".wig"


# For annotating with scan_region.pl
refGeneUrl="http://hgdownload.soe.ucsc.edu/goldenPath/equCab2/database/refGene.txt.gz"


##############################################################
# Output files:
##############################################################

gc_file_model=${prj_prefix}".gcmodel"
signal_outdir=${prj_prefix}"_signals/"
signal_outsuffix=".txt"
signal_file_list="signal_fltr_files.txt"

##############################################
# PennCNV: Detect CNV parameters
##############################################

pfb_file=${prj_prefix}"_AxiomGT1.pfb"
pfb_file_cleaned=${prj_prefix}"_AxiomGT1.cleaned.pfb"
penncnv_hmm1=${penncnv_root}"lib/hhall.hmm"
penncnv_results_dir="penncnv_output/"

# Call CNVs containing more than or equal to 3 SNPs
minsnp=3
# Last chromosome nr
lastchr=31

############ For main HMM file

log_file1=penncnv_hmm1.minsnp_"$minsnp".log
raw_file1=penncnv_hmm1.minsnp_"$minsnp".rawcnv
qc_log_file1=penncnv_hmm1_minsnp_"$minsnp".log
qclrrsd1=0.3
qc_passout1=penncnv_hmm1_minsnp_"$minsnp".qcpass
qc_sumout1=penncnv_hmm1_minsnp_"$minsnp".qcsum
qc_goodcnv1=penncnv_hmm1_minsnp_"$minsnp".goodcnv


##############################################
# Compile PFB
##############################################

# Clean previous run and create output directories
#rm -frv $signal_outdir
rm -fv ${signal_file_list}
#mkdir $signal_outdir

# No .pfb files? Use compile_pfb.pl
# ls='ls -lkF --color=auto'
echo "Creating signal file list"
# ls is aliased? 
#unalias ls
signal_files=$(ls -1 $signal_outdir)

# The output is a new file with directory/signal_file_name in each line
for f in ${signal_files}; do
	echo ${signal_outdir}/${f} >> ${signal_file_list}
done
echo "Created signal file list: $signal_files"

# The output is a .pfb file
echo "About compiling PFB..."
mkdir -p ${penncnv_results_dir}
${perlexec} ${penncnv_root}compile_pfb.pl --listfile ${signal_file_list} --snpposfile ${ax_penncnv_map} --output ${penncnv_results_dir}/${pfb_file}
echo "done compile PFB in : "${penncnv_results_dir}/${pfb_file}

# Sort and clean PFB
echo "Sort and clean PFB file..."
sort -k 2,2n -k 3,3n < ${penncnv_results_dir}/${pfb_file} | awk '{ if (NR == 1 || $2 >= 1 && $2 <= '${lastchr}') print }' > > ${penncnv_results_dir}/${pfb_file_cleaned}
echo "done cleaning PFB in :"${penncnv_results_dir}/${pfb_file_cleaned}

##############################################
# Obtain or generate GC model file
##############################################

if [ ${ucsc_gen} = true ]; then
	# Download and sort GC file if not found
	[ -f $gc_file_sorted ] || { wget $gc_file_url; gunzip $gc_file_gz; sort -k 2,2 -k 3,3n < $gc_file > $gc_file_sorted; }
else
	# GC 5 Base = GC content per 5120bp
	bigWigToWig ${gc_bigWig} stdout | gzip -c > ${gc_file_prefix}.varStep.txt.gz
	wigEncode ${gc_file_prefix}.varStep.txt.gz ${gc_wig} ${gc_file_prefix}.wib
	# Add dummy column to the wiggle file.
	# See https://groups.google.com/a/soe.ucsc.edu/forum/?utm_medium=email&utm_source=footer#!topic/genome/2Jp8cGCsKBU for details
	mv ${gc_wig} ${gc_file_prefix}.tmp
	sed 's/^/2112\t/' ${gc_file_prefix}.tmp > ${gc_wig}
	gc_file_sorted=${gc_wig}
fi

# Make GC model
echo "About creating GC Model..."
cal_gc_snp.pl ${gc_file_sorted} ${pfb_file} -output ${gc_file_model}
echo "done"

#
# 5. Generate canonical genotype clustering file
#
# The following step will generate canonical genotype clusters.
# file_sex = two-column file that annotates the sex information for each CEL file
# 10918.CEL       male
# 11321_2.CEL     female

$perlexec $gen_cluster_exec ${ax_calls_file} ${ax_confs_file} ${ax_qnpmpes_file} -locfile ${pfb_file} -sexfile ${computed_gender_file} -out ${genoclust}

#
# 6. Calculate LRR and BAF
#
# Finally, we can extract allele-specific signal intensity measures to calculate the Log R Ratio (LRR) values and the B Allele Frequency (BAF) values for each marker in each individual. The cluster file gw6.genocluster is a file which specifies the chromosome position of each SNP or CN probe.

# $perlexec $norm_cluster_exec ${genoclust} ${ax_qnpmpes_file} -locfile ${pfb_file} -out ${ax_results_dir}.lrr_baf.txt
