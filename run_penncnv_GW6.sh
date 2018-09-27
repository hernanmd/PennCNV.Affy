#!/bin/sh


library_dir='/c/Users/Public/Documents/AxiomAnalysisSuite/Library/'
probeset_file=$library_dir'Axiom_GW_Bos_SNP_1.r3.cdf'		
spsnps_file=$library_dir'Axiom_GW_Bos_SNP_1.r3.specialSNPs'
models_file=$library_dir'Axiom_GW_Bos_SNP_1.r3.AxiomGT1.models'
sketch_file=$library_dir'Axiom_GW_Bos_SNP_1.r3.AxiomGT1.sketch'

# Genome-wide 5.0 array
apt-probeset-genotype \
	-c lib/GenomeWideSNP_5.Full.r2.cdf \
	-a brlmm-p \
	--read-models-brlmmp $models_file \
	--out-dir apt \
	--cel-files $list_file \
	--chrX-snps lib/GenomeWideSNP_5.Full.chrx \

# Genome-wide 6.0 array
apt-probeset-genotype \
	-o $results_dir \
	-c $probeset_file \
	-a birdseed \
	--read-models-birdseed $models_file \
	--special-snps $spsnps_file \
	--cel-files $list_file