#!/bin/bash

# Shell script to find the PennCNV root directory on current platform
# Meant to be sourced from other scripts.

# Usually does not need to be changed
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
		echo "unknown OS: $OSTYPE"
		exit 1
		;;
esac