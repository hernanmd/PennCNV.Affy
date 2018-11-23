#!/bin/bash

# Shell script to find the Axiom root library directory on current platform
# Meant to be sourced from other scripts.

case "$OSTYPE" in
	linux*|darwin*)
		ax_library_root="/usr/local/src/AxiomAnalysisSuite/Library/"
		# APT commands should be already available in PATH as indicated in INSTALL file
		PATH=$PATH:$penncnv_root:$penncnv_axdir
		;;
	msys*)
		ax_library_root="/c/Users/Public/Documents/AxiomAnalysisSuite/Library/"
		PATH=$PATH:$penncnv_root:$penncnv_axdir:/c/Program\ Files/Thermo\ Fisher\ Scientific/Analysis\ Power\ Tools/APT-2.10.0/bin/
		;;
	*)
		echo "unknown: $OSTYPE"
		exit 1
		;;
esac
