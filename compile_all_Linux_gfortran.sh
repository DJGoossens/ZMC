#! /bin/bash
#
# Manual compile of ZMC, DZMC and toolbox programs.
#
# This is for Linux/Mac (POSIX-type) environments.
# I suggest using 'compile_all.bat' for windows.
# (If you don't have the file, ask me.)
#
# Dependencies are accounted for by the order of the
# compilations -- some modules depend on others.
# This means that this could be broken by upstream changes
# to the modules by Aidan Heerdegen (aidan@rsc.anu.edu.au,
# though this email may not be current.)
#
# However, if this does not happen, compilation is simple.
#
# NOTE!!: Some manual customisation of this file will
# probably be needed, but it is very simple:
#
# Namely the -I and -L flags to gfortran need to point
# to the directory where the compilation is happening,
# and that it is easiest for a relatively simple project
# like this to use the flat structure, then just copy out
# the executables to a directory in the path.
#
# No guarantees/warranties are given or implied.  Use at
# your own risk.  Caveat emptor.  Buyer beware.
#
# And you get what you pay for...
#
# d.goossens@adfa.edu.au Nov 2015  
# darren.goossens@gmail.com  2026  

# Set compilation flags for everything
#GFFLAGS="-static-libgfortran -O2 -fallow-argument-mismatch -ffree-line-length-0"
GFFLAGS="-static-libgfortran -O2  -ffree-line-length-0"

# Use these for testing
#GFFLAGS="-static-libgfortran -Wall -Wextra -fcheck=all -ffree-line-length-0"


#
echo ---------------------------
COMPILEDIR=`pwd`
echo Compiling in $COMPILEDIR
# Clean up previous; both cases in case of files leftover from a non-POSIX compile
rm *.o *.O *.mod *.MOD
bindir=ZMCLinux_gfortran
echo ---------------------------
mkdir -v $bindir
echo ---------------------------
#
# Echo everything
set -x verbose
gfortran $GFFLAGS -c precision.f90
gfortran $GFFLAGS -c fundamental_constants.f90
gfortran $GFFLAGS -c varying_string.f90
gfortran $GFFLAGS -c globals.f90
gfortran $GFFLAGS -c cartesian_class.f90
gfortran $GFFLAGS -c array_functions.f90

gfortran $GFFLAGS -c sort_functions.f90
gfortran $GFFLAGS -c variable_array.f90
gfortran $GFFLAGS -c binary.f90

gfortran $GFFLAGS -c string_functions.f90
gfortran $GFFLAGS -c statistics.f90

gfortran $GFFLAGS -c binary_io.f90

gfortran $GFFLAGS -c file_functions.f90

gfortran $GFFLAGS -c hash_table.f90
gfortran $GFFLAGS -c polysample.f90
gfortran $GFFLAGS -c vector_class.f90

gfortran $GFFLAGS -c mol2_class.f90
gfortran $GFFLAGS -c cmdline_arguments.f90
gfortran $GFFLAGS -c keyword_class.f90
gfortran $GFFLAGS -c rotmatrix_class.f90

gfortran $GFFLAGS -c crystallography_class.f90
gfortran $GFFLAGS -c image_transforms.f90
gfortran $GFFLAGS -c quaternion_class.f90
gfortran $GFFLAGS -c zmatrix_class.f90

gfortran $GFFLAGS -c pnm_class.f90
gfortran $GFFLAGS -c superimpose.f90



gfortran $GFFLAGS   -o ZMC zmc_06_March_2026.f90 rannum_2026.f ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o
gfortran $GFFLAGS       -o DZMC  readat_zmc_Aug29_2014.f90 diffuse_allocatable_March20_2014.f90 rannum_2026.f  ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o
gfortran $GFFLAGS       -o bin2gray   bin2gray.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o 
gfortran $GFFLAGS       -o zmat_maker  zmat_maker.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o 
#gfortran $GFFLAGS       -o zmat2xyz  zmat2xyz.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o 
echo =======================================================
echo NOT COMPILING ZMAT2XYZ -- IT IS BUGGY AS OF MARCH 2026
echo =======================================================
gfortran $GFFLAGS       -o zmatchk  zmatchk.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o 
gfortran $GFFLAGS       -o zmat_anim  zmat_anim.f90 rannum_2026.f  ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o 
#gfortran $GFFLAGS       -o zmat2mol2  zmat2mol2.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o 
echo =======================================================
echo NOT COMPILING ZMAT2MOL2 -- IT IS BUGGY AS OF MARCH 2026
echo =======================================================
gfortran $GFFLAGS       -o pgmave  pgmave.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o 
gfortran $GFFLAGS       -o pgmcombine  pgmcombine.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o 
gfortran $GFFLAGS       -o make_random_occ  make_random_occ.f90 rannum_2026.f  ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o 
gfortran $GFFLAGS       -o mol2xyz  mol2xyz.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o
gfortran $GFFLAGS       -o catmol2  catmol2.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o
gfortran $GFFLAGS       -o chkmol2  chkmol2.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o
gfortran $GFFLAGS       -o pgm2mask  pgm2mask.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o
gfortran $GFFLAGS       -o pgm2ni  pgm2ni.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o
gfortran $GFFLAGS       -o ni2pgm  ni2pgm.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o
gfortran $GFFLAGS       -o raw2pgm  raw2pgm.f90   ps_routines.f fundamental_constants.o varying_string.o globals.o precision.o cartesian_class.o array_functions.o sort_functions.o variable_array.o binary.o binary_io.o string_functions.o statistics.o file_functions.o hash_table.o polysample.o vector_class.o mol2_class.o cmdline_arguments.o keyword_class.o rotmatrix_class.o crystallography_class.o image_transforms.o quaternion_class.o zmatrix_class.o superimpose.o pnm_class.o
mv raw2pgm $bindir
mv ni2pgm $bindir
mv  pgm2ni $bindir
mv  pgm2mask  $bindir
mv  chkmol2   $bindir
mv  catmol2   $bindir
mv  mol2xyz   $bindir
mv  make_random_occ $bindir
mv  pgmcombine $bindir
mv  pgmave    $bindir
mv  zmat2mol2 $bindir
mv  zmat_anim $bindir
mv  zmatchk   $bindir
mv   zmat2xyz $bindir
mv   zmat_maker $bindir
mv   bin2gray   $bindir
mv   DZMC       $bindir
mv   ZMC       $bindir
cp DISCLAIMER.txt $bindir
cp COPYRIGHT.txt $bindir
cp READ* $bindir
