#!/usr/bin/env bash


BASEDIR=$(readlink -f "$0")

#Current directory
CURRDIR=$(dirname "$BASEDIR")

GSL='/gsl'
# GSL directory 
GSLDIR=$CURRDIR$GSL

echo $CURRDIR

#GSL library directory 
GSLLIB="$GSLDIR/lib"
#GSL header files directory 
GSLINCLUDE="$GSLDIR/include" 

export LD_LIBRARY_PATH=$GSLLIB 


# Link GSL header files 
cmd="g++ -Wall -I$GSLINCLUDE -c main.cpp"
echo $cmd
$cmd

# Link GSL libraries 
cmd="g++ -L$GSLLIB main.o -lgsl -lgslcblas -lm"
echo $cmd 
$cmd 

# Showing output
cmd="./a.out"
echo $cmd 
$cmd 




