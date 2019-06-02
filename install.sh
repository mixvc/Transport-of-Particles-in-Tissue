#!/usr/bin/env bash

# Unzip gsl
cmd="tar -zxvf gsl-2.1.tar.gz" 
$cmd 

cmd="rm -rf gsl" 
$cmd 

mkdir gsl



BASEDIR=$(readlink -f "$0")

#Current directory
CURRDIR=$(dirname "$BASEDIR")

GSL='/gsl'
# GSL directory 
GSLDIR=$CURRDIR$GSL

echo $CURRDIR
cmd="./gsl-2.1/configure --prefix=$GSLDIR"
$cmd 
make
sudo make install

