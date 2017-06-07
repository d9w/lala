#!/bin/bash
# source me
module purge
module load gcc/5.3.0
module load intel/14.0.2.144
export CC=gcc
export CXX=g++
export PATH=/usr/local/julia/0.5.0_stand/julia-3c9d75391c/bin:$PATH
