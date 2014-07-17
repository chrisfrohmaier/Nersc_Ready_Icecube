#!/bin/bash
#PBS -S /bin/bash
#PBS -N main_XXX 
#PBS -j oe 
#PBS -l nodes=1:ppn=8,walltime=04:00:00
#PBS -q regular
#PBS -A m937


cd /project/projectdirs/deepsky/rates/icecube/fakes

./subit.first ptf_XXX
./subit.last ptf_XXX

