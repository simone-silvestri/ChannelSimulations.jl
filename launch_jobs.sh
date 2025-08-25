#!/bin/bash

# TSP is the timestepper (0 => AB2, 1 => RK3)
# ZST is the z coord     (0 => Zstar, 1 => Zcoordinate)
# MOM is the mom advect  (0 => WENOVectorInvariant, 
#                         1 => VectorInvariant, 
#                         2 => WENOVectorInvariant Centered advecting velocity)
# CLO is the z closure   (0 => CA, 1 => CATKE, 2 => RiBased)
# TRA is the tracer adv  (0 => WENO7, 
#                         1 => (x, y => WENO5, z => Centered), 
# 			  2 => (x, y => WENO9, z => Centered), 
# 			  3 => Upwind3, 
# 			  4 => (x, y => Upwind3, z => Centered), 
# 			  5 => WENO5)

# EXP = string(CLO) * string(MOM) * string(TRA) * string(TSP) * string(ZST)

sbatch --export=CASE=00001 --ouput="out00001.txt" --error="err00001.txt" chan.sh 
sbatch --export=CASE=00011 --ouput="out00011.txt" --error="err00011.txt" chan.sh 
sbatch --export=CASE=01001 --ouput="out01001.txt" --error="err01001.txt" chan.sh 
sbatch --export=CASE=01011 --ouput="out01011.txt" --error="err01011.txt" chan.sh 

sbatch --export=CASE=10001 --ouput="out10001.txt" --error="err10001.txt" chan.sh 
sbatch --export=CASE=10011 --ouput="out10011.txt" --error="err10011.txt" chan.sh 
sbatch --export=CASE=11001 --ouput="out11001.txt" --error="err11001.txt" chan.sh 
sbatch --export=CASE=11011 --ouput="out11011.txt" --error="err11011.txt" chan.sh 

sbatch --export=CASE=00301 --ouput="out00301.txt" --error="err00301.txt" chan.sh 
sbatch --export=CASE=00311 --ouput="out00311.txt" --error="err00311.txt" chan.sh 
sbatch --export=CASE=01301 --ouput="out01301.txt" --error="err01301.txt" chan.sh 
sbatch --export=CASE=01311 --ouput="out01311.txt" --error="err01311.txt" chan.sh 
