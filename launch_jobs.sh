#!/bin/bash

# TSP is the timestepper (0 => AB2, 1 => RK3)
# ZST is the z coord     (0 => Zstar, 1 => Zcoordinate)
# MOM is the mom advect  (0 => WENOVectorInvariant, 
#                         1 => VectorInvariant, 
#                         2 => WENOVectorInvariant Centered advecting velocity)
# CLO is the z closure   (0 => CA, 1 => CATKE, 2 => RiBased)
# TRA is the tracer adv  (0 => WENO7 - Buffeted, 
#                         1 => (x, y => WENO5, z => Centered), 
# 			  2 => (x, y => WENO9, z => Centered), 
# 			  3 => Upwind3, 
# 			  4 => (x, y => Upwind3, z => Centered), 
# 			  5 => WENO5)

# EXP = string(CLO) * string(MOM) * string(TRA) * string(TSP) * string(ZST)
sbatch -J j10000 --export=CASE=10000 --output="out10000.txt" --error="err10000.txt" chan.sh 
sbatch -J j10010 --export=CASE=10010 --output="out10010.txt" --error="err10010.txt" chan.sh 

# sbatch --export=CASE=00020 --output="out00020.txt" --error="err00020.txt" chan.sh
# 
# sbatch --export=CASE=00000 --output="out00000.txt" --error="err00000.txt" chan.sh 
# sbatch --export=CASE=00010 --output="out00010.txt" --error="err00010.txt" chan.sh 
# sbatch --export=CASE=01000 --output="out01000.txt" --error="err01000.txt" chan.sh 
# sbatch --export=CASE=01010 --output="out01010.txt" --error="err01010.txt" chan.sh 
# 
# sbatch --export=CASE=10000 --output="out10000.txt" --error="err10000.txt" chan.sh 
# sbatch --export=CASE=10010 --output="out10010.txt" --error="err10010.txt" chan.sh 
# sbatch --export=CASE=11000 --output="out11000.txt" --error="err11000.txt" chan.sh 
# sbatch --export=CASE=11010 --output="out11010.txt" --error="err11010.txt" chan.sh 
# 
# sbatch --export=CASE=00300 --output="out00300.txt" --error="err00300.txt" chan.sh 
# sbatch --export=CASE=00310 --output="out00310.txt" --error="err00310.txt" chan.sh 
# sbatch --export=CASE=01300 --output="out01300.txt" --error="err01300.txt" chan.sh 
# sbatch --export=CASE=01310 --output="out01310.txt" --error="err01310.txt" chan.sh 
