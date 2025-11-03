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

export JOB=10000; sbatch -J j${JOB} --export=CASE=${JOB} --output="out${JOB}.txt" --error="err${JOB}.txt" chan.sh 
export JOB=10010; sbatch -J j${JOB} --export=CASE=${JOB} --output="out${JOB}.txt" --error="err${JOB}.txt" chan.sh 
export JOB=10030; sbatch -J j${JOB} --export=CASE=${JOB} --output="out${JOB}.txt" --error="err${JOB}.txt" chan.sh 

