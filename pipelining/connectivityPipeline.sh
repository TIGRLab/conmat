#!/bin/bash

. $MODULESHOME/init/bash
module load camino/1886d5
module load FSL/5.0.10
module load python

python /home/lliu/Desktop/conmat/pipelining/connectivityPipeline.py SPINS
