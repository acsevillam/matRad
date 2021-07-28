#!/bin/sh
#
# job1.sh run matRad with test macro via APOLO 
# Copyright Andrés Camilo Sevilla
#
# !!!! Must be changed for each user !!! #

# Defining variables
matLab_DIR=/Applications/MATLAB_R2021a.app/bin
matRad_DIR=/Users/asevilla/workspace/matRad/dev_contourPropagation
JOB_DIR=matRad_DIR/jobs

DataNo=1
i=0

# Creating DataNo directory

while [ $i -lt 1 ]
do
	if [ -d "./${DataNo}" ]
	then
		echo " “./${DataNo}” has been created"
		DataNo=$(($DataNo+1))
		i=0
	else
		mkdir ./${DataNo}
		i=1
fi
done

# Creating log directory

if [ -d "./${DataNo}/log" ]
then
	echo " “./${DataNo}/log” has been created"
else
	mkdir ./${DataNo}/log
fi

# Creating data directory

if [ -d "./${DataNo}/data" ]
then
	echo " “./${DataNo}/data” has been created"
else
	mkdir ./${DataNo}/data
fi

# simple script to run matlab script
if [ $# -eq 0 ]
  then
    echo "please pass m script"
fi

$matLab_DIR/matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('${matRad_DIR}'),'-end');run('${1}');exit;"

cd $JOB_DIR/