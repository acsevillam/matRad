matRad_DIR=../
FileName=$matRad_DIR/examples_interval/matRad_breastPhoton.m
matLab_DIR=/Applications/MATLAB_R2021a.app/bin

## se corren los programas 
${matLab_DIR}/matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('${matRad_DIR}'),'-end');run('${FileName}');exit;"