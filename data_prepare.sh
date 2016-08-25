#!/bin/bash
#$ -pe omp 8 
#$ -l h_rt=24:00:00
#$ -N prepare
#$ -V
printf "JOB_NAME: $JOB_NAME\n"
printf "JOB_ID: $JOB_ID\n" 
ML="/usr/local/bin/matlab -nodisplay -r "
$ML "addpath('~/ccdc');task=1;ntask=7;data_prep" &
$ML "addpath('~/ccdc');task=2;ntask=7;data_prep" &
$ML "addpath('~/ccdc');task=3;ntask=7;data_prep" &
$ML "addpath('~/ccdc');task=4;ntask=7;data_prep" &
$ML "addpath('~/ccdc');task=5;ntask=7;data_prep" &
$ML "addpath('~/ccdc');task=6;ntask=7;data_prep" &
$ML "addpath('~/ccdc');task=7;ntask=7;data_prep" &
wait