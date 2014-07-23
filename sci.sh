#!/bin/bash

#$ -t 1-2207 #This will start the job programname 2000 times with an incrementing SGE_TASK-ID for every run, starting with 1 until 2000
#$ -pe para 6
#$ -e /scr/kronos/egg/temp/RfamStat2/error
#$ -l hostname="archer|tc00|tc01|tc02|tc03|tc04"
#$ -N area51
/home/mescalin/egg/Tools/bin/mlocarna --threads=7 --local-progressive /scr/kronos/egg/temp/RfamStat2/$SGE_TASK_ID.fa > /scr/kronos/egg/temp/RfamStat2/$SGE_TASK_ID.mlocarna
/home/mescalin/egg/exec/bin/RNAz /scr/kronos/egg/temp/RfamStat2/$SGE_TASK_ID.mlocarna /scr/kronos/egg/temp/RfamStat2/$SGE_TASK_ID.rnazmlocarna
