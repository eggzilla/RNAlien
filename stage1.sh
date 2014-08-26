#Rfam Statistics - auxiliary script - computes mlocarna alignments for all Rfam families
#!/bin/bash
#$ -t 1-2207 #This will start the job for each Rfam family
#$ -pe para 6
#$ -o /scr/kronos/egg/temp/RfamStat/output  
#$ -e /scr/kronos/egg/temp/RfamStat/error
#$ -l hostname="archer|tc00|tc01|tc02|tc03|tc04"
#$ -N Rfamstat-stage1
#compute structural alignments with locarna
/home/mescalin/egg/Tools/bin/mlocarna --threads=6 --local-progressive /scr/kronos/egg/temp/RfamStat2/$SGE_TASK_ID.fa > /scr/kronos/egg/temp/RfamStat2/$SGE_TASK_ID.mlocarna
