#Rfam Statistics - auxiliary script - computes clustalw2 alignments and RNAz results for all Rfam families
#!/bin/bash
#$ -t 1-2207 #This will start the job for each Rfam family
#$ -o /scr/kronos/egg/temp/RfamStat/output  
#$ -e /scr/kronos/egg/temp/RfamStat/error
#$ -l hostname="archer|tc00|tc01|tc02|tc03|tc04"
#$ -N Rfamstat-stage2
#compute alignments with locarna
clustalw2 -INFILE=/scr/kronos/egg/temp/RfamStat2/$SGE_TASK_ID.fa -OUTFILE=/scr/kronos/egg/temp/RfamStat/$SGE_TASK_ID.clustalw2 > /scr/kronos/egg/temp/RfamStat/$SGE_TASK_ID.alnsum
#Run RNAz for clustalw2 alignment results
/home/mescalin/egg/src/exec/bin/RNAz /scr/kronos/egg/temp/RfamStat/$SGE_TASK_ID.clustalw2 > /scr/kronos/egg/temp/RfamStat/$SGE_TASK_ID.rnaz
#run RNAz for locarna alignment results in clustalw2 format
/home/mescalin/egg/src/exec/bin/RNAz /scr/kronos/egg/temp/RfamStat2/$SGE_TASK_ID.out/results/result.aln  /scr/kronos/egg/temp/RfamStat2/$SGE_TASK_ID.rnazmlocarna
