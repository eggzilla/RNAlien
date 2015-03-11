#Alien Benchmark
#!/bin/bash
#$ -t 1-373 #This will start the job for each sRNA Rfam family
#$ -pe para 7
#$ -l mem_free=34.9G
#$ -j yes
#$ -o /scratch/egg/temp/
#$ -e /scratch/egg/temp/
#$ -l hostname="xc00|xc01|xc02|xc03|xc04|xc05|xc06|xc07|xc08"
#$ -N area54
#alienrun
if [ ! -f /scratch/egg/AlienResultsCollected/$SGE_TASK_ID/done ]; then
	/home/mescalin/egg/current/Projects/Haskell/RNAlien/dist/build/RNAlien/RNAlien -i /scr/kronos/egg/AliensRNATestSet/$SGE_TASK_ID.fa -c 7 -t "$(</scr/kronos/egg/AliensRNATestSet/$SGE_TASK_ID.tax)" -n /scratch/egg/nodes.dmp -o /scratch/egg/temp/ -d $SGE_TASK_ID > /scratch/egg/temp/$SGE_TASK_ID.alienout
	sleep 1
	echo "File not found!"
fi
