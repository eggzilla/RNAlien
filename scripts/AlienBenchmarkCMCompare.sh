#Alien Benchmark
#!/bin/bash
#$ -t 1-373 #This will start the job for each sRNA Rfam family
#$ -l mem_free=10G
#$ -j yes
#$ -o /scratch/egg/temp/
#$ -e /scratch/egg/temp/
#$ -l hostname="xc00|xc01|xc02|xc03|xc04|xc05|xc06|xc07|xc08"
#$ -N area54
#alienrun
if [ -f /scr/kronos/egg/AlienResultsCollected/$SGE_TASK_ID/done ]; then
		cmComparevsRfam.pl $SGE_TASK_ID 
		sleep 1
		echo "File not found!"
fi
