#Alien Benchmark
#!/bin/bash
#$ -t 1-56 #This will start the job for each sRNA Rfam family
#$ -l mem_free=10G
#$ -j yes
#$ -o /scratch/egg/temp/
#$ -e /scratch/egg/temp/
#$ -l hostname="tc00|tc01|tc02|tc03|tc04"
#$ -N area54
#alienrun
if [ -f /scr/kronos/egg/AlienStructuredResultsCollected4/$SGE_TASK_ID/done ]; then
		cmComparevsRfam.pl $SGE_TASK_ID 
		sleep 1
		echo "File not found!"
fi
