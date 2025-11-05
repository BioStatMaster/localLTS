#!/bin/sh
jobName=run_exps
#for jobId in $(seq 1 5)
for jobId in 1
do
    logDir=/pghbio/dbmi/batmanlab/mgong1/${jobName}
    if [ ! -d "${logDir}" ]; then
        mkdir ${logDir}
    fi
    outFile=${logDir}/output_${jobName}${jobId}.log
    errorFile=${logDir}/error_${jobName}${jobId}.log
    echo submit ${jobName}${jobId}
    echo ${outFile}
    echo ${errorFile}
    sbatch -A bi561ip -p DBMI -n 10 -o ${outFile} -e ${errorFile} -t 24:00:00  ${jobName}.sh
    #sbatch -p GPU-shared --gres=gpu:p100:1  -o ${outFile} -e ${errorFile} -t 48:00:00  ${jobName}.sh
done
squeue -u mgong1
