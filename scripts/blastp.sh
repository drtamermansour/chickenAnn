#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=4,mem=32Gb
#mdiag -A ged
#PBS -m abe
#PBS -N blastp

module load BLAST+/2.2.29

cd $PBS_O_WORKDIR

blastp -query "$input" -db "$DB" -num_threads 4 -max_target_seqs 20 -outfmt 5 -seg yes -evalue 1e-3 > "$input"."$label".blastp.xml

qstat -f ${PBS_JOBID}
