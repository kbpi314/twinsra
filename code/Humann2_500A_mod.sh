#!/bin/bash

LOGDIR=/sc/arion/projects/clemej05a/kevin/twinsra/outputs/ensemble_humann3/twinsra_h3/tool_logging/Humann2
OUTPUTDIR=/sc/arion/projects/clemej05a/kevin/twinsra/outputs/ensemble_humann3/twinsra_h3/tool_output/Humann2/500A
JOB=humann2_500A
JOBFILE=/sc/arion/projects/clemej05a/kevin/twinsra/outputs/ensemble_humann3/twinsra_h3/tool_scripts/Humann2/jobs/${JOB}.lsf

cat <<EOF > $JOBFILE
#!/bin/bash
#BSUB -P acc_clemej05a
#BSUB -n 1
#BSUB -q premium
#BSUB -W 24:00
#BSUB -J $JOB
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=30000]
#BSUB -o $LOGDIR/%J.stdout
#BSUB -eo $LOGDIR/%J.stderr
#BSUB -L /bin/bash

# new
mkdir /sc/arion/projects/clemej05a/kevin/twinsra/outputs/ensemble_humann3/twinsra_h3/temp_data/Humann2/500A/

# change last path to Humann2/500A/500A.fastq
cat /sc/arion/projects/clemej05a/PsA_IL17/inputs/Metagenomes/combined/twinsra_kneaddata/500A_R1.fastq /sc/arion/projects/clemej05a/PsA_IL17/inputs/Metagenomes/combined/twinsra_kneaddata/500A_R2.fastq > /sc/arion/projects/clemej05a/kevin/twinsra/outputs/ensemble_humann3/twinsra_h3/temp_data/Humann2/500A/500A.fastq

module purge
ml anaconda3
source activate /sc/arion/projects/clemej05a/adam/conda/envs/humann3 

# modify input fastq path
humann --input /sc/arion/projects/clemej05a/kevin/twinsra/outputs/ensemble_humann3/twinsra_h3/temp_data/Humann2/500A/500A.fastq --output $OUTPUTDIR  --metaphlan-options "-t rel_ab_w_read_stats --bowtie2db /sc/arion/projects/CVDlung/databases/ensemble-metaphlan" --nucleotide-database /sc/arion/projects/clemej05a/adam/databases/humann3/chocophlan --protein-database /sc/arion/projects/clemej05a/adam/databases/humann3/uniref

# delete ONLY dir pertaining to sample
rm -rf /sc/arion/projects/clemej05a/kevin/twinsra/outputs/ensemble_humann3/twinsra_h3/temp_data/Humann2/500A

EOF
bsub < $JOBFILE
