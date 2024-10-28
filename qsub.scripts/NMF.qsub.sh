## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -r yes
#$ -l mem_requested=20G
#$ -pe smp 8
#$ -M lacgra@garvan.org.au
#$ -m ae
#$ -N out.NMF_UC
#$ -t 1-28

conda activate Aim_1

#study_list=(UC_GSE125527 CD_Kong/colon CD_Kong/TI lupus_Chun MS_GSE193770)
# study=${study_list[$SGE_TASK_ID-1]}

study=UC_GSE182270

cd /directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/$study

Rscript /directflow/SCCGGroupShare/projects/lacgra/PhD/functions/NMF.R pbmc.female.RDS $SGE_TASK_ID