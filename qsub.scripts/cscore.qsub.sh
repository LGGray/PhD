## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -r yes
#$ -l mem_requested=20G
#$ -pe smp 8
#$ -M lacgra@garvan.org.au
#$ -m ae
#$ -N out.cscore.pSS
#$ -j 1:20

#conda activate Aim_1
conda activate cscore

#study_list=(UC_GSE125527 CD_Kong/colon CD_Kong/TI lupus_Chun MS_GSE193770)
# study=${study_list[$SGE_TASK_ID-1]}

study=pSS_GSE157278

cd /directflow/SCCGGroupShare/projects/lacgra/autoimmune.datasets/$study

mkdir -p cscore

if [[ $study == "lupus_Chun" ]]; then
  Rscript /directflow/SCCGGroupShare/projects/lacgra/PhD/functions/cscore.R pbmc.female.control-managed.RDS $SGE_TASK_ID
else
  Rscript /directflow/SCCGGroupShare/projects/lacgra/PhD/functions/cscore.R pbmc.female.RDS $SGE_TASK_ID
fi