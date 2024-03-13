#!/bin/bash
#SBATCH --nodes=3
#SBATCH --ntasks=300
#SBATCH --mem=300000
#SBATCH --time=1-00:00:00
#SBATCH --partition=regular
#SBATCH --account=hy299_0001
#SBATCH --chdir=/home/yz2296/slurm/
#SBATCH --job-name=week8
#SBATCH --output=/home/yz2296/slurm/logs/week8.%j
#SBATCH --mail-user=yz2296@cornell.edu
#SBATCH --mail-type=ALL
set -o errexit;

# RUN THIS SCRIPT: sbatch submit.sh ${job_name}
# CHECK JOB STATUS: squeue -u yz2296
# CANCEL JOB: scancel 1564
# If the job duration is longer than 24 hours, use --partition=long7. Otherwise, use --partition=regular

# # init conda env
# source /home/yc2553/miniconda3_new/etc/profile.d/conda.sh;
# conda activate cy;

# create tmp folder for this job
mkdir -p /workdir/$USER/$SLURM_JOB_ID;
cd /workdir/$USER/$SLURM_JOB_ID;

# # mount drives
/programs/bin/labutils/mount_server cbsuhy01 /storage
/programs/bin/labutils/mount_server cbsuhyfs1 /storage
/programs/bin/labutils/mount_server cbsuhyfs1 /storage1
/programs/bin/labutils/mount_server cbsuhy02 /storage

# source /home/yc2553/slurm/_load_NFS.sh;

# print parameters
echo $@;

       
# permanent_folder=$1;      # move results to this folder
# job_name=$1; # job name
# pl_bw=$1;                  # path to the input bam (pl)
# mn_bw=$2;                  # path to the input bam (mn)
# cp ${pl_bw} .;
# cp ${mn_bw} .;

# python /fs/cbsuhy01/storage/yz2296/cancer_hotspot_new/website/website/revision/revision.py
# python /fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/week8_localrun.py -o ./
# python /fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/week8_localrun.py 
python /fs/cbsuhy01/storage/yz2296/cancer_hotspot_new/website/website/revision/pval_calibration.py
# if [ -z ${5+x} ];
# then
# 	fdr_target=0.1;
# else
# 	fdr_target=$5;
# fi;

# pints_caller --bw-pl "$(basename -- $pl_bw)" \
# 	--bw-mn "$(basename -- $mn_bw)" \
# 	--save-to . --file-prefix ${job_name} \
# 	--thread 16 --fdr-target ${fdr_target} \
# 	--min-lengths-opposite-peaks 5;


# if [ ! -d ${permanent_folder} ];
# then
# 	mkdir ${permanent_folder};
# fi;
# # rm "$(basename -- $pl_bw)";
# # rm "$(basename -- $mn_bw)";
# mv * ${permanent_folder};
