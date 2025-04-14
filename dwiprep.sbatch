#!/bin/bash
#SBATCH --time 2:00:00              # Runtime in minutes.
#SBATCH --mem 8G                       # Reserve 32 GB RAM for the job
#SBATCH --output SlurmFiles/Preprocessing_output_%A_%a.txt
#SBATCH --error SlurmFiles/Preprocessing_error_%A_%a.txt
#SBATCH --nodes=1               # Request one node
#SBATCH --array=0-114

module load fsl/6.0.7.11
source $FSLDIR/etc/fslconf/fsl.sh
source activate /scratch/nbe/istim/Simona/code/NeuroEnv

PARENT_DIR=/scratch/nbe/istim/Simona/data/BIDS_scaffold

# Get the list of all subject folders
FOLDERS=($(find $PARENT_DIR -maxdepth 1 -name "sub-R*" -type d -printf '%f\n' | sort ))

# Select the folder for this array task
FOLDER=${FOLDERS[$SLURM_ARRAY_TASK_ID]}

srun python dwiprep.py $PARENT_DIR/$FOLDER/dwi $PARENT_DIR/$FOLDER/anat/${FOLDER}_T1w.nii.gz $PARENT_DIR/derivatives/$FOLDER/dwi_preproc $FSLDIR
