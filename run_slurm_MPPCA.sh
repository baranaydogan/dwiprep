#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=45
#SBATCH --output SlurmFiles/Denoising_output_%A.txt
#SBATCH --error SlurmFiles/Denoising_error_%A.txt
##SBATCH -o /dev/null

set -euo pipefail


inp=$1
basenameout=$2
extent=$3

mrtrix=$4 #/scratch/work/aydogad1/tools/build/mrtrix3/bin

${mrtrix}/dwidenoise \
${inp} \
${basenameout}_dwidenoise.nii.gz \
-noise ${basenameout}_dwidenoise_noise.nii.gz \
-extent ${extent} \
-force

