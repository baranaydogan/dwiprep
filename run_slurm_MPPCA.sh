#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=45
#SBATCH -o /dev/null

set -euo pipefail


inp=$1
basenameout=$2
extent=$3

mrtrix=/scratch/work/aydogad1/tools/build/mrtrix3/bin

${mrtrix}/dwidenoise \
${inp} \
${basenameout}_dwidenoise.nii.gz \
-noise ${basenameout}_dwidenoise_noise.nii.gz \
-extent ${extent} \
-force

