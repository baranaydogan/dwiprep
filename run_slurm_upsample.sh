#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=45
#SBATCH -o SlurmFiles/upsample_output_%A.txt
#SBATCH -e SlurmFiles/upsample_error_%A.txtx

set -euo pipefail

outdir=$1

mrtrix= $2#/scratch/work/aydogad1/tools/build/mrtrix3/bin

cp ${outdir}/Step4_eddyOutput/prepped_dMRI.bval ${outdir}/Step5_preprocessedDMRI/dMRI.bval
cp ${outdir}/Step4_eddyOutput/prepped_dMRI.bvec ${outdir}/Step5_preprocessedDMRI/dMRI.bvec

${mrtrix}/mrgrid -force ${outdir}/Step4_eddyOutput/prepped_dMRI.nii.gz -force regrid -voxel 1.25 -datatype float32 ${outdir}/Step5_preprocessedDMRI/dMRI.nii.gz
${mrtrix}/mrgrid -force ${outdir}/Step4_eddyOutput/prepped_dMRI_mask.nii.gz -force regrid -voxel 1.25 -datatype uint16 ${outdir}/Step5_preprocessedDMRI/dMRI_mask.nii.gz
${mrtrix}/mrcalc -force ${outdir}/Step5_preprocessedDMRI/dMRI.nii.gz ${outdir}/Step5_preprocessedDMRI/dMRI_mask.nii.gz -multiply ${outdir}/Step5_preprocessedDMRI/dMRI.nii.gz

dMRI=${outdir}/Step5_preprocessedDMRI/dMRI

# Filter implausible values
${mrtrix}/mrconvert -force -fslgrad ${dMRI}.bvec ${dMRI}.bval ${dMRI}.nii.gz ${dMRI}.mif
${mrtrix}/dwiextract -force -bzero ${dMRI}.mif - | ${mrtrix}/mrmath -force - median ${dMRI}_b0.nii.gz -axis 3
${mrtrix}/mrcalc -force ${dMRI}_b0.nii.gz 0.001 -lt 0.001 ${dMRI}_b0.nii.gz -if ${dMRI}_b0.nii.gz
${mrtrix}/dwiextract -force -no_bzero ${dMRI}.mif ${dMRI}_nob0.mif
${mrtrix}/mrcalc -force ${dMRI}_nob0.mif ${dMRI}_b0.nii.gz -gt ${dMRI}_b0.nii.gz ${dMRI}_nob0.mif -if ${dMRI}_nob0.mif
${mrtrix}/mrconvert -force -export_grad_fsl ${dMRI}_nob0.bvec ${dMRI}_nob0.bval ${dMRI}_nob0.mif ${dMRI}_nob0.nii.gz
${mrtrix}/mrcat -axis 3 -force ${dMRI}_b0.nii.gz ${dMRI}_nob0.nii.gz ${dMRI}.nii.gz

echo 0 `sed '1q;d' ${dMRI}_nob0.bval` >  ${dMRI}.bval
echo 1 `sed '1q;d' ${dMRI}_nob0.bvec` >  ${dMRI}.bvec
echo 0 `sed '2q;d' ${dMRI}_nob0.bvec` >> ${dMRI}.bvec
echo 0 `sed '3q;d' ${dMRI}_nob0.bvec` >> ${dMRI}.bvec

rm ${dMRI}_b0.*
rm ${dMRI}_nob0.*
rm ${dMRI}.mif

echo "Done"

