#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:32:15 2018

@author: baran
"""

# Denoise dwi data with MPPCA
# mrtrix3's dwidenoise command is used for this
def run_denoiseWithMPPCA(prepEnv, mrtrix_dwinoise_extent,inp,basenameOut,basenameInp,mrtrix,report):
    
    import sys
    import json
    import subprocess
    from glob import glob
    from aux import nifti_subtract   
    import nibabel as nib 
    import numpy as np

    fid     = open((basenameInp + '.json'));
    inpinfo = json.load(fid);
    fid.close();
    
    if (inpinfo['number of volumes']<2):
       sys.stdout.write('   Skipped denoising since number of volumes = 1...\n');
       report.write('"result": "skipped denoising since number of volumes = 1"');
       nii  = nib.load(inp);
       img  = nii.get_data();
       if (len(img.shape)==3):
           img = np.expand_dims(img,3);
       if (nii.header.sizeof_hdr==348):
          nio   = nib.Nifti1Image(img,nii.affine);
       else:
          nio   = nib.Nifti2Image(img,nii.affine);

       nib.save(nio,(basenameOut + '_dwidenoise.nii.gz'));
       return True;
    
    sys.stdout.write('   Denoising... ');
    
    report.write('"software": "mrtrix, '+ str(subprocess.check_output(['dwidenoise', '--version'],env=prepEnv).splitlines()[0])[2:-1] + '"');
    report.write(',\n"parameters": { \n');
    report.write('"-extent": "' + mrtrix_dwinoise_extent +'"');
    report.write('\n}');
    
    subprocess.call(['sbatch', '--wait', 'run_slurm_MPPCA.sh',inp,basenameOut,mrtrix_dwinoise_extent, mrtrix], env=prepEnv);
    
    # subprocess.call(['dwidenoise', \
    #                  inp, \
    #                  (basenameOut + '_dwidenoise.nii.gz'), \
    #                  '-noise', \
    #                  (basenameOut + '_dwidenoise_noise.nii.gz'), \
    #                  '-extent', \
    #                  mrtrix_dwinoise_extent, \
    #                  '-force'], \
    #                  env=prepEnv);
        
    
    if not glob((basenameOut + '_dwidenoise.nii.gz')):
        report.write('\n}');
        report.write(',\n"result": "FAILED, denoising failed"');
        report.write('\n}');
        report.close();
        sys.exit('denoising failed');            
    if not glob((basenameOut + '_dwidenoise_noise.nii.gz')):
        report.write('\n}');
        report.write(',\n"result": "FAILED, denoising failed"');
        report.write('\n}');
        report.close();
        sys.exit('denoising failed');        
    sys.stdout.write('DONE \n');

    # Setting sub-zero values to zero
    subprocess.call(['fslmaths', \
                     (basenameOut + '_dwidenoise.nii.gz'), \
                     '-thr', \
                     '0', \
                     (basenameOut + '_dwidenoise.nii.gz')], \
                     env=prepEnv);
    
    # Calculate residual between the input image and the denoised image using fslmaths
    sys.stdout.write('   Writing residual... ');
    nifti_subtract(inp,(basenameOut + '_dwidenoise.nii.gz'),(basenameOut + '_dwidenoise_residual.nii.gz'));    
    
    if not glob((basenameOut + '_dwidenoise_residual.nii.gz')):
        report.write('\n}');
        report.write(',\n"result": "FAILED, cannot calculate denoising residual"');
        report.write('\n}');
        report.close();
        sys.exit('cannot calculate denoising residual');        
    sys.stdout.write('DONE (please visually check the residual) \n');
    
    out     = (basenameOut + '_dwidenoise.nii.gz');
    
    report.write(',\n"output": "'+ out +'"');
    report.write(',\n"residual": "'+ (basenameOut + '_dwidenoise_residual.nii.gz') +'"');
    report.write(',\n"result": "success, please visually check the denoised output and residual"');
        
    return True;
