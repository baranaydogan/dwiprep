#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:56:04 2021

@author: BARAN AYDOGAN
"""

def main():
    import os
    import sys
    from glob import glob
    import subprocess
    from aux import nifti_extractB0s
    from verifyInput import run_verifyInput
    from parseInput import run_parseInput,run_mergePhaseEncodedDwi
    from synb0Functions import synb0ifneeded
    from denoiseWithMPPCA import run_denoiseWithMPPCA    

    inpDir 	= sys.argv[1];
    T1path      = sys.argv[2];
    outDir 	= sys.argv[3];
    
    # SET THESE PATHS
    fsldir      = os.environ['FSLDIR'];
    mrtrix      = '/scratch/work/aydogad1/tools/build/mrtrix3/bin';
    synb0       = '/scratch/work/aydogad1/tools/barans/Synb0-DISCO_scripts';
    dwiprep     = '/scratch/work/aydogad1/tools/barans/dwiprep';             # Path to this script 
    
    if not os.path.exists(outDir): os.makedirs(outDir)
    if not os.path.exists(outDir + '/Step0_parsedInputImages'): os.mkdir(outDir + '/Step0_parsedInputImages')
    if not os.path.exists(outDir + '/Step1_beforeTopupAndEddy'): os.mkdir(outDir + '/Step1_beforeTopupAndEddy')
    if not os.path.exists(outDir + '/Step2_topupAndEddyInputs'): os.mkdir(outDir + '/Step2_topupAndEddyInputs')
    if not os.path.exists(outDir + '/Step3_topupOutput'): os.mkdir(outDir + '/Step3_topupOutput')
    if not os.path.exists(outDir + '/Step4_eddyOutput'): os.mkdir(outDir + '/Step4_eddyOutput')
    if not os.path.exists(outDir + '/Step5_preprocessedDMRI'): os.mkdir(outDir + '/Step5_preprocessedDMRI')
    if not os.path.exists(outDir + '/Step6_microstructureAnalysis'): os.mkdir(outDir + '/Step6_microstructureAnalysis')
    if not os.path.exists(outDir + '/Step7_registerDMRI2T1'): os.mkdir(outDir + '/Step7_registerDMRI2T1')

    report      = open((outDir + '/dwiprep_report.json'),'w');
    report.write('{\n');
    report.write('"inpDir": "' + inpDir + '"');
    report.write(',\n"outDir": "' + outDir + '"');
      
    prepEnv                     = os.environ.copy();
    prepEnv["PATH"]             = mrtrix + ":" + fsldir + ":" + synb0 + ":" + dwiprep + ":" + prepEnv["PATH"];
    
    
    # Step1 - Prepare for topup and eddy
    
    # Verify and parse input
    # - seperate if different phase encoded images exists in the same nifti
    # - merge all same phase encoded images in one volume
    
    dwis     = glob((inpDir+'/*.nii*'));

    RL_dwis  = [];
    LR_dwis  = [];
    AP_dwis  = [];
    PA_dwis  = [];
    SI_dwis  = [];
    IS_dwis  = [];
    
    reference = [];
    imageNo   = 1;
    for dwi in dwis:
        
        if (dwi[-3:]=='.gz'):
            basename    = os.path.basename(dwi[0:-7]);
        else:
            basename    = os.path.basename(dwi[0:-4]);
            
        basenameInp = inpDir + '/' + basename;
        basenameOut = outDir + '/Step0_parsedInputImages/' + basename;
        
        sys.stdout.write('Verifying and parsing input image: ' + basenameInp + '\n');
        
        # Verify input
        run_verifyInput(dwi,basenameOut,reference);
        report.write(',\n"Image ' + str(imageNo) + ' verified": "' + dwi + '"');
        
        if not reference : reference=(basenameOut + '.json');
        
        # Parse input image, i.e. seperate if multiple phase encodings exists
        (_RL_dwis,_LR_dwis,_PA_dwis,_AP_dwis,_IS_dwis,_SI_dwis) = run_parseInput(dwi,basenameOut);
        report.write(',\n"Image ' + str(imageNo) + ' parsed": "' + dwi + '"');
        
        RL_dwis  =  RL_dwis +  _RL_dwis;
        LR_dwis  =  LR_dwis +  _LR_dwis;
        PA_dwis  =  PA_dwis +  _PA_dwis;
        AP_dwis  =  AP_dwis +  _AP_dwis;
        IS_dwis  =  IS_dwis +  _IS_dwis;
        SI_dwis  =  SI_dwis +  _SI_dwis;
        
        imageNo = imageNo + 1;

    # Merge all dwi with same phase encoding
    dwis = run_mergePhaseEncodedDwi(RL_dwis,LR_dwis,PA_dwis,AP_dwis,IS_dwis,SI_dwis,(outDir + '/Step0_parsedInputImages'));

    # Use Synb0 if a reverse phase encoded b0 image is not available
    (RL_dwis,LR_dwis,PA_dwis,AP_dwis,IS_dwis,SI_dwis) = synb0ifneeded(prepEnv,T1path,RL_dwis,LR_dwis,PA_dwis,AP_dwis,IS_dwis,SI_dwis,(outDir + '/Step0_parsedInputImages'));
    dwis = RL_dwis + LR_dwis + PA_dwis + AP_dwis + IS_dwis + SI_dwis 

    # Report what is found
    if RL_dwis: report.write(',\n"Found RL image": "' + (outDir + '/Step0_parsedInputImages/RL.nii.gz') + '"');
    if LR_dwis: report.write(',\n"Found LR image": "' + (outDir + '/Step0_parsedInputImages/LR.nii.gz') + '"');
    if PA_dwis: report.write(',\n"Found PA image": "' + (outDir + '/Step0_parsedInputImages/PA.nii.gz') + '"');
    if AP_dwis: report.write(',\n"Found AP image": "' + (outDir + '/Step0_parsedInputImages/AP.nii.gz') + '"');
    if IS_dwis: report.write(',\n"Found IS image": "' + (outDir + '/Step0_parsedInputImages/IS.nii.gz') + '"');
    if SI_dwis: report.write(',\n"Found SI image": "' + (outDir + '/Step0_parsedInputImages/SI.nii.gz') + '"');

    # Preprocess
    name_base   = [];
    name_b0s    = [];
    name_dwi    = [];
    
    for dwi in dwis:
        
        if (dwi[-3:]=='.gz'):
            basename    = os.path.basename(dwi[0:-7]);
        else:
            basename    = os.path.basename(dwi[0:-4]);
            
        basenameInp = os.path.dirname(dwi) + '/' + basename;
        basenameOut = outDir + '/Step1_beforeTopupAndEddy/' + basename;
        
        name_base.append(basenameInp);
        
        sys.stdout.write('Processing: ' + basenameInp + '\n');
        
        report.write(',\n"Processing ' + basename + '": {');
        
        # Denoise dwi data with MPPCA
        out = dwi;

        report.write('\n"denoise with MPPCA": "true"');
        report.write(',\n"MPPCA info": { \n');
        denoiseResult=run_denoiseWithMPPCA(prepEnv, '5,5,5',out,basenameOut,basenameInp,report);
        report.write('\n}');
        if (denoiseResult==True):
            out=(basenameOut + '_dwidenoise.nii.gz');

        name_dwi.append(out);
        out=nifti_extractB0s(out,basenameInp,basenameOut,50,report);
        out=(basenameOut + '_b0.nii.gz');
        name_b0s.append(out);
        
        report.write('\n}');
        

    # Step2 - Apply topup and eddy
    from prepareTopupAndEddyInputs import run_prepareTopupAndEddyInputs
    run_prepareTopupAndEddyInputs(prepEnv,outDir,name_base,name_b0s,name_dwi,True);
    
    subprocess.call(['sbatch', '--wait', (outDir + '/Step2_topupAndEddyInputs/run_topup.sh')], env=prepEnv);
    subprocess.call(['sbatch', '--wait', (outDir + '/Step2_topupAndEddyInputs/run_eddy.sh')],  env=prepEnv);
    
    subprocess.call(['sbatch', '--wait', 'run_slurm_upsample.sh',outDir],  env=prepEnv);
    
    report.write(',\n"result": "SUCCESS"');
    report.write('\n}');
    report.close();
    
    
    
    
if __name__== "__main__":
  main();


    
    


