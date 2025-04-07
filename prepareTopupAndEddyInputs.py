#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 13:29:08 2018

@author: baran
"""

def run_prepareTopupAndEddyInputs(prepEnv,outPath,name_base,name_b0s,name_dwi,fsldir,average_bvecs_path, doTopup):
    
    import os
    import sys
    import numpy as np
    from aux import nifti_concatenateImages,nifti_extractVolumes,nifti_extract_a_b0
    import json

    # Combine all b0 images
    sys.stdout.write('Combining: \n');
    for i in range(len(name_b0s)):
        sys.stdout.write(str(i+1) + '. ' + name_b0s[i] + '\n');
    nv_b0s          = nifti_concatenateImages((outPath+'/Step2_topupAndEddyInputs/combined_b0s.nii.gz'),name_b0s);
    
    # Combine all dwi images
    sys.stdout.write('Combining: \n');
    for i in range(len(name_dwi)):
        sys.stdout.write(str(i+1) + '. ' + name_dwi[i] + '\n');
    nv_dwi          = nifti_concatenateImages((outPath+'/Step2_topupAndEddyInputs/combined_dwi.nii.gz'),name_dwi);

    # Combine bvals and bvecs
    
    b0_threshold = 50;
    
    sys.stdout.write('Combining bvals and bvecs \n');

    combined_bvals = np.empty((0,),float);
    combined_bvecs = np.empty((3,0),float);

    blipVolCounts  = [0];
    combined_b0_inds=[];
    for i in range(len(name_base)):
        cur_bvals       = np.atleast_1d(np.loadtxt((name_base[i] + '.bval')));
        cur_bvecs       = np.atleast_2d(np.loadtxt((name_base[i] + '.bvec')));
        
        if (np.shape(cur_bvecs)[0]==1):
            cur_bvecs = cur_bvecs.transpose();

        combined_bvals  = np.concatenate((combined_bvals,cur_bvals),axis=0);
        combined_bvecs  = np.concatenate((combined_bvecs,cur_bvecs),axis=1);
        
        combined_b0_inds.append(np.nonzero(cur_bvals<b0_threshold));
        blipVolCounts.append(cur_bvals.shape[0]);
        
        np.savetxt((outPath + '/Step2_topupAndEddyInputs/' + os.path.basename(name_base[i]) + '_SeriesVolNum.txt'),[blipVolCounts[i+1],blipVolCounts[i+1]], delimiter=' ',   fmt='%d', newline=' ');
    
    np.savetxt((outPath+'/Step2_topupAndEddyInputs/combined_dwi.bval'), combined_bvals.astype(int), delimiter='',   fmt='%d', newline=' ');
    np.savetxt((outPath+'/Step2_topupAndEddyInputs/combined_dwi.bvec'), combined_bvecs.astype(float), delimiter=' ', fmt='%1.6f', newline='\n');

    # Prepare acqparams
    sys.stdout.write('Extracting index files \n');

    # Here a maximum of images can be set used to limit pairs of blup up-down b0 couples
    max_b0_count    = min(nv_b0s);
    if (max_b0_count>4):
       max_b0_count=4

    nvc_b0          = np.append([-1],np.cumsum(nv_b0s)-1);
    nvc_dwi         = np.append([0],np.cumsum(nv_dwi));
    nvs_b0          = np.zeros_like(nv_b0s);
    nvs_inds        = np.zeros((max_b0_count*nv_b0s.shape[0],),dtype=int);
    selected_b0_inds= np.zeros((max_b0_count*nv_b0s.shape[0],),dtype=int);
    for i in range(0,nv_b0s.shape[0]):
        nvs_b0[i]   = max_b0_count;
        nvs_inds[i*max_b0_count:(i+1)*max_b0_count,]=np.unique(np.round(np.linspace(nvc_b0[i]+1,nvc_b0[i+1],max_b0_count)));
        selected_b0_inds[i*max_b0_count:(i+1)*max_b0_count,]=combined_b0_inds[i][0][nvs_inds[i*max_b0_count:(i+1)*max_b0_count,]-nvs_inds[i*max_b0_count]]+nvc_dwi[i];

    if ( (selected_b0_inds.shape[0]==1) & (selected_b0_inds[0]==0) ):
        indices         = np.ones((np.sum(nv_dwi),),dtype=int);
    else:
        indices         = np.zeros((np.sum(nv_dwi),),dtype=int);
        indices[0:selected_b0_inds[0]]=selected_b0_inds[0];
        for i in np.arange(0,selected_b0_inds.shape[0]-1):
            indices[selected_b0_inds[i]:selected_b0_inds[i+1]]=i+1;
        indices[selected_b0_inds[i+1]:,]=i+2;
    
    indexFile       = open((outPath + '/Step2_topupAndEddyInputs/index.txt'),'w');
    for i in indices:
        indexFile.write(str(i) + '\n');
    indexFile.close();
#    BE CAREFUL WITH THIS ONE!    
    nifti_extractVolumes((outPath+'/Step2_topupAndEddyInputs/combined_b0s.nii.gz'),(outPath+'/Step2_topupAndEddyInputs/combined_b0s.nii.gz'),nvs_inds);


    # Prepare acqparams
    sys.stdout.write('Extracting acqparams.txt \n');
    
    acqparams       = open((outPath + '/Step2_topupAndEddyInputs/acqparams.txt'),'w');
    
    hifi_b0_fnames  = [];
    hifi_b0_inds    = [];
    indCounter      = 0;
    
    for i in range(len(name_base)):
    
        fid = open((name_base[i] + '.json'));
        dwiInfo = json.load(fid);
        fid.close();


        hifi_b0_fname = (outPath+'/Step2_topupAndEddyInputs/b0_at_ind_' + str(indCounter+1) + '.nii.gz');
        hifi_b0_ind   = indCounter;
        
        nifti_extractVolumes(hifi_b0_fname,(outPath+'/Step2_topupAndEddyInputs/combined_b0s.nii.gz'),hifi_b0_ind);
        
        hifi_b0_fnames.append(hifi_b0_fname);
        hifi_b0_inds.append(hifi_b0_ind+1);
        
        trt = str(float(dwiInfo['TotalReadoutTime']));
        

        try:
           phaseEncodingInfo=dwiInfo['PhaseEncodingDirection'];
        except:
           phaseEncodingInfo=dwiInfo['PhaseEncodingAxis'];

        
        if   (phaseEncodingInfo=='i' ): pe = '1 0 0 ';
        elif (phaseEncodingInfo=='i-'): pe = '-1 0 0 ';
        elif (phaseEncodingInfo=='j' ): pe = '0 1 0 ';
        elif (phaseEncodingInfo=='j-'): pe = '0 -1 0 ';
        elif (phaseEncodingInfo=='k' ): pe = '0 0 1 ';
        elif (phaseEncodingInfo=='k-'): pe = '0 0 -1 ';
        
        for j in np.arange(0,nvs_b0[i]):
            indCounter = indCounter + 1;
            acqparams.write(pe + trt + '\n');
            
    acqparams.close();
    
    # Prepare slspec
    slspec = np.array([]);
    try:
        sys.stdout.write('Extracting slspec.txt\n');
        mbfactor = int(dwiInfo['MultibandAccelerationFactor']);
        imgDims = dwiInfo['dcmmeta_shape'];
        numberOfExcitations = int(imgDims[2]/mbfactor);
        sliceTimings = np.asarray(dwiInfo['global']['slices']['CsaImage.TimeAfterStart'][0:imgDims[2]]);
        sortInd = np.argsort(sliceTimings);
        slspec=np.reshape(sortInd,[numberOfExcitations,mbfactor]);
        np.savetxt((outPath+'/Step2_topupAndEddyInputs/slspec.txt'), slspec.astype(int), delimiter=' ',   fmt='%d', newline='\n');
    except:
        slspec = np.array([]);
        sys.stdout.write('Eddy will not do slice to volume registration \n');
    
    
    # Prepare topup script
    sys.stdout.write('Compiling run_topup.sh script\n');
    
    topupScript   = open((outPath + '/Step2_topupAndEddyInputs/run_topup.sh'),'w');
    

    topupScript.write('#!/bin/bash \n');
    if (doTopup==True) :     topupScript.write('#SBATCH -N 1\n');
    if (doTopup==True) :     topupScript.write('#SBATCH --time=600 \n');
    if (doTopup==True) :     topupScript.write('#SBATCH --mem-per-cpu=2G \n');
    if (doTopup==True) :     topupScript.write('#SBATCH --error SlurmFiles/Topup_error_%A.txt \n');
    if (doTopup==True) :     topupScript.write('#SBATCH --output  SlurmFiles/Topup_output_%A.txt \n');

    if (doTopup==False) :    topupScript.write('#SBATCH -N 1\n');
    if (doTopup==False) :    topupScript.write('#SBATCH --time=30 \n');
    if (doTopup==False) :    topupScript.write('#SBATCH --mem-per-cpu=4G \n');
    if (doTopup==False) :     topupScript.write('#SBATCH --error  SlurmFiles/Topup_error_%A.txt \n');
    if (doTopup==False) :     topupScript.write('#SBATCH --output  SlurmFiles/Topup_output_%A.txt \n');

    fsl_dir="fsl"+fsldir.split("fsl",1)[1]
    topupScript.write('\n');
    topupScript.write('module load '+fsl_dir+' \n');
    topupScript.write('source $FSLDIR/etc/fslconf/fsl.sh\n');
    topupScript.write('\n');
    topupScript.write('\nderPath=' + outPath);
    topupScript.write('\n');
    topupScript.write('\n');
    topupScript.write('input=${derPath}/Step2_topupAndEddyInputs/combined_b0s.nii.gz' + '\n');
    topupScript.write('acqparams=${derPath}/Step2_topupAndEddyInputs/acqparams.txt' + '\n');
    topupScript.write('\n');
    topupScript.write('out=${derPath}/Step3_topupOutput/topup_results \n');
    topupScript.write('iout=${derPath}/Step3_topupOutput/topup_unwarped \n');
    topupScript.write('fout=${derPath}/Step3_topupOutput/topup_field \n');
    topupScript.write('\n');
    topupScript.write('hifib0=${derPath}/Step3_topupOutput/hifib0' + '\n');
    topupScript.write('nodif_brain=${derPath}/Step3_topupOutput/nodif_brain.nii.gz' + '\n');
    topupScript.write('nodif_brain_mask=${derPath}/Step3_topupOutput/nodif_brain_mask.nii.gz' + '\n');
    topupScript.write('nodif_brain_mask_dilated=${derPath}/Step3_topupOutput/nodif_brain_mask_dilated.nii.gz' + '\n');
    topupScript.write('\n');

    if (doTopup==True):
       topupScript.write('topup \\');
       topupScript.write('\n--verbose \\');
       topupScript.write('\n--imain=${input} \\');
       topupScript.write('\n--datain=${acqparams} \\');
       topupScript.write('\n--config=b02b0.cnf \\');
#    topupScript.write('\n--warpres=20,20,16,16,12,10,8,4 \\');
#    topupScript.write('\n--subsamp=8,8,4,4,2,2,2,2 \\');
#    topupScript.write('\n--fwhm=16,16,8,8,8,6,4,2 \\');
#    topupScript.write('\n--miter=4,4,4,4,5,5,5,5 \\');
#    topupScript.write('\n--lambda=0.5,0.1,0.05,0.01,0.005,0.001,0.0001,0.000015 \\');
#    topupScript.write('\n--ssqlambda=1 \\');
#    topupScript.write('\n--regmod=bending_energy \\');
#    topupScript.write('\n--estmov=1,1,1,1,1,1,1,0 \\');
#    topupScript.write('\n--minmet=0,0,0,0,0,0,0,0 \\');
#    topupScript.write('\n--splineorder=3 \\');
#    topupScript.write('\n--numprec=double \\');
#    topupScript.write('\n--interp=spline \\');
#    topupScript.write('\n--scale=0 \\');
       topupScript.write('\n--out=$out \\');
       topupScript.write('\n--iout=$iout \\');
       topupScript.write('\n--fout=$fout');
       topupScript.write('\n');
       topupScript.write('\n');
       topupScript.write('\n');
       topupScript.write('\n');
       topupScript.write('applytopup \\');
       topupScript.write('\n--verbose \\');
       topupScript.write('\n--imain=');
       for i in range(0,len(hifi_b0_fnames)):
           tmp = os.path.basename(hifi_b0_fnames[i]);
           topupScript.write('${derPath}/Step2_topupAndEddyInputs/' + tmp);
           if (i!=(len(hifi_b0_fnames)-1)): topupScript.write(',');
       topupScript.write(' \\');

       topupScript.write('\n--topup=$out \\');
       topupScript.write('\n--datain=${acqparams} \\');
    
       topupScript.write('\n--inindex=');
       for i in range(0,len(hifi_b0_inds)):
           topupScript.write(str(hifi_b0_inds[i]));
           if (i!=(len(hifi_b0_inds)-1)): topupScript.write(',');
       topupScript.write(' \\');
       topupScript.write('\n--out=${hifib0}');
       topupScript.write('\n');
       topupScript.write('\n');
       topupScript.write('\n');
       topupScript.write('\n');
    else:
       nifti_extract_a_b0((outPath+'/Step2_topupAndEddyInputs/combined_dwi.nii.gz'),(outPath+'/Step2_topupAndEddyInputs/combined_dwi'),(outPath + '/Step3_topupOutput/hifib0'),50)
    
    topupScript.write('bet $hifib0 ${nodif_brain} -m -R -f 0.1'); #R flag removes neck
    topupScript.write('\n');
    topupScript.write('\n');
    topupScript.write('\n');
    topupScript.write('\n');
    topupScript.write('fslmaths ${nodif_brain_mask} -kernel box 10 -dilM ${nodif_brain_mask_dilated}');
    topupScript.write('\n');
    topupScript.write('\n');
    topupScript.close();
    


    # Prepare eddy script
    sys.stdout.write('Compiling run_eddy.sh script\n');
    eddyScript   = open((outPath + '/Step2_topupAndEddyInputs/run_eddy.sh'),'w');
    eddyScript.write('#!/bin/bash \n');
    eddyScript.write('#SBATCH --gpus=1 \n');
    eddyScript.write('#SBATCH --partition=gpu-v100-16g \n');
    eddyScript.write('#SBATCH --cpus-per-task=12 \n');
    eddyScript.write('#SBATCH --time=120 \n');
    eddyScript.write('#SBATCH --error SlurmFiles/eddy_error_%A.txt \n');
    eddyScript.write('#SBATCH --output SlurmFiles/eddy_output_%A.txt \n');
    eddyScript.write('#SBATCH --mem-per-cpu=1G \n \n \n');
    

    eddyScript.write('export OMP_PROC_BIND=TRUE \n');
    #eddyScript.write('export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \n');


#    eddyScript.write('#SBATCH --time=600 \n');
#    eddyScript.write('#SBATCH --gres=gpu:1 \n');
#    eddyScript.write('#SBATCH --constraint=\'pascal|kepler\' \n');
#    eddyScript.write('#SBATCH --mem=8G \n');
    eddyScript.write('\n');
    eddyScript.write('\n');
    eddyScript.write('\nrawPath=' + os.path.dirname(name_base[0]));
    eddyScript.write('\nderPath=' + outPath);
    eddyScript.write('\n');
    eddyScript.write('\n');
    eddyScript.write('imagePath=${derPath}/Step2_topupAndEddyInputs\n');
    eddyScript.write('topupOutputPath=${derPath}/Step3_topupOutput\n');
    eddyScript.write('eddyOutputPath=${derPath}/Step4_eddyOutput\n');
    eddyScript.write('\n');
    eddyScript.write('\n');
    # eddyScript.write('module load OpenBLAS\n');
    # eddyScript.write('module load CUDA/8.0.61\n');
    eddyScript.write('\n');
    eddyScript.write('module load '+fsl_dir+'\n');
    eddyScript.write('source $FSLDIR/etc/fslconf/fsl.sh\n');
    eddyScript.write('\n');
    eddyScript.write('\n');

    if (doTopup==True):
       eddyScript.write('eddy_cuda10.2 \\'); #DEPRECATED eddy_openmp https://neurostars.org/t/eddy-openmp-deprecated-v6-0-6/28189
       # eddyScript.write('eddy_cuda8.0 \\');
       eddyScript.write('\n--verbose \\');
       eddyScript.write('\n--imain=${imagePath}/combined_dwi.nii.gz \\');
       eddyScript.write('\n--mask=${topupOutputPath}/nodif_brain_mask_dilated.nii.gz \\');
       eddyScript.write('\n--acqp=${imagePath}/acqparams.txt \\');
       eddyScript.write('\n--index=${imagePath}/index.txt \\');
       eddyScript.write('\n--bvecs=${imagePath}/combined_dwi.bvec \\');
       eddyScript.write('\n--bvals=${imagePath}/combined_dwi.bval \\');
       eddyScript.write('\n--topup=${topupOutputPath}/topup_results \\');
       eddyScript.write('\n--residuals \\');
       eddyScript.write('\n--data_is_shelled \\');
       eddyScript.write('\n--cnr_maps \\');
       eddyScript.write('\n--repol \\'); # outlier replacement
       # if (slspec.size>0) : eddyScript.write('\n--ol_type=both \\'); # outlier replacement takes into account both slice groups (due to multi band imaging) and individual slices
       # if (slspec.size>0) : eddyScript.write('\n--mporder=6 \\');
       # if (slspec.size>0) : eddyScript.write('\n--slspec=${imagePath}/slspec.txt \\'); # slice-to-volume movement correction
       #eddyScript.write('\n----dont_sep_offs_move \\'); # Do not attempt to separate subject movement from field DC component - DONT USE THIS
       eddyScript.write('\n--dont_peas \\'); # Do not end with an alignment of shells to each other - USE THIS
       eddyScript.write('\n--out=${eddyOutputPath}/eddyResult '); # Do not end with an alignment of shells to each other - USE THIS
       eddyScript.write('\n');
       eddyScript.write('\n');
	    
       averageDirections=True;
       
       for i in range(1,len(name_base)):
           if (blipVolCounts[i]!=blipVolCounts[1]):
               averageDirections=False;

       if os.path.exists(outPath+'/Step0_parsedInputImages/synb0'):
               averageDirections=False;

       if (averageDirections==True):
           eddyScript.write('\n# Spliting volumes for eddy_combine');
           eddyScript.write('\n');
           for i in range(len(name_base)):
              tmp = os.path.basename(name_base[i]);
              eddyScript.write('\nfslroi ${eddyOutputPath}/eddyResult ${eddyOutputPath}/' + tmp + '_eddied ' + str(blipVolCounts[i]) +  ' ' + str(blipVolCounts[i+1]) );
           eddyScript.write('\n');
           eddyScript.write('\n');    
           eddyScript.write('\n# Combine volumes');
           eddyScript.write('\n');
           eddyScript.write( '\neddy_combine \\');
           for i in range(len(name_base)):
              tmp     = os.path.basename(name_base[i]);
              f_bval  = (tmp + '.bval');
              f_bvec  = (tmp + '.bvec');
              svolnum = ('${imagePath}/' + tmp + '_SeriesVolNum.txt');
              eddyScript.write('\n${eddyOutputPath}/' + tmp + '_eddied ${rawPath}/' + f_bval + ' ${rawPath}/' + f_bvec + ' ' + svolnum + ' \\');
           eddyScript.write('\n${eddyOutputPath} 1');
           eddyScript.write('\n');
           eddyScript.write('\n');    
           eddyScript.write('\n# Split rotated directions');
           eddyScript.write('\n');
           eddyScript.write('\nline1=`awk \'NR==1 {print; exit}\' ${eddyOutputPath}/eddyResult.eddy_rotated_bvecs`');
           eddyScript.write('\nline2=`awk \'NR==2 {print; exit}\' ${eddyOutputPath}/eddyResult.eddy_rotated_bvecs`');
           eddyScript.write('\nline3=`awk \'NR==3 {print; exit}\' ${eddyOutputPath}/eddyResult.eddy_rotated_bvecs`');
           eddyScript.write('\n');
           for i in range(len(name_base)):
              tmp = os.path.basename(name_base[i]);
              eddyScript.write('\noutLine1=""');
              eddyScript.write('\noutLine2=""');
              eddyScript.write('\noutLine3=""');
              eddyScript.write('\nfor ((i=' + str(blipVolCounts[i]+1) + '; i<' + str(np.cumsum(blipVolCounts)[i+1]+1 ) + '; i++)); do');
              eddyScript.write('\n   outLine1="$outLine1 `echo $line1 | awk -v N=$i \'{print $N}\'`"');
              eddyScript.write('\n   outLine2="$outLine2 `echo $line2 | awk -v N=$i \'{print $N}\'`"');
              eddyScript.write('\n   outLine3="$outLine3 `echo $line3 | awk -v N=$i \'{print $N}\'`"');
              eddyScript.write('\ndone');
              eddyScript.write('\necho $outLine1 > ${eddyOutputPath}/' + tmp + '_rotated.bvec');
              eddyScript.write('\necho $outLine2 >> ${eddyOutputPath}/' + tmp + '_rotated.bvec');
              eddyScript.write('\necho $outLine3 >> ${eddyOutputPath}/' + tmp + '_rotated.bvec');
              eddyScript.write('\n');
           eddyScript.write('\n');
           eddyScript.write('\n');
           eddyScript.write('\n# Average rotated directions');
           eddyScript.write('\naverage_bvecs=' + average_bvecs_path);
           eddyScript.write('\n');
           eddyScript.write('\n${average_bvecs} \\');
           for i in range(len(name_base)):
              tmp     = os.path.basename(name_base[i]);
              f_bval  = ('${rawPath}/' + tmp + '.bval');
              eddyScript.write('\n' + f_bval + ' ');
              eddyScript.write('${eddyOutputPath}/' + tmp + '_rotated.bvec \\');	
           eddyScript.write('\n${eddyOutputPath}/avg_data');
           eddyScript.write('\n');
           eddyScript.write('\n');
           eddyScript.write('\n# Clean and prep output');
           eddyScript.write('\n');
           eddyScript.write('\nmv ${eddyOutputPath}/data.nii.gz ${eddyOutputPath}/eddied.nii.gz');
           eddyScript.write('\nmv ${eddyOutputPath}/avg_data.bval ${eddyOutputPath}/eddied.bval');
           eddyScript.write('\nmv ${eddyOutputPath}/avg_data.bvec ${eddyOutputPath}/eddied.bvec');
           eddyScript.write('\n');
           eddyScript.write('\n');
           eddyScript.write('\n');
           eddyScript.write('\n# Eddy post processing');
           eddyScript.write('\n');
           eddyScript.write('\nfslmaths ${eddyOutputPath}/eddied.nii.gz -thr 0 ${eddyOutputPath}/eddied.nii.gz');
           eddyScript.write('\n');
           eddyScript.write('\n');
           eddyScript.write('\n');
           eddyScript.write('\nmv ${eddyOutputPath}/eddied.nii.gz ${eddyOutputPath}/prepped_dMRI.nii.gz');
           eddyScript.write('\nmv ${eddyOutputPath}/eddied.bval ${eddyOutputPath}/prepped_dMRI.bval');
           eddyScript.write('\nmv ${eddyOutputPath}/eddied.bvec ${eddyOutputPath}/prepped_dMRI.bvec');
           eddyScript.write('\ncp ${topupOutputPath}/nodif_brain_mask_dilated.nii.gz ${eddyOutputPath}/prepped_dMRI_mask.nii.gz');                   
           eddyScript.write('\n');
           eddyScript.write('\n');
           eddyScript.write('\n');
       else:
           eddyScript.write('\nmv ${eddyOutputPath}/eddyResult.nii.gz ${eddyOutputPath}/prepped_dMRI.nii.gz');
           eddyScript.write('\ncp ${imagePath}/combined_dwi.bval ${eddyOutputPath}/prepped_dMRI.bval');
           eddyScript.write('\ncp ${eddyOutputPath}/eddyResult.eddy_rotated_bvecs ${eddyOutputPath}/prepped_dMRI.bvec');
           eddyScript.write('\ncp ${topupOutputPath}/nodif_brain_mask_dilated.nii.gz ${eddyOutputPath}/prepped_dMRI_mask.nii.gz');
           
    else:
       eddyScript.write('eddy_cuda10.2 \\');
       eddyScript.write('\n--verbose \\');
       eddyScript.write('\n--imain=${imagePath}/combined_dwi.nii.gz \\');
       eddyScript.write('\n--mask=${topupOutputPath}/nodif_brain_mask_dilated.nii.gz \\');
       eddyScript.write('\n--acqp=${imagePath}/acqparams.txt \\');
       eddyScript.write('\n--index=${imagePath}/index.txt \\');
       eddyScript.write('\n--bvecs=${imagePath}/combined_dwi.bvec \\');
       eddyScript.write('\n--bvals=${imagePath}/combined_dwi.bval \\');
       eddyScript.write('\n--residuals \\');
       eddyScript.write('\n--data_is_shelled \\');
       eddyScript.write('\n--cnr_maps \\');
       eddyScript.write('\n--repol \\'); # outlier replacement
       eddyScript.write('\n--dont_peas \\'); # Do not end with an alignment of shells to each other - USE THIS
       eddyScript.write('\n--out=${eddyOutputPath}/eddyResult '); # Do not end with an alignment of shells to each other - USE THIS
       eddyScript.write('\n');
       eddyScript.write('\n');
       eddyScript.write('\n');
       eddyScript.write('\nmv ${eddyOutputPath}/eddyResult.nii.gz ${eddyOutputPath}/prepped_dMRI.nii.gz');
       eddyScript.write('\ncp ${imagePath}/combined_dwi.bval ${eddyOutputPath}/prepped_dMRI.bval');
       eddyScript.write('\ncp ${eddyOutputPath}/eddyResult.eddy_rotated_bvecs ${eddyOutputPath}/prepped_dMRI.bvec');
       eddyScript.write('\ncp ${topupOutputPath}/nodif_brain_mask_dilated.nii.gz ${eddyOutputPath}/prepped_dMRI_mask.nii.gz');

    eddyScript.close();
    # Don't use this flag --dont_sep_offs_move
    # Use this flag --dont_peas
