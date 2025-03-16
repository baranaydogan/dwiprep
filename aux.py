#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:31:50 2018

@author: baran
"""

# source and update environment
def source(prepEnv, script):
    """
    http://pythonwise.blogspot.fr/2010/04/sourcing-shell-script.html (Miki Tebeka)
    http://stackoverflow.com/questions/3503719/#comment28061110_3505826 (ahal)
    """
    import subprocess
    
    proc = subprocess.Popen(
        ['bash', '-c', 'set -a && source {} && env -0'.format(script)], 
        stdout=subprocess.PIPE, 
        shell=False, 
        env=prepEnv);
    
    output, err = proc.communicate();
    output = output.decode('utf8');
    env = dict((line.split("=", 1) for line in output.split('\x00') if line));
    return env;


# A,B,out are file paths for nifti images
# writes out=A-B
def nifti_subtract(A,B,out):
    
    import nibabel as nib
    import numpy as np
    
    A       = nib.load(A);
    B       = nib.load(B);
    
    img_A   = np.squeeze(A.get_data());
    img_B   = np.squeeze(B.get_data());
    
    diff    = img_A - img_B;
    
    if (A.header.sizeof_hdr==348):
        C   = nib.Nifti1Image(diff,A.affine);
    else:
        C   = nib.Nifti2Image(diff,A.affine);
    
    nib.save(C,out);


# Write new nifti image like the provided nifti template
def nifti_write_like(data,nii_template,file_name):

    import nibabel as nib
    
    if (nii_template.header.sizeof_hdr==348):
        nio = nib.Nifti1Image(data,nii_template.affine);
    else:
        nio = nib.Nifti2Image(data,nii_template.affine);
    
    nib.save(nio,file_name);
    
# Extract b0s
def nifti_extractB0s(out,basenameInp,basenameOut,b0_threshold,report):
    
    import sys
    import nibabel as nib
    import numpy as np
    from glob import glob
    
    sys.stdout.write('   Extracting B0 images... ');
    
    nii  = nib.load(out);
    bval = np.loadtxt((basenameInp + '.bval'));
    
    img  = nii.get_data();
    
    img_b0 = img[:,:,:,bval<b0_threshold];
    ave_b0 = np.mean(np.mean(np.mean(img_b0,axis=0),axis=0),axis=0);
    
#    if (img_b0.shape[3]==0):
#        report.write('\n}');
#        report.write(',\n"result": "FAILED, no b0 found"');
#        report.write('\n}');
#        report.close();
#        return [];

    if (len(np.shape(img_b0))==5):
        if (img_b0.shape[4]==1):
            img_b0=np.squeeze(img_b0,4)
        else:
            report.write('\n}');
            report.write(',\n"result": "FAILED, b0 can not be 5 dimensionsal"');
            report.write('\n}');
            report.close();
            return [];

    if (nii.header.sizeof_hdr==348):
        nio   = nib.Nifti1Image(img_b0,nii.affine);
    else:
        nio   = nib.Nifti2Image(img_b0,nii.affine);
    
    out = (basenameOut + '_b0.nii.gz');
    
    nib.save(nio,out);
    
    if not glob(out):
        report.write('\n}');
        report.write(',\n"result": "FAILED, extracting b0 images failed"');
        report.write('\n}');
        report.close();
        sys.exit('extracting b0 images failed');
        
    sys.stdout.write('DONE \n');
    
    report.write(',\n"individual means of b0 images":' + str(np.array2string(ave_b0,                 separator=', ')))
    report.write(',\n"mean of all b0 images":'         + str(np.array2string(np.mean(ave_b0),        separator=', ')))
    report.write(',\n"number of b0 images":'           + str(np.array2string(sum(np.atleast_1d(bval<b0_threshold)), separator=', ')))
    report.write(',\n"output": "'+ out +'"');
    
    return out;


# Extract a b0
def nifti_extract_a_b0(out,basenameInp,basenameOut,b0_threshold):
    
    import sys
    import nibabel as nib
    import numpy as np
    from glob import glob
    
    sys.stdout.write('   Extracting a B0 image... ');
    
    nii  = nib.load(out);
    bval = np.loadtxt((basenameInp + '.bval'));
    
    img  = nii.get_data();
    
    img_b0 = img[:,:,:,np.argmax(bval<b0_threshold)];
    
    if (nii.header.sizeof_hdr==348):
        nio   = nib.Nifti1Image(img_b0,nii.affine);
    else:
        nio   = nib.Nifti2Image(img_b0,nii.affine);
    
    out = (basenameOut + '.nii.gz');
    
    nib.save(nio,out);

    np.savetxt((basenameOut+'.bvec'), np.zeros((3,1),'float'), delimiter=' ',  fmt='%1.6f', newline='\n');
    np.savetxt((basenameOut+'.bval'), np.zeros((1,1),'int'),   delimiter=' ',  fmt='%d',    newline=' ');
    
    if not glob(out):
        sys.exit('extracting b0 images failed');
        
    sys.stdout.write('DONE \n');
    
    return out;


# Merge all images
def nifti_concatenateImages(out,imageList):
    
    import sys
    import nibabel as nib
    import numpy as np

    nii  = nib.load(imageList[0]);
    img  = nii.get_data();
    
#    if (img.ndim<4):
#        img = np.expand_dims(img,axis=4);
    
    nv   = np.zeros((len(imageList)),dtype=int);
    nv[0]= img.shape[3];
    
    for i in range(1,len(imageList)):
        tmp     = nib.load(imageList[i]).get_data();
        nv[i]   = tmp.shape[3];
        img     = np.concatenate((img,tmp),axis=3);
        
    if (nii.header.sizeof_hdr==348):
        nio   = nib.Nifti1Image(img,nii.affine);
    else:
        nio   = nib.Nifti2Image(img,nii.affine);
    
    nib.save(nio,out);
    
    return nv;

# Extract given volumes from input image
def nifti_extractVolumes(out,inp,indices):
    
    import nibabel as nib    

    nii  = nib.load(inp);
    img  = nii.get_data();
    
    img  = img[:,:,:,indices];
        
    if (nii.header.sizeof_hdr==348):
        nio   = nib.Nifti1Image(img,nii.affine);
    else:
        nio   = nib.Nifti2Image(img,nii.affine);
    
    nib.save(nio,out);

# Merge all dwi and their bvals and bvecs
def nifti_concatenateDwi(out,dwiList):

    if (len(dwiList)<1):
        return ([],[]);
    
    import os
    import nibabel as nib
    import numpy as np
        
    # Merge dwi data
    nii  = nib.load(dwiList[0]);
    img  = nii.get_data();
    if (len(img.shape)==3):
        img = np.expand_dims(img,3);
    
    for i in range(1,len(dwiList)):
        tmp     = nib.load(dwiList[i]).get_data();
        if (len(tmp.shape)==3):
            tmp = np.expand_dims(tmp,3);
        img     = np.concatenate((img,tmp),axis=3);
        
    if (nii.header.sizeof_hdr==348):
        nio   = nib.Nifti1Image(img,nii.affine);
    else:
        nio   = nib.Nifti2Image(img,nii.affine);

    
    # Merge bvals and bvecs
    merged_bvals = np.empty((0,),float);
    merged_bvecs = np.empty((3,0),float);

    for i in range(len(dwiList)):
        
        dwi = dwiList[i];
        
        if (dwi[-3:]=='.gz'):
            base_inp    = os.path.basename(dwi[0:-7]);
        else:
            base_inp    = os.path.basename(dwi[0:-4]);

        basename        = os.path.dirname(dwi) + '/' + base_inp;
        
        cur_bvals       = np.loadtxt((basename + '.bval'));
        cur_bvecs       = np.loadtxt((basename + '.bvec'));
        
        if (cur_bvals.ndim==0):
            tmp         = np.empty((1,),float);
            tmp[0]      = cur_bvals;
            cur_bvals   = tmp;
            
            tmp         = np.empty((3,1),float);
            tmp[:,0]    = cur_bvecs;
            cur_bvecs   = tmp;

        merged_bvals    = np.concatenate((merged_bvals,cur_bvals),axis=0);
        merged_bvecs    = np.concatenate((merged_bvecs,cur_bvecs),axis=1);

        
    # Write output
    if (out[-3:]=='.gz'):
        base_out = os.path.basename(out[0:-7]);
    else:
        base_out = os.path.basename(out[0:-4]);
    
    nib.save(nio,out);
    
    np.savetxt((os.path.dirname(out) + '/' + base_out+'.bval'), merged_bvals.astype(int),   delimiter='',   fmt='%d',    newline=' ' );
    np.savetxt((os.path.dirname(out) + '/' + base_out+'.bvec'), merged_bvecs.astype(float), delimiter=' ',  fmt='%1.6f', newline='\n');
    
        
    return (out,merged_bvals.shape[0]);
    
    
    
    
    
    
    
