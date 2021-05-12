import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
import glob
import nibabel as nib
import shutil
plt.style.use('ggplot')
from recon_register import transformMRI, applyTransformMRI, transformMRIRetro, applyTransformMRIRetro
from img_quality_metrics import Comp_Metrics




def FullAnalysis(sub, dicom_dir, nifti_dir, SUBJECTS_DIR, outDir, save, convert, recon_all, transform, metrics) :
    
    '''
    Performs all steps of analysis on prospectively corrected and uncorrected
    scans (if convert, recon_all, transform and metrics are all True)
    
    Parameters
    ----------
    sub : str
        subject ID.
    dicom_dir : str
        path to dicom directory.
    nfiti_dir : str
        path to nifti directory.
    SUBJECTS_DIR : str
        path to the FreeSurfer subject's directory.
    outDir : str
        path to the output directory for registered scans.
    save : str
        string appended to the filename of the output, as e.g. date _05_11
    convert : bool
        whether to convert dicom files to nifti files.
    recon_all : bool
        whether to run the still MPRAGE without MoCo though FreeSurfer's 
        recon_all.
    transform : bool
        whether to register and transform the images.
    metrics : bool
        whether to calculate metrics.  
        

    Returns
    -------
    int
        returns 0, when completed.
    
    '''

    
    '''Convert dicoms into nifti: ''' 
    if convert: 
        list_sequ = os.listdir(dicom_dir)
        sequences = [x for x in list_sequ if x.startswith('TCLMOCO')] 
        
        # check that nifti directory exists, otherwise make a new nifti directory 
        # for this subject: 
        if not os.path.exists(nifti_dir):
            os.makedirs(nifti_dir)
            
        for seq in sequences:
            # sort TRACEW files into B0 and B1000 volumes
            if 'TRACEW_0' in seq:
                files = os.listdir(dicom_dir+seq+'/')
                if len(files)>0:
                    for i in range(0, 27):
                        new = dicom_dir+seq[:-5]+'_B0'+seq[-5:]+'/'+files[i]
                        if not os.path.exists(dicom_dir+seq[:-5]+'_B0'+seq[-5:]):
                            os.makedirs(dicom_dir+seq[:-5]+'_B0'+seq[-5:])
                        shutil.move(dicom_dir+seq+'/'+files[i], new)
                    for i in range(27,54):
                        new = dicom_dir+seq[:-5]+'_B1000'+seq[-5:]+'/'+files[i]
                        if not os.path.exists(dicom_dir+seq[:-5]+'_B1000'+seq[-5:]):
                            os.makedirs(dicom_dir+seq[:-5]+'_B1000'+seq[-5:])
                        shutil.move(dicom_dir+seq+'/'+files[i], new)
        
        
        # check that Dicom files in directory are complete
        error=0
        list_sequ = os.listdir(dicom_dir)
        sequences = [x for x in list_sequ if x.startswith('TCLMOCO')]
        for seq in sequences:
            nr_files = np.array ([192, 25, 21, 160, 28, 27, 54, 27, 27])
            tag = np.array(['T1_MPR_', 'T2_TSE_', 'T1_TIRM_', 'T2_FLAIR_', 'T2STAR_', 'ADC', 'TRACEW_0', 'TRACEW_B0', 'TRACEW_B1000'])
            for i in range(len(tag)):
                check = seq.find(tag[i])
                if check > 0:
                    length = len(os.listdir(dicom_dir + seq))   
                    if length != nr_files[i]:
                        if tag[i] == 'TRACEW_0':
                            if length != 0:
                                print('ERROR: number of files for '+seq+' is not correct. Found '+str(length)+', but expected '+str(nr_files[i])+'or 0')
                                error=1
                        else:
                            print('ERROR: number of files for '+seq+' is not correct. Found '+str(length)+', but expected '+str(nr_files[i]))
                            error=1
        
        
        if error == 0:
            print('All DICOM files complete.')
            print('#######')
            print('#######')
                            
        
            for seq in sequences:
                if len(os.listdir(dicom_dir+seq))>0:
                    dcm = os.listdir(dicom_dir + seq)[0]
                    in_volume = dicom_dir + seq + '/' + dcm
                    out_volume = nifti_dir + 'TCL'+sub+'_' + seq + '.nii'
                    if seq.startswith('TCLMOCO_OFF_STILL_T1_MPR_3D_SAG'):
                        nifti_path = out_volume         
                              
                    for tag in ['T1_MPR_', 'T2_TSE_', 'T1_TIRM_', 'T2_FLAIR_', 'T2STAR_', 'SWI']:
                        if tag in seq:
                            subprocess.run('mri_convert ' + in_volume + ' ' + out_volume, shell=True)
                            
                    
                    for tag in ['TRACEW_B0', 'TRACEW_B1000', 'ADC']:
                        if tag in seq:
                            subprocess.run('mri_convert ' + in_volume + ' ' + out_volume+' --no-dwi', shell=True)
                            
        
        else:
            print('check number of files and try again')
    
    print('//////////////////////////////////////////////////////')
    print('//                                                  //')
    print('//                Conversion done                   //')
    print('//                                                  //')
    print('//////////////////////////////////////////////////////')
    
    
    '''Recon_all:'''
    # Run recon all:
    if recon_all:
        list_sequ = os.listdir(dicom_dir)
        sequences = [x for x in list_sequ if x.startswith('TCLMOCO')]
        for seq in sequences:
            if len(os.listdir(dicom_dir+seq))>0:                
                out_volume = nifti_dir + 'TCL'+sub+'_' + seq + '.nii'
                if seq.startswith('TCLMOCO_OFF_STILL_T1_MPR_3D_SAG'):
                    nifti_path = out_volume 
        
        
        subj_id = sub
        print(subj_id, nifti_path)
        subprocess.run('recon-all -i ' + nifti_path + ' -s ' + subj_id + ' -sd '+SUBJECTS_DIR+' -all -parallel', 
                       shell=True)
        print('//////////////////////////////////////////////////////')
        print('//                                                  //')
        print('//                  ReconAll done                   //')
        print('//                                                  //')
        print('//////////////////////////////////////////////////////')
    
    
    
    '''Registration:''' 
    
    if transform:
        # get all transforms:
        transformMRI(sub, nifti_dir, outDir)
        print('//////////////////////////////////////////////////////')
        print('//                                                  //')
        print('//                TransformMRI done                 //')
        print('//                                                  //')
        print('//////////////////////////////////////////////////////')
        
        
        # brainMask all images:
        if os.path.exists(SUBJECTS_DIR + sub + '/mri/brainmask_edit.mgz'):
            brainmask = SUBJECTS_DIR + sub + '/mri/brainmask_edit.mgz'
        else:
            brainmask = SUBJECTS_DIR + sub + '/mri/brainmask.mgz'
        applyTransformMRI(nifti_dir, brainmask, outDir)
        print('//////////////////////////////////////////////////////')
        print('//                                                  //')
        print('//             applyTransformMRI done               //')
        print('//                                                  //')
        print('//////////////////////////////////////////////////////')
        
        
        # check brainmasks correctly registered:
        plt.figure(figsize=(10,7))
        i = 1
        for tag in ['T1_MPR_', 'T2_FLAIR_', 'T2_TSE_', 'T1_TIRM_', 'T2_STAR_', 'ADC', 'TRACEW_B0', 'TRACEW_B1000']:   #'EPI_SWI',
            if len(glob.glob(outDir+'/*TCLMOCO_OFF_STILL_*'+tag+'*moved.nii')) > 0:
                img_file = glob.glob(outDir+'/*TCLMOCO_OFF_STILL_*'+tag+'*moved.nii')[0]
                img = nib.load(img_file).get_fdata().astype(np.uint16)
                if tag == 'T1_MPR_':
                    bm_file = outDir+'brainmask_bin.nii'
                else:
                    bm_file = glob.glob(outDir+'bm*TCLMOCO_OFF_STILL_*'+tag+'*.nii')[0]
                bm = nib.load(bm_file).get_fdata().astype(np.uint16)
                
                sh = np.shape(img)
                sl_dir = np.argmin(sh)
                
                slice_dim = sl_dir
                indx = [slice(None)]*img.ndim
                indx[slice_dim] = int(sh[sl_dir]/2)
            
                plt.subplot(2,4,i)
                plt.imshow(img[indx], cmap='gray')   
                plt.imshow(bm[indx], cmap='Reds', alpha=0.4)
                plt.axis('off')
                plt.title('Brainmask '+tag, fontsize = 10)
                i+=1
        
        plt.savefig(outDir+'Check_Brainmask', dpi=300)
        plt.show()
    
    
    '''Calculate Metrics:'''   
    if metrics:
        folder = '/data1/hannah/Registrations/'+sub+'/'
        all_img = glob.glob(folder +'*_moved.nii')
        out_dir = '/data1/hannah/Metrics_Results/'+sub+'/'
        
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        # sort different sequences/contrasts:
        test = 0
        if len(all_img) == 0:
            print('No moved images in directory. Check subdirectories.')
        for i in range(len(all_img)):
            tmp, filename = os.path.split(all_img[i])
            # search for different sequence types:
            for tag in ['T1_MPR_', 'T2_TSE_', 'T1_TIRM_', 'T2_FLAIR_', 'T2STAR_', 'EPI_SWI', 'ADC', 'TRACEW_B0', 'TRACEW_B1000']:
                test = filename.find(tag)
                if test > 0:
                    newpath = tmp + '/' + tag
                    if not os.path.exists(newpath):
                        os.makedirs(newpath)
                    shutil.move(all_img[i], tmp + '/' + tag + '/' + filename)
                    test = 1
            if test == 0:
                print('ERROR: Filename ' + filename + ' does not contain the right sequence type description')
        
        
        for tag in ['T1_MPR_', 'T2_TSE_', 'T1_TIRM_', 'T2_FLAIR_', 'T2STAR_', 'EPI_SWI', 'ADC', 'TRACEW_B0', 'TRACEW_B1000']:
            if os.path.isdir(folder + tag):
                print('######## ' + tag)
                img_mov = glob.glob(folder + tag + '/' +'*_moved.nii')
                check = 0
                for file in img_mov:
                    test = file.find('TCLMOCO_OFF_STILL_')
                    if test > 0:
                        ref = nib.load(file).get_fdata().astype(np.uint16)
                        check = 1
                if check == 0:
                    print('ERROR: no reference Scan in ' + tag)
        
                if tag == 'T1_MPR_':
                    brainmask_bin = nib.load('/data1/hannah/Registrations/'+sub+'/brainmask_bin.nii').get_fdata().astype(np.uint16)
                else:
                    bm_mov = glob.glob('/data1/hannah/Registrations/'+sub+'/bm_mov_*'+tag+'*.nii')[0]
                    brainmask_bin = nib.load(bm_mov).get_fdata().astype(np.uint16)
        
                # use function Comp_Metrics to calculate the metrics:
                psnr, ssim, tg, names = Comp_Metrics(img_mov, ref,   
                                                     brainmask=brainmask_bin,
                                                     whole_image=False, compute_gr=True,
                                                     normal=True)
        
                print(names)
        
                save_arr = np.array([names, ssim, psnr, tg]).T
                np.savetxt(out_dir+'Values_'+save+tag, save_arr, fmt='%s', header='Values of the metrics SSIM, PSNR and TG for acquisitions specified in first column')
    
    return 0
       

def RetroAnalysis(sub, orig_dir, corr_dir, nifti_dir, outDir, copy, transform, calculate):
    '''
    Performs all steps of analysis for MPRAGE and FLAIR scans 
    reconstructed retrospectively scans (if convert, recon_all, 
    transform and metrics are all True)

    Parameters
    ----------
    sub : str
        subject ID.
    orig_dir : str
        path to directory with original niftie files 
        from retrospective reconstruction.
    nifti_dir : str
        path to direcotry with renamed niftie files.
    corr_dir : str
        path to directory with bias field corrected nfiti files.
    outDir : str
        path to the output directory for registered scans.
    copy : bool
        whether to copy nifti files with different filename from
        orig_dir to nifti_dir.
    transform : bool
        whether to register and transform the images.
    calculate : bool
        whether to calculate metrics.

    Returns
    -------
    int
        returns 0, when completed.

    '''
    
    
    ''' Copy data and change description to Capital letters: '''
    if copy:
        for file in glob.glob(orig_dir+'TCLmoco*'):
            filename = os.path.basename(file)

            # filename in capital letters and new name:
            filename = filename[:-4].upper()+'.nii'
            if 'T1_MPR' in filename:
                ind = filename.find('ISO')
                newname = filename[:ind+3]
            elif 'FLAIR' in filename:
                ind = filename.find('CHANGED')
                newname = filename[:ind+7]


            # sort after retro:
            if 'TOS_' in filename:
                ind = newname.find('OFF')
                newname = newname[:ind]+'RETRO'+newname[ind+3:]

            if 'REACQ' not in filename and 'STILL' not in filename:
                newname = newname+'_RR'

            if not os.path.exists(nifti_dir):
                print('New nifti_dir created.')
                os.makedirs(nifti_dir)

            shutil.copy(file, nifti_dir+newname+'.nii')


    '''Registration:'''
    if transform:
        # get all transforms:
        transformMRIRetro(sub, corr_dir, outDir, corr=True)    
        print('//////////////////////////////////////////////////////')
        print('//                                                  //')
        print('//                TransformMRI done                 //')
        print('//                                                  //')
        print('//////////////////////////////////////////////////////')


        # brainMask all images:
        brainmask = SUBJECTS_DIR + sub + '/mri/brainmask.mgz'
        applyTransformMRIRetro(corr_dir, brainmask, outDir, corr=True)   
        print('//////////////////////////////////////////////////////')
        print('//                                                  //')
        print('//             applyTransformMRI done               //')
        print('//                                                  //')
        print('//////////////////////////////////////////////////////')


        # check brainmasks correctly registered:
        plt.figure(figsize=(10,10))
        i = 1
        for tag in ['T1_MPR_', 'T2_FLAIR']:
            if len(glob.glob(outDir+'retro_scans/*corr_TCLMOCO_OFF_STILL_'+tag+'*_moved.nii')) > 0:
                img_file = glob.glob(outDir+'retro_scans/*corr_TCLMOCO_OFF_STILL_'+tag+'*_moved.nii')[0] 
                img = nib.load(img_file).get_fdata().astype(np.uint16)
                if tag == 'T1_MPR_':
                    bm_file = outDir+'brainmask_bin.nii'
                    bm = nib.load(bm_file).get_fdata().astype(np.uint16)
                else:
                    bm_file = glob.glob(outDir+'bm_mov_retro*TCLMOCO_OFF_STILL_'+tag+'*.npy')[0]
                    bm = np.load(bm_file).astype(np.uint16)

                sh = np.shape(img)
                sl_dir = np.argmin(sh)

                slice_dim = sl_dir
                indx = [slice(None)]*img.ndim
                indx[slice_dim] = int(sh[sl_dir]/2)

                plt.subplot(1,2,i)
                plt.imshow(img[indx], cmap='gray')
                plt.imshow(bm[indx], cmap='Reds', alpha=0.4)
                plt.axis('off')
                plt.title('Brainmask '+tag)
                i+=1
        plt.savefig(outDir+'Check_Brainmask_retro', dpi=300)
        plt.show()


    '''Calculate Metrics:'''
    if calculate:
        folder = '/data1/hannah/Registrations/'+sub+'/retro_scans/'
        all_img = glob.glob(folder +'*_moved.nii')
        out_dir = '/data1/hannah/Metrics_Results/'+sub+'/'

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # sort different sequences/contrasts:
        test = 0
        if len(all_img) == 0:
            print('No moved images in directory. Check subdirectories.')
        for i in range(len(all_img)):
            tmp, filename = os.path.split(all_img[i])
            # search for different sequence types:
            for tag in ['T1_MPR_', 'T2_FLAIR_']:
                test = filename.find(tag)
                if test > 0:
                    newpath = tmp + '/' + tag
                    if not os.path.exists(newpath):
                        os.makedirs(newpath)
                    shutil.move(all_img[i], tmp + '/' + tag + '/' + filename)
                    test = 1
            if test == 0:
                print('ERROR: Filename ' + filename + ' does not contain the right sequence type description')


        for tag in ['T1_MPR_', 'T2_FLAIR_']:
            if os.path.isdir(folder + tag):
                print('######## ' + tag + ' RetroMoCo')
                img_mov = glob.glob(folder  + tag + '/' +'*_corr_*_moved.nii')
                check = 0
                for file in img_mov:
                    test = file.find('TCLMOCO_OFF_STILL_')
                    if test > 0:
                        ref = nib.load(file).get_fdata().astype(np.uint16)
                        check = 1
                if check == 0:
                    print('ERROR: no reference Scan in ' + tag)

                if tag == 'T1_MPR_':
                    brainmask_bin = nib.load('/data1/hannah/Registrations/'+sub+'/brainmask_bin.nii').get_fdata().astype(np.uint16)
                else:
                    bm_mov = glob.glob('/data1/hannah/Registrations/'+sub+'/bm_mov_retro*'+tag[:-1]+'*.npy')[0]
                    brainmask_bin = np.load(bm_mov).astype(np.uint16) 

                # use function Comp_Metrics to calculate the metrics:
                psnr, ssim, tg, names = Comp_Metrics(img_mov, ref,
                                                     brainmask=brainmask_bin,
                                                     whole_image=False, compute_gr=True,
                                                     normal=True) 

                print(names)

                save_arr = np.array([names, ssim, psnr, tg]).T  
                np.savetxt(out_dir+'Values_Retro_'+save+tag, save_arr, fmt='%s', header='Values of the metrics SSIM, PSNR and TG for acquisitions specified in first column')    
    
    return 0


def PerformGCut(sub, SUBJECTS_DIR):
    '''
    Performs graph cut brainmask editing.

    Parameters
    ----------
    sub : str
        subject ID.
    SUBJECTS_DIR : str
        path to FreeSurfer's subjects directory.

    Returns
    -------
    int
        returns 0, when completed.

    '''
    
    shutil.copy(SUBJECTS_DIR+sub+'/mri/brainmask.mgz', SUBJECTS_DIR+sub+'/mri/brainmask_tmp.mgz')
    
    
    subprocess.run('recon-all -skullstrip -clean-bm -gcut -subjid ' + sub, 
                   shell=True)
    
    os.rename(SUBJECTS_DIR+sub+'/mri/brainmask.mgz', SUBJECTS_DIR+sub+'/mri/brainmask_edit.mgz')
    os.rename(SUBJECTS_DIR+sub+'/mri/brainmask_tmp.mgz', SUBJECTS_DIR+sub+'/mri/brainmask.mgz')
    
    print('fcut done - PLEASE CONTROL OUTPUT: load brainmask.gcuts.mgz overlayed on T1.mgz in Freeview')
    return 0



''' (1) Run analysis on prospectively corrected and uncorrected data:'''
# define specific input parameters for the current run:
subjs = ['Protoype_01']
for sub in subjs:
    dicom_dir = '/mnt/mocodata1/Data_Analysis/DICOMS/'+sub+'/'   
    nifti_dir = '/data1/hannah/NIFTIS/'+sub+'/'
    SUBJECTS_DIR = '/mnt/mocodata1/Data_Analysis/Data_Recon_All/'
    outDir = '/data1/hannah/Registrations/'+sub+'/'
    
    save = '05_11_'
    
    # which steps to perform:
    convert = False
    recon_all = False
    transform = False
    metrics = False
    
    FullAnalysis(sub, dicom_dir, nifti_dir, SUBJECTS_DIR, outDir, save, convert, recon_all, transform, metrics)


    # Perform gcuts:
    #PerformGCut(sub, SUBJECTS_DIR)



''' (2) Run analysis on retrospectively corrected data'''
for sub in subjs:
    SUBJECTS_DIR = '/mnt/mocodata1/Data_Analysis/Data_Recon_All/'
    outDir = '/data1/hannah/Registrations/'+sub+'/'
    nifti_dir = '/mnt/mocodata1/Data_Analysis/NIFTIS_Retro/'+sub
    corr_dir = '/mnt/mocodata1/Data_Analysis/NIFTIS_Retro_Corr/'+sub
    orig_dir = '/mnt/mocodata1/MoCoHealthy/RetroMoCoData/'+sub

    save = '05_11_'
    
    # which steps to perform
    copy = False
    transform = False
    calculate = False
    
    RetroAnalysis(sub, orig_dir, corr_dir, nifti_dir, outDir, copy, transform, calculate)

