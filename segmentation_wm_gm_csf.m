% Freesurfer labels: https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT

PATH   = 'C:\Users\yangq\OneDrive - Queensland University of Technology\MGH Connectome';
SUBJECT    = 'sub_001';
datatype = '_real'; % '_real' or empty string ''
FOLDER = [PATH '\' SUBJECT '\dwi' datatype '\'];
seg = double(niftiread([FOLDER SUBJECT '_dwi' datatype '_aparc+aseg.nii.gz'])); % load the segmentations

l = 71; % select an axial slice

% white matter
mask_wm_L = seg==2; % left cerebral white matter
mask_wm_R = seg==41; % right cerebral white matter
mask_wm = mask_wm_L + mask_wm_R; % white matter mask
figure, imagesc(mask_wm(:,:,l)),axis image, view(-90,90), axis off, title('white matter')

% cortical gray matter
mask_gm_L = seg>1000 & seg<2000; % left cerebral cortex
mask_gm_R = seg>2000 & seg<3000; % right cerebral cortex
mask_gm = mask_gm_L + mask_gm_R; % gray matter mask
figure, imagesc(mask_gm(:,:,l)),axis image, view(-90,90), axis off, title('gray matter')

% CSF
mask_csf_L = seg==4; % Left-Lateral-Ventricle
mask_csf_R = seg==43; % Right-Lateral-Ventricle
mask_csf = mask_csf_L + mask_csf_R; % CSF mask
figure, imagesc(mask_csf(:,:,l)),axis image, view(-90,90), axis off, title('CSF')