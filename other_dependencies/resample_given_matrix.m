function resample_given_matrix

root_dir = '/media/barbara/hdd2/DATA/FIL/7T/M700153_analysis/MORSE_Opt_phase/SEPIA/pdw_RR_54_55_rot/';

% Original image
im = nifti([root_dir 'mask.nii']);

% Copy of original
im1 = im ;

% 30 degree rotation about x to be applied
R = spm_matrix([0,0,0,deg2rad(30)]);

% Combined mat field, 30deg rotation + original
M = inv(R)*spm_get_space(im.dat.fname);

% Update mat field of the copy to this new oritentation
spm_get_space(im1.dat.fname, M);

% Reslice it relative to the original, will add prefix 'r'
spm_reslice({im.dat.fname, im1.dat.fname});

disp(['File ' root_dir 'rim1.nii will be in correct voxel space (x,z,y)'])
