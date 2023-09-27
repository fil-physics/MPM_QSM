%%% Description: MPM SWI pipeline
% main steps:
% 1) rotation of the magnitude and phase into scanner space
%   (in complex domain so that phase wraps are not blurred)
% 2) CLEAR SWI 
%   - SNR-optimal magnitude combination,
%   - phase mask calculation and multiplication with combined magnitude
%   - minimum intensity projection

%%% Publications:
% Please remember to give credit to the authors of the methods used:
% 1. CLEARSWI paper: Eckstein, Korbinian, et al. Neuroimage 237 (2021): 118175.


%%% Inputs:
% clearswi_command       : path to  followed by clearswi command, i.e. (in linux) '/your_path/bin/clearswi' or (in windows) 'D:\your_path\bin\clearswi'
% in_root_dir            : root directory to input nifti files
% out_root_dir           : root directory to output nifti files
% sensitivity_corr       : corrects for bias field, 'on' or 'off'
% filter_size            : high-pass filter kernel size, for 7T MPM [2 2 0] is adviced

%%%% Inputs - directories, parameters and files specific to given contrast
% ATTENTION: ensure only niftis you want to use are in that folder, with increasing echo numbering:
% mag_dir                : folder with magnitude niftis
% ph_dir                 : folder with phase inftis
% TEs                    : echo time in ms
% output_dir             : output SWI directory for a specific MPM contrast

%%% Outputs:
%%%% combined final results in output_dir:
% clearswi.nii          : SWI image
% mip.nii               : Minimum Intensity Projection over 7 axial slices

%%%% additional outputs:
% mag_rot.nii               : phase rotated into scanner space (if FOV was angulated)
% ph_rot.nii                : magnitude rotated into scanner space (if FOV was angulated)

% script created by Barbara Dymerska
% @ UCL FIL Physics

function MPM_SWI(para)

mag_fulldir = fullfile(para.in_root_dir, para.mag_dir) ;
ph_fulldir = fullfile(para.in_root_dir, para.ph_dir) ;

output_fulldir = fullfile(para.out_root_dir, para.output_dir) ;

if ~exist(output_fulldir, 'dir')
    mkdir(output_fulldir)
end
cd(output_fulldir)

TEs = para.TEs ;
if isempty(para.ph_file) && isempty(para.mag_file)
    ph_files = spm_select('FPList', ph_fulldir, '^.*\.(nii|img)$');
    mag_files = spm_select('FPList', mag_fulldir, '^.*\.(nii|img)$');

    ph_1tp = nifti(ph_files(1,:));
    ph = zeros([size(ph_1tp.dat) size(TEs,2)]);
    mag = zeros([size(ph_1tp.dat) size(TEs,2)]);
    for t = 1:size(TEs,2)

        ph_1tp = nifti(ph_files(t,:));
        ph(:,:,:,t) = ph_1tp.dat(:,:,:) ;

        mag_1tp = nifti(mag_files(t,:));
        mag(:,:,:,t) = mag_1tp.dat(:,:,:) ;

    end
    ph(~isfinite(ph))=0;
    mag(~isfinite(mag))=0;
    mat_image = mag_1tp.mat ;
    
else
    ph_merged_file = fullfile(ph_fulldir, para.ph_file) ;
    mag_merged_file = fullfile(mag_fulldir, para.mag_file) ;
    mag = nifti(mag_merged_file).dat(:,:,:,:);
    ph = nifti(ph_merged_file).dat(:,:,:,:);
    mat_image = nifti(mag_merged_file).mat;
end
    data_dim = size(mag) ;
    cmplx=mag.*exp(1i*ph);

    real_merged_file = fullfile(output_fulldir, 'real.nii') ;
    imag_merged_file = fullfile(output_fulldir, 'imag.nii') ;

    createNifti(real(cmplx), real_merged_file, mat_image)
    createNifti(imag(cmplx), imag_merged_file, mat_image)
clear cmplx mag ph
%% rotation of real and imaginary images
% defining affine matrix in scanner space for data rotation to scanner
% space with mantaining the same image origin (i.e. no translation)
disp('data rotation to scanner space')
real_V = spm_vol(real_merged_file) ;
imag_V = spm_vol(imag_merged_file) ;

Z = spm_imatrix(mat_image) ;
pixdim = Z(7:9);
O = mat_image\[0 0 0 1]' ;
O = O(1:3)' ;

mat_scanner(1,:) = [pixdim(1) 0 0 -pixdim(1)*O(1)] ;
mat_scanner(2,:) = [0 pixdim(2) 0 -pixdim(2)*O(2)] ;
mat_scanner(3,:) = [0 0 pixdim(3) -pixdim(3)*O(3)] ;
mat_scanner(4,:) = [0 0 0 1] ;


img2scanner_mat = mat_image\mat_scanner ;
real_rot = zeros(data_dim) ;
imag_rot = zeros(data_dim) ;
data_dim_xy = data_dim(1:2);
for echo = 1:data_dim(4)
for slice = 1 : data_dim(3)
    real_rot(:,:,slice,echo) = spm_slice_vol(real_V(echo), img2scanner_mat*spm_matrix([0 0 slice]), data_dim_xy, -7) ;
    imag_rot(:,:,slice,echo) = spm_slice_vol(imag_V(echo), img2scanner_mat*spm_matrix([0 0 slice]), data_dim_xy, -7) ;
end
end
real_rot(~isfinite(real_rot)) = 0 ;
imag_rot(~isfinite(imag_rot)) = 0 ;
cmplx_rot = real_rot + 1i*imag_rot;

% saving rotated magnitude and phase
mag_rot_merged_file = fullfile(output_fulldir, 'mag_rot.nii') ;
ph_rot_merged_file = fullfile(output_fulldir, 'ph_rot.nii') ;
createNifti(abs(cmplx_rot), mag_rot_merged_file, mat_scanner)
createNifti(angle(cmplx_rot), ph_rot_merged_file, mat_scanner)
clear cmplx_rot real_rot imag_rot

%%
disp('CLEAR SWI calculation')
status = system(sprintf('%s -p %s -m %s -t [%s] --filter-size [%s] --mag-sensitivity-correction %s', para.clearswi_command, ph_rot_merged_file, mag_rot_merged_file, num2str(TEs), num2str(para.filter_size), para.sensitivity_corr)) ;

if status == 1
    error('CLEARSWI did not run properly - check your installation path')
end

if para.data_cleanup
    delete(real_merged_file,imag_merged_file,mag_rot_merged_file,ph_rot_merged_file)
end