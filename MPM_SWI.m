%%% Description: MPM SWI pipeline
% main steps:
% 1) rotation of the magnitude and phase into scanner space
%   (in complex domain so that phase wraps are not blurred)
% 2) CLEAR SWI
%   - SNR-optimal magnitude combination,
%   - phase mask calculation and multiplication with combined magnitude
%   - minimum intensity projection
% 3) rotation  of CLEAR-SWI and MIP back into image space

%%% Publications:
% Please remember to give credit to the authors of the methods used:
% 1. CLEARSWI paper: Eckstein, Korbinian, et al. Neuroimage 237 (2021): 118175.


%%% Inputs in para. - structure:
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
%%% in the scanner space:
% clearswi.nii          : SWI image
% mip.nii               : Minimum Intensity Projection over 7 axial slices
%%% in the original image space:
% clearswi_invrot.nii   : SWI image
% mip_invrot.nii        : Minimum Intensity Projection over 7 axial slices

%%%% optional extra outputs:
% mag_rot.nii               : magnitue rotated into scanner space (if FOV was angulated)
% ph_rot.nii                : phase rotated into scanner space (if FOV was angulated)

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
[ph_min,ph_max]  = bounds( ph(:) , "all" );
ph = 2*pi*(ph - ph_min)/(ph_max-ph_min) ;

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

mat_scanner_read_x(1,:) = [0 0 pixdim(3) -pixdim(3)*O(3)] ;
mat_scanner_read_x(2,:) = [pixdim(1) 0 0 -pixdim(1)*O(1)] ;
mat_scanner_read_x(3,:) = [0 pixdim(2) 0 -pixdim(2)*O(2)] ;
mat_scanner_read_x(4,:) = [0 0 0 1] ;

img2scanner_mat = mat_image\mat_scanner_read_x ;
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

% readout direction moved to z, i.e. 3rd dimension
% to get minimal intensity projection for axial slices using clearswi:
cmplx_rot = permute(cmplx_rot, [3 1 2 4]);
mat_scanner_read_z(1,:) = [pixdim(1) 0 0 -pixdim(1)*O(1)] ;
mat_scanner_read_z(2,:) = [0 pixdim(2) 0 -pixdim(2)*O(2)] ;
mat_scanner_read_z(3,:) = [0 0 pixdim(3) -pixdim(3)*O(3)] ;
mat_scanner_read_z(4,:) = [0 0 0 1] ;

% saving rotated magnitude and phase - inputs for clearswi:
mag_rot_merged_file = fullfile(output_fulldir, 'mag_rot.nii') ;
ph_rot_merged_file = fullfile(output_fulldir, 'ph_rot.nii') ;
createNifti(abs(cmplx_rot), mag_rot_merged_file, mat_scanner_read_z)
createNifti(angle(cmplx_rot), ph_rot_merged_file, mat_scanner_read_z)
clear cmplx_rot real_rot imag_rot

%%
disp('CLEAR SWI calculation')
if para.data_cleanup
    status = system(sprintf('%s -p %s -m %s -t [%s] --filter-size [%s] --mag-sensitivity-correction %s', para.clearswi_command, ph_rot_merged_file, mag_rot_merged_file, num2str(TEs), num2str(para.filter_size), para.sensitivity_corr)) ;
else
    status = system(sprintf('%s -p %s -m %s -t [%s] --filter-size [%s] --mag-sensitivity-correction %s --writesteps %s', para.clearswi_command, ph_rot_merged_file, mag_rot_merged_file, num2str(TEs), num2str(para.filter_size), para.sensitivity_corr, output_fulldir)) ;
end
if status == 1
    error('CLEARSWI did not run properly - check your installation path')
end

% rotating SWI and MIP to original MPM space
clearswi = nifti('clearswi.nii').dat(:,:,:);
mip = nifti('mip.nii').dat(:,:,:);
clearswi = ipermute(clearswi,[3 1 2]);
mip = ipermute(mip,[3 1 2]);
clearswi_readx_file = fullfile(output_fulldir, 'clearswi_readx.nii');
mip_readx_file = fullfile(output_fulldir, 'mip_readx.nii');
createNifti(clearswi, clearswi_readx_file, mat_scanner_read_x)
createNifti(mip, mip_readx_file, mat_scanner_read_x)

clearswi_V = spm_vol(clearswi_readx_file);
mip_V = spm_vol(mip_readx_file);
mip_rot = zeros(data_dim(1:3));
clearswi_rot = zeros(data_dim(1:3));
for slice = 1 : data_dim(3)
    clearswi_rot(:,:,slice) = spm_slice_vol(clearswi_V, inv(img2scanner_mat)*spm_matrix([0 0 slice]), data_dim_xy, -7) ;
    mip_rot(:,:,slice) = spm_slice_vol(mip_V, inv(img2scanner_mat)*spm_matrix([0 0 slice]), data_dim_xy, -7) ;
end
createNifti(clearswi_rot,'clearswi_invrot.nii', mat_image);
createNifti(mip_rot,'mip_invrot.nii', mat_image);

% data cleanup if interested only in the final result
if para.data_cleanup
    delete(real_merged_file,imag_merged_file,mag_rot_merged_file,ph_rot_merged_file, clearswi_readx_file, mip_readx_file)
end