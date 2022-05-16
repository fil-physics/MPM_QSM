%%% Description: MPM QSM pipeline
% main steps:
% 1) complex-fit over echoes for pdw and t1w images,
%    simple phase difference for mtw images
%    for odd and even echoes done separately
% 2) ROMEO phase unwrapping
% 3) masking based on ROMEO quality map
% 4) rotation to scanner space
% 5) PDF background field removal
% 6) star QSM for dipole inversion as default (optional: non-linear dipole inversion)


%%% Publications:
% Please remember to give credit to the authors of the methods used:
% 1. SEPIA toolbox:
% Chan, K.-S., Marques, J.P., 2021. Neuroimage 227, 117611.
% 2. SPM12 - rigid body registration:
% Friston KJ, et al. Magnetic Resonance in Medicine 35 (1995):346-355
% 3. complex fit of the phase:
% Liu, Tian, et al. MRM 69.2 (2013): 467-476.
% 4. ROMEO phase uwnrapping:
% Dymerska, Barbara, and Eckstein, Korbinian et al. Magnetic Resonance in Medicine (2020).
% 5. PDF background field removal:
% Liu, Tian, et al. NMR in Biomed. 24.9 (2011): 1129-1136.
% 6. starQSM:
% Wei, Hongjiang, et al. NMR in Biomed. 28.10 (2015): 1294-1303.

%%% Inputs:
% romeo_command          : path to romeo phase uwnrapping followed by romeo command, i.e. (in linux) '/your_path/bin/romeo' or (in windows) 'D:\your_path\bin\romeo'
% in_root_dir            : root directory to input nifti files
% out_root_dir           : root directory to output nifti files
% B0                     : magnetic field strength, in Tesla
% dipole_inv             : dipole inversion method, either 'Star-QSM' or 'ndi'
%                          'ndi'      - non-linear dipole inversion
%                                       (also known as iterative Tikhonov),
%                                       may give more contrast than Star-QSM but is less robust to noise
%                          'Star-QSM' - is very robust to noise and quick

%%%% Inputs - directories, parameters and files specific to given contrast
% ATTENTION: ensure only niftis you want to use are in that folder, with increasing echo numbering:
% mag_dir                : % folder with magnitude niftis
% ph_dir                 : % folder with phase inftis
% TEs                    : % echo time in ms
% output_dir             : % output QSM directory for a specific MPM contrast
% calc_mean_qsm          : % 'yes' or 'no' , if 'yes' it calculates mean QSM from all contrasts

%%% Outputs:
%%%% combined final results in out_root_dir:
% QSM_all_mean.nii             : mean QSM over all contrasts in scanner space (3rd dimension is along B0-axis)
% QSM_all_invrot_mean.nii      : mean QSM over all contrasts in image space (as acquired, for comparison with MPM quantitative maps)
% QSM_pdw_t1w_mean.nii         : mean QSM over PDw and T1w contrasts (without noisy MTw) in scanner space
% QSM_pdw_t1w_invrot_mean.nii  : mean QSM over PDw and T1w contrasts in image space

%%%% final results - per contrast in subfolders in out_root_dir:
% sepia_QSM.nii             : QSM in scanner space
% sepia_QSM_invrot.nii      : QSM in image space

%%%% additional outputs:
% ph.nii                    : two volumes (odd and even) of fitted phase
% ph_romeo.nii              : ph.nii unwrapped with ROMEO
% quality.nii               : quality map calculated by ROMEO algorithm and used for masking
% mask.nii                  : binary mask in image space
% mask_rot.nii              : binary mask in scanner space
% B0.nii                    : field map in Hz in image space
% B0_rot.nii                : field map in Hz in scanner space
% sepia_local-field.nii.gz  : map of local field variations (after background field removal using PDF)
% settings_romeo.txt        : settings used for ROMEO unwrapping (useful if unwrapping again outside MPM QSM the pipeline)
% header_sepia.mat          : header used for SEPIA toolbox (useful when exploring SEPIA GUI)

% script created by Barbara Dymerska
% @ UCL FIL Physics

totstart = tic ;

%%%%% USER PARAMETERS %%%%%
para.romeo_command = '/your_path/romeo_linux_3.2.0/bin/romeo' ;
para.in_root_dir = '/your/root/path' ;
para.out_root_dir =   '/your/output/path';

para.B0 = 3;
para.B0_dir = [0;1;0];	% main magnetic field direction after reslicing the data
para.dipole_inv = 'Star-QSM' ;

para.data_cleanup = 'big' ; % 'small' leaves B0 maps and QSMs, 'big' leaves only QSMs

for run = 1:3
    
    switch run
        case 1 % PDw
            para.mag_dir = 'pdw/mag/' ; % folder with magnitude niftis
            para.ph_dir = 'pdw/ph/' ; % folder with phase inftis
            para.TEs =  2.3*[1 2 3 4 5 6 7 8] ;  % echo time in ms
            para.output_dir = 'QSM_pdw' ; % output directory for a specific submeasurement from MPM
            para.mask_thr = 0.15 ; % larger threshold smaller mask
            
        case 2 % T1w
            para.mag_dir = 't1w/mag/' ;
            para.ph_dir = 't1w/ph/' ;
            para.TEs = 2.3*[1 2 3 4 5 6 7 8] ;
            para.output_dir = 'QSM_t1w' ;
            para.mask_thr = 0.15 ; 
            
        case 3 % MTw
            para.mag_dir = 'mtw/mag/' ;
            para.ph_dir = 'mtw/ph/' ;
            para.TEs = 2.3*[1 2 3 4 5 6] ; 
            para.output_dir = 'QSM_mtw' ;
            para.mask_thr = 0.1 ;
 
    end
    %%%%% END OF USER PARAMETERS %%%%%
    
MPM_QSM(para) ;
    
end

sprintf('total processing finished after %s' , secs2hms(toc(totstart)))
clear