%%% Description: MPM QSM pipeline
% main steps:
%  1) phase unwrapping and B0 map calculation using ROMEO
%  2) masking based on ROMEO quality map
%  3) rotation to scanner space for oblique acquisitions using SPM
%  4) PDF background field removal within SEPIA toolbox
%  5) star QSM for dipole inversion as default (optional: non-linear dipole inversion) within SEPIA toolbox
%  6) rotation back of QSM results to image space (for comparisons with PD, R2*, R1 and MT maps) using SPM (optional: non-linear dipole inversion)


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
% sepia_QSM.nii OR sepia_Chimap.nii            : QSM in scanner space (name depends on SEPIA toolbox version)
% sepia_QSM_invrot.nii                         : QSM in image space

%%%% additional outputs:
% ph.nii                    : two volumes (odd and even) of fitted phase
% ph_romeo.nii              : ph.nii unwrapped with ROMEO
% quality.nii               : quality map calculated by ROMEO algorithm and used for masking
% mask.nii                  : binary mask in image space
% mask_rot.nii              : binary mask in scanner space
% B0.nii                    : field map in Hz in image space
% B0_rot.nii                : field map in Hz in scanner space
% sepia_local-field.nii.gz OR  sepia_localfield.nii.gz : map of local field variations (after background field removal using PDF)
% settings_romeo.txt        : settings used for ROMEO unwrapping (useful if unwrapping again outside the MPM_QSM pipeline)
% header_sepia.mat          : header used for SEPIA toolbox (useful when exploring SEPIA GUI)

% script created by Barbara Dymerska
% @ UCL FIL Physics

function [QSM_V, QSM , QSMinvrot_V, QSMinvrot] = MPM_QSM(para)

tstart = tic ;

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
    
    ph_merged_file = fullfile(output_fulldir, 'ph.nii') ;
    mag_merged_file = fullfile(output_fulldir, 'mag.nii') ;

    createNifti(ph, ph_merged_file, ph_1tp.mat)
    createNifti(mag, mag_merged_file, mag_1tp.mat)
    else 
        ph_merged_file = fullfile(ph_fulldir, para.ph_file) ;
        mag_merged_file = fullfile(mag_fulldir, para.mag_file) ;
    end

    disp('phase unwrapping with ROMEO and:')
    disp('...removing global mean value')
    disp('......field map calculation')
    disp('.........saving quality map for masking')
    ph_romeo_file = fullfile(output_fulldir, 'ph_romeo.nii') ;

    status = system(sprintf('%s -p %s -m %s -o %s -t [%s] -q -B --weights 111111 --phase-offset-correction bipolar', para.romeo_command, ph_merged_file, mag_merged_file, ph_romeo_file, num2str(TEs))) ;

    if status == 1
        error('ROMEO did not run properly - check your installation path')
    end
    
    %% field map rotation 
    % defining affine matrix in scanner space for data rotation to scanner
    % space with mantaining the same image origin (i.e. no translation)
    FM_V = spm_vol('B0.nii') ;
    FM = nifti('B0.nii');
    disp('field map rotation to scanner space')
    data_dim = size(FM.dat) ;
    Z = spm_imatrix(FM.mat) ;
    pixdim = Z(7:9);
    
    mat_image = FM.mat ;
    O = mat_image\[0 0 0 1]' ;
    O = O(1:3)' ;
    
    mat_scanner(1,:) = [0 0 pixdim(3) -pixdim(3)*O(3)] ;
    mat_scanner(2,:) = [pixdim(1) 0 0 -pixdim(1)*O(1)] ;
    mat_scanner(3,:) = [0 pixdim(2) 0 -pixdim(2)*O(2)] ;
    mat_scanner(4,:) = [0 0 0 1] ;
    
    
    img2scanner_mat = mat_image\mat_scanner ;
    FMrot = zeros(data_dim) ;
    data_dim_xy = data_dim(1:2);
    
    for slice = 1 : data_dim(3)
        FMrot(:,:,slice) = spm_slice_vol(FM_V, img2scanner_mat*spm_matrix([0 0 slice]), data_dim_xy, -7) ;
    end
    
    FMrot(~isfinite(FMrot)) = 0 ;
    createNifti(FMrot, 'B0_rot.nii', mat_scanner)
    
    %% creating mask for QSM calculation
    
    disp('quality masking')
    qmask = nifti('quality.nii') ;
    qmask = qmask.dat(:,:,:) ;
    qmask(~isfinite(qmask)) = 0;
    qmask(qmask>para.mask_thr) = 1 ;
    qmask(qmask<=para.mask_thr) = 0 ;
    % filling holes in the mask
    qmask = imfill(qmask,6,'holes') ;
    qmask = smooth3(qmask, 'gaussian') ;
    qmask(qmask>0.6) = 1 ;
    qmask(qmask<=0.6) = 0 ;

    qmask = int16(qmask) ;
    createNifti(qmask, 'mask.nii', mat_image)
    clear qmask
    
    mask_V = spm_vol('mask.nii') ;
    qmask_rot = zeros(data_dim) ;
    for slice = 1 : data_dim(3)
        qmask_rot(:,:,slice) = spm_slice_vol(mask_V, img2scanner_mat*spm_matrix([0 0 slice]), data_dim_xy, -7) ;
    end
    qmask_rot(~isfinite(qmask_rot)) = 0 ;
    createNifti(int16(qmask_rot), 'mask_rot.nii', mat_scanner)    
    
    %% SEPIA - background field removal and dipole inversion yielding final QSM
    
    disp('creating SEPIA header')
    B0 = para.B0 ;
    CF = B0*42.58*1e6;	% imaging frequency, in Hz (B0*gyromagnetic_ratio*1e6)
    delta_TE = 1;	    % echo spacing, in second - we have already combined data, in such situation set to 1
    TE = 1 ;
    B0_dir = para.B0_dir;	% main magnetic field direction, it's always [0,1,0] because the images are resliced so that 2nd dimension is aligned with B0
    matrixSize = data_dim ;	    % image matrix size
    voxelSize = pixdim ;	% spatial resolution of the data, in mm
    header_fullfile = fullfile(output_fulldir, 'header_sepia.mat') ;
    save(header_fullfile, 'B0', 'B0_dir', 'CF', 'TE', 'delta_TE', 'matrixSize', 'voxelSize')
    
    % general SEPIA parameters
    sepia_addpath
    algorParam.general.isBET = 0 ;
    algorParam.general.isInvert = 1 ;
    algorParam.general.isGPU = 0 ;
    output_basename = fullfile(output_fulldir, 'sepia') ;
    
    % inputs for background field removal
    input(1).name = 'B0_rot.nii' ;
    input(2).name = 'mask_rot.nii' ;
    input(4).name = header_fullfile ;
    
    algorParam.bfr.refine = 0 ;
    algorParam.bfr.erode_radius = 1 ;
    algorParam.bfr.method = 'pdf' ;
    algorParam.bfr.tol = 0.1 ;
    algorParam.bfr.iteration = 50 ;
    algorParam.bfr.padSize = 30 ;
    
    
    % inputs for dipole inversion
    algorParam.qsm.method = para.dipole_inv ;
    if strcmp(algorParam.qsm.method , 'ndi')
        
        algorParam.qsm.method = 'ndi' ;
        algorParam.qsm.tol = 1 ;
        algorParam.qsm.maxiter = 200 ;
        algorParam.qsm.stepSize = 1 ;
        
    elseif strcmp(algorParam.qsm.method , 'Star-QSM')
        
        algorParam.qsm.padsize = ones(1,3)*6 ;
        
    end
    
    
    disp('background field removal using PDF')
    BackgroundRemovalMacroIOWrapper(input,output_basename,input(2).name,algorParam);
    
    % added for back-compatibility to older SEPIA versions
    if exist('sepia_local-field.nii.gz','file') == 2
        input(1).name = 'sepia_local-field.nii.gz' ;
    elseif exist('sepia_localfield.nii.gz','file') == 2
        input(1).name = 'sepia_localfield.nii.gz' ;
    else
        error('no local field file found in output dir')
    end

    fprintf('dipole inversion using %s', algorParam.qsm.method)
    QSMMacroIOWrapper(input,output_basename,input(2).name,algorParam);
    
    
    disp('rotation of QSM back to the original image space')
    % added for back-compatibility to older SEPIA versions
    if exist('sepia_QSM.nii.gz') == 2
        gunzip('sepia_QSM.nii.gz')
        QSM = nifti('sepia_QSM.nii') ;
        QSM_V = spm_vol('sepia_QSM.nii');
    elseif exist('sepia_Chimap.nii.gz') == 2
        gunzip('sepia_Chimap.nii.gz')
        QSM = nifti('sepia_Chimap.nii') ;
        QSM_V = spm_vol('sepia_Chimap.nii');
    else
        error('no QSM maps in output dir')
    end

    QSM = QSM.dat(:,:,:) ;
    
    scanner2img_mat = mat_scanner\mat_image ;
    QSMinvrot = zeros(data_dim) ;
    
    for slice = 1 : data_dim(3)
        QSMinvrot(:,:,slice) = spm_slice_vol(QSM_V, scanner2img_mat*spm_matrix([0 0 slice]), data_dim_xy, -7) ;
    end
    
    QSMinvrot_V = FM_V ;
    QSMinvrot_V.fname = 'sepia_QSM_invrot.nii';
    spm_write_vol(QSMinvrot_V, QSMinvrot);
    
    warning('off'); 
    delete sepia_mask-qsm.nii.gz sepia_QSM.nii.gz sepia_mask_QSM.nii.gz sepia_Chimap.nii.gz
    if strcmp(para.data_cleanup,'small') || strcmp(para.data_cleanup,'big')
        delete mag.nii corrected_phase.nii mask.nii mask_rot.nii ph.nii ph_romeo.nii quality.nii sepia_local-field.nii.gz sepia_localfield.nii.gz
    end
    if strcmp(para.data_cleanup,'big')
        delete B0.nii B0_rot.nii
    end
    warning('on');
    sprintf('finished after %s' , secs2hms(toc(tstart)))
end


