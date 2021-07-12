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
% 2. complex fit of the phase:
% Liu, Tian, et al. MRM 69.2 (2013): 467-476.
% 3. ROMEO phase uwnrapping:
% Dymerska, Barbara, and Eckstein, Korbinian et al. Magnetic Resonance in Medicine (2020).
% 4. PDF background field removal:
% Liu, Tian, et al. NMR in Biomed. 24.9 (2011): 1129-1136.
% 5. starQSM:
% Wei, Hongjiang, et al. NMR in Biomed. 28.10 (2015): 1294-1303.

%%% Inputs:
% romeo_command          : path to romeo phase uwnrapping followed by romeo command, i.e. (in linux) '/your_path/bin/romeo' or (in windows) 'D:\your_path\bin\romeo'
% B0                     : magnetic field strength, in Tesla
% algorParam.qsm.method  : dipole inversion method, either 'Star-QSM' or 'ndi'
%                          'ndi' - non-linear dipole inversion (also known as iterative Tikhonov), may give more contrast than Star-QSM but is less robust to noise
%                          'Star-QSM' is very robust to noise and quick
% in_root_dir            : root directory to input nifti files
% out_root_dir           : root directory to output nifti files
%%%% Inputs - directories, parameters and files specific to given contrast
% mag_dir                : % folder with magnitude niftis
% ph_dir                 : % folder with phase inftis
% TEs                    : % echo time in ms
% output_dir             : % output directory for a specific submeasurement from MPM
% mag_file               : % magnitude reference nifti file for ROMEO unwrapping and masking


%%% Outputs:
%%%% combined final results in out_root_dir:
% QSM_all_mean.nii             : mean QSM over all contrasts in scanner space (3rd dimension is along B0-axis)
% QSM_all_invrot_mean.nii      : mean QSM over all contrasts in image space (as acquired, for comparison with MPM quantitative maps)
% QSM_pdw_t1w_mean.nii         : mean QSM over PDw and T1w contrasts (without noisy MTw) in scanner space
% QSM_pdw_t1w_invrot_mean.nii  : mean QSM over PDw and T1w contrasts in image space

%%%% final results - per contrast in subfolders in out_root_dir:
% sepia_QSM.nii.gz          : QSM in scanner space
% sepia_QSM_invrot.nii.gz   : QSM in image space

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
% last modifications 12/07/2021

tstart = tic ;
%%%%% USER PARAMETERS %%%%%

romeo_command = 'C:\wtcnapps\romeo_win_3.2.0\bin\romeo' ;
B0 = 7;
algorParam.qsm.method = 'Star-QSM' ;

in_root_dir = 'D:\Users\bdymerska\data\7T\2021\20210603.M700159_FIL_analysis' ;
out_root_dir = 'D:\Users\bdymerska\data\7T\2021\20210603.M700159_FIL_analysis\SEPIA';

% directories, parameters and files specific to given contrast:
for run = 1:3
    
    switch run
        case 1 %pdw
            mag_dir = '59' ; % folder with magnitude niftis
            ph_dir = '60' ; % folder with phase inftis
            TEs = [2.2 4.58 6.96 9.34 11.72 14.1] ; % echo time in ms
            output_dir = 'pdw_RR_59_60' ; % output directory for a specific submeasurement from MPM
            mag_file = 'sM700159-0059-00001-001728-06-001.nii' ; % magnitude reference nifti file for ROMEO unwrapping and masking
            
        case 2 % t1w
            mag_dir = '57' ;
            ph_dir = '58' ;
            TEs = [2.3 4.68 7.06 9.44 11.82 14.2] ;
            output_dir = 't1w_RR_57_58' ;
            mag_file = 'sM700159-0057-00001-001728-06-001.nii' ;
            
        case 3 % mtw
            mag_dir = '55' ;
            ph_dir = '56' ;
            TEs = [2.2 4.58 6.96 9.34] ; % echo time in ms
            output_dir = 'mtw_RR_55_56' ;
            mag_file = 'sM700159-0055-00001-001152-04-001.nii' ; % magnitude reference nifti file for ROMEO unwrapping and masking
    end
    
    
    %%%%% END OF USER PARAMETERS %%%%%
    
    mag_fulldir = fullfile(in_root_dir, mag_dir) ;
    ph_fulldir = fullfile(in_root_dir, ph_dir) ;
    
    output_fulldir = fullfile(out_root_dir, output_dir) ;
    if ~exist(output_fulldir, 'dir')
        mkdir(output_fulldir)
    end
    cd(output_fulldir)
    
    TEs = TEs/10^3 ;
    
    ph_files = dir(ph_fulldir);
    mag_files = dir(mag_fulldir);
    
    for t = 1:size(TEs,2)
        
        ph_1tp = load_untouch_nii(fullfile(ph_fulldir, ph_files(t+2).name));
        ph(:,:,:,t) = ph_1tp.img ;
        
        if size(TEs,2) >= 6 % for mtw acquisition we cannot perform complex fit so we don't need magnitude data for each TE
            
            mag_1tp = load_untouch_nii(fullfile(mag_fulldir, mag_files(t+2).name)).img ;
            mag(:,:,:,t) = mag_1tp ;
            
        end
        
    end
    
    
    clear ph_1tp.img mag_1tp
    
    % rescaling phase into [0,2*pi] phase range
    ph = 2*pi*single(ph - min(vector(ph)))/single(max(vector(ph))-min(vector(ph))) ;
    
    for read_dir = 1:2
        
        if size(TEs,2) < 6 % complex fit is only possible if at least 3 echoes per readout direction available
            disp('calculating phase difference')
            FM = angle(exp(1i*(ph(:,:,:,read_dir+2)-ph(:,:,:,read_dir)))) ;
        else
            disp('complex fitting phase')
            compl = single(mag).*exp(-1i*ph);
            [FM, ~, ~, ~] = Fit_ppm_complex_TE(compl(:,:,:,read_dir:2:end),TEs(read_dir:2:end));
        end
        
        FM_both(:,:,:,read_dir) = FM ;
        
    end
    clear mag ph FM
    
    % saving odd and even echoes as one file for ROMEO unwrapping
    FM_file = 'ph.nii' ;
    FM_both = make_nii(FM_both) ;
    FM_both.hdr.hist = ph_1tp.hdr.hist ;
    centre_and_save_nii(FM_both, FM_file, ph_1tp.hdr.dime.pixdim);
    clear FM_both
    
    % creating magnitude reference for ROMEO unwrapping
    mag_fullfile = fullfile(mag_fulldir, mag_file);
    mag_ref = load_untouch_nii(mag_fullfile) ;
    mag_ref = repmat(mag_ref.img,[1 1 1 2]) ;
    mag_fullfile = sprintf('mag_TE%i.nii',size(TEs,2)) ;
    centre_and_save_nii(make_nii(mag_ref), mag_fullfile, ph_1tp.hdr.dime.pixdim);
    clear mag_ref
    
    disp('phase unwrapping with ROMEO and:')
    disp('...removing global mean value')
    disp('......field map calculation')
    disp('.........saving quality map for masking')
    TE = (TEs(3)-TEs(1))*10^3; % effective echo time difference after phase complex fitting in seconds
    [~, FM_name,~] = fileparts(FM_file) ;
    FM_romeo_file = sprintf('%s_romeo.nii',FM_name) ;
    if isunix
        unix(sprintf('%s %s -m %s -o %s -t [%i,%i] -k nomask -g -q -B', romeo_command, FM_file, mag_fullfile, FM_romeo_file, TE, TE)) ;
    elseif ispc
        system(sprintf('%s %s -m %s -o %s -t [%i,%i] -k nomask -g -q -B', romeo_command, FM_file, mag_fullfile, FM_romeo_file, TE, TE)) ;
    end
    delete('mag_TE6.nii');
    reslice_nii('B0.nii', 'B0_rot.nii',ph_1tp.hdr.dime.pixdim(2:4), 1, 0)
    FM = load_nii('B0_rot.nii') ;
    FM.img(isnan(FM.img)) = 0;
    FM.img = changeImageSize(FM.img, circshift(ph_1tp.hdr.dime.dim(2:4),1)) ;
    FM.hdr.dime.dim(2:4) = circshift(ph_1tp.hdr.dime.dim(2:4),1) ;
    centre_and_save_nii(FM, 'B0_rot.nii', ph_1tp.hdr.dime.pixdim);
    clear FM
    
    disp('quality masking')
    qmap = load_untouch_nii('quality.nii') ;
    qmap_bin = qmap.img ;
    qmap_bin(qmap.img>0.3) = 1 ;
    qmap_bin(qmap.img<=0.3) = 0 ;
    clear qmap
    qmap_bin(isnan(qmap_bin)) = 0 ;
    qmap_bin = imfill(qmap_bin,6,'holes') ;
    qmask = smooth3(qmap_bin, 'gaussian') ;
    qmask(qmask>0.6) = 1 ;
    qmask(qmask<=0.6) = 0 ;
    clear qmap_bin
    
    qmask = make_nii(int16(qmask)) ;
    qmask.hdr.hist = ph_1tp.hdr.hist ;
    qmask.hdr.dime.pixdim = ph_1tp.hdr.dime.pixdim ;
    qmask_file = fullfile(output_fulldir, 'mask.nii') ;
    save_nii(qmask, qmask_file);
    reslice_nii('mask.nii', 'mask_rot.nii', ph_1tp.hdr.dime.pixdim(2:4), 1 , 0)
    qmask_file = fullfile(output_fulldir, 'mask_rot.nii') ;
    qmask = load_nii(qmask_file) ;
    qmask.img = round(qmask.img) ;
    
    qmask.img = changeImageSize(qmask.img, circshift(ph_1tp.hdr.dime.dim(2:4),1)) ;
    qmask.hdr.dime.dim(2:4) = circshift(ph_1tp.hdr.dime.dim(2:4),1) ;
    centre_and_save_nii(qmask, qmask_file, ph_1tp.hdr.dime.pixdim);
    clear qmask
    
    %% SEPIA - calculates QSM
    
    % create SEPIA header
    CF = B0*42.58*1e6;	% imaging frequency, in Hz (B0*gyromagnetic_ratio*1e6)
    delta_TE = 1;	    % echo spacing, in second - we have already combined data, in such situation set to 1
    B0_dir = [0;0;1];	% main magnetic field direction, it's always [0,0,1] because the images are resliced so that 3rd dimension is aligned with B0
    hdr = load_nii_hdr('B0_rot.nii') ;
    matrixSize = hdr.dime.dim(2:4) ;	    % image matrix size
    voxelSize = ph_1tp.hdr.dime.pixdim(2:4) ;	% spatial resolution of the data, in mm
    header_fullfile = fullfile(output_fulldir, 'header_sepia.mat') ;
    save(header_fullfile, 'B0', 'B0_dir', 'CF', 'TE', 'delta_TE', 'matrixSize', 'voxelSize')
    
    % general SEPIA parameters
    sepia_addpath
    
    algorParam.general.isBET = 0 ;
    algorParam.general.isInvert = 1 ;
    algorParam.general.isGPU = 0 ;
    
    % inputs for background field removal
    input(1).name = 'B0_rot.nii' ;
    input(2).name = 'mask_rot.nii' ;
    input(4).name = header_fullfile ;
    
    algorParam.bfr.refine = 0 ;
    algorParam.bfr.erode_radius = 0 ;
    algorParam.bfr.method = 'pdf' ;
    algorParam.bfr.tol = 0.1 ;
    algorParam.bfr.iteration = 50 ;
    algorParam.bfr.padSize = 30 ;
    
    
    % inputs for dipole inversion
    if strcmp(algorParam.qsm.method , 'ndi')
        
        algorParam.qsm.method = 'ndi' ;
        algorParam.qsm.tol = 1 ;
        algorParam.qsm.maxiter = 200 ;
        algorParam.qsm.stepSize = 1 ;
        
    elseif strcmp(algorParam.qsm.method , 'Star-QSM')
        
        algorParam.qsm.padsize = ones(1,3)*12 ;
        
    end
    
    output_basename = fullfile(output_fulldir, 'sepia') ;
    
    disp('background field removal')
    BackgroundRemovalMacroIOWrapper(input,output_basename,input(2).name,algorParam);
    
    disp('dipole inversion')
    input(1).name = fullfile(output_fulldir, 'sepia_local-field.nii.gz') ;
    QSMMacroIOWrapper(input,output_basename,input(2).name,algorParam);
    
    sprintf('run %i finished after %s' ,run, secs2hms(toc(tstart)))
    
    QSM = load_nii(fullfile(output_fulldir, 'sepia_QSM.nii.gz'));
    
    disp('rotation of QSM back to the original image space')
    QSM_invrot = flip(permute(QSM.img, [2 3 1]),1);
    alpha = rad2deg(-asin(ph_1tp.hdr.hist.srow_y(2)/ph_1tp.hdr.dime.pixdim(2))) ;
    
    M_aff(1,:) = [cos(deg2rad(alpha))    -sin(deg2rad(alpha))     0        0    ];
    M_aff(2,:) = [sin(deg2rad(alpha))    cos(deg2rad(alpha))      0        0    ];
    M_aff(3,:) = [0                      0                        1.0000   0    ];
    M_aff(4,:) = [0                      0                        0        1.0000];
    
    [QSM_invrot, ~] = affine(QSM_invrot, M_aff);
    QSM_invrot = changeImageSize(QSM_invrot, ph_1tp.hdr.dime.dim(2:4)) ;
    QSM_invrot = make_nii(QSM_invrot) ;
    QSM_invrot.hdr.hist = ph_1tp.hdr.hist ;
    centre_and_save_nii(QSM_invrot, 'sepia_QSM_invrot.nii.gz', ph_1tp.hdr.dime.pixdim)
    
    QSM_all(:,:,:,run) = QSM.img ;
    QSM_all_invrot(:,:,:,run) = QSM_invrot ;
    clear QSM QSM_invrot
    delete('sepia_mask-qsm.nii.gz')
    
end

QSM_all_mean = mean(QSM_all, 4) ;
QSM_pdw_t1w_mean = mean(QSM_all(:,:,:,1:2), 4) ;

QSM_all_invrot_mean = mean(QSM_all_invrot, 4) ;
QSM_pdw_t1w_invrot_mean = mean(QSM_all_invrot(:,:,:,1:2), 4) ;

QSM_all_mean = make_nii(QSM_all_mean) ;
QSM_all_mean.hdr.hist = ph_1tp.hdr.hist ;

QSM_pdw_t1w_mean = make_nii(QSM_pdw_t1w_mean) ;
QSM_pdw_t1w_mean.hdr.hist = ph_1tp.hdr.hist ;

QSM_all_invrot_mean = make_nii(QSM_all_invrot_mean) ;
QSM_all_invrot_mean.hdr.hist = ph_1tp.hdr.hist ;

QSM_pdw_t1w_invrot_mean = make_nii(QSM_pdw_t1w_invrot_mean) ;
QSM_pdw_t1w_invrot_mean.hdr.hist = ph_1tp.hdr.hist ;

centre_and_save_nii(QSM_all_mean, fullfile(out_root_dir,'QSM_all_mean.nii'), ph_1tp.hdr.dime.pixdim);
centre_and_save_nii(QSM_pdw_t1w_mean, fullfile(out_root_dir,'QSM_pdw_t1w_mean.nii'), ph_1tp.hdr.dime.pixdim);

centre_and_save_nii(QSM_all_invrot_mean, fullfile(out_root_dir,'QSM_all_invrot_mean.nii'), ph_1tp.hdr.dime.pixdim);
centre_and_save_nii(QSM_pdw_t1w_invrot_mean, fullfile(out_root_dir,'QSM_pdw_t1w_invrot_mean.nii'), ph_1tp.hdr.dime.pixdim);


sprintf('total processing finished after %s' , secs2hms(toc(tstart)))
