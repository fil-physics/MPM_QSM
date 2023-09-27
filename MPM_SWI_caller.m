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

totstart = tic ;

%%%%% USER PARAMETERS %%%%%
para.clearswi_command = '/your_path/mritools_Linux_3.6.4/bin/clearswi' ;
para.in_root_dir = '/your/root/path' ;
para.out_root_dir =   '/your/output/path';
para.sensitivity_corr = 'on' ; % corrects for bias field, 'on' or 'off'
para.filter_size = [2 2 0] ; % for 7T MPM [2 2 0] is adviced
para.data_cleanup = true ; % true - cleans intermediate outputs, false - leaves them

for run = 1
    
    switch run
        case 1 % PDw
            para.mag_dir = 'pdw/mag' ; % folder with magnitude niftis
            para.ph_dir = 'pdw/ph' ; % folder with phase inftis
            para.ph_file = '' ; % specify if data saved as 4D
            para.mag_file = ''; % specify if data saved as 4D
            para.TEs =  [2.2 4.58 6.96 9.34 11.72 14.1] ;  % echo time in ms
            para.output_dir = 'pdw/SWI' ; % output directory for a specific submeasurement from MPM
           
        case 2 % T1w
            para.mag_dir = 't1w/mag' ;
            para.ph_dir = 't1w/ph' ;
            para.ph_file = '' ; 
            para.mag_file = ''; 
            para.TEs = [2.3 4.68 7.06 9.44 11.82 14.2] ;
            para.output_dir = 't1w/SWI' ;
            
        case 3 % MTw
            para.mag_dir = 'mtw/mag' ;
            para.ph_dir = 'mtw/ph' ;
            para.ph_file = '' ; 
            para.mag_file = '';
            para.TEs = [2.2 4.58 6.96 9.34] ;
            para.output_dir = 'mtw/SWI' ;

    end
    %%%%% END OF USER PARAMETERS %%%%%
    
MPM_SWI(para) ;
    
end

sprintf('total processing finished after %s' , secs2hms(toc(totstart)))
clear