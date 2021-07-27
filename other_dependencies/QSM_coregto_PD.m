function QSM_coregto_PD(ref, source, other)
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {'/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_0009/s2021-06-23_10-18-103458-00001-00288-1.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {'/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_0025/s2021-06-23_10-18-112654-00001-00288-1.nii,1'};
%%
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_0025/s2021-06-23_10-18-112654-00001-00576-2.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_0025/s2021-06-23_10-18-112654-00001-00864-3.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_0025/s2021-06-23_10-18-112654-00001-01152-4.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_0025/s2021-06-23_10-18-112654-00001-01440-5.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_0025/s2021-06-23_10-18-112654-00001-01728-6.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0044_scaled10/s2021-06-23_10-18-112654-00001-00288-1.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0044_scaled10/s2021-06-23_10-18-112654-00001-00576-2.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0044_scaled10/s2021-06-23_10-18-112654-00001-00864-3.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0044_scaled10/s2021-06-23_10-18-112654-00001-01152-4.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0044_scaled10/s2021-06-23_10-18-112654-00001-01440-5.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0044_scaled10/s2021-06-23_10-18-112654-00001-01728-6.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0058/s2021-06-23_10-18-112654-00001-00288-1.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0058/s2021-06-23_10-18-112654-00001-00576-2.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0058/s2021-06-23_10-18-112654-00001-00864-3.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0058/s2021-06-23_10-18-112654-00001-01152-4.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0058/s2021-06-23_10-18-112654-00001-01440-5.nii,1'
                                                   '/media/barbara/hdd2/DATA/FIL/7T/20210623.M700198_FIL_analysis/pdw_mfc_3dflash_v1k_RR_0058/s2021-06-23_10-18-112654-00001-01728-6.nii,1'
                                                   };
%%
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 7;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
