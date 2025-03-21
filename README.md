# MPM_QSM
QSM pipeline for Multi Parametric Mapping acquisitions

# MPM_SWI   
CLEAR-SWI pipeline for Multi Parametric Mapping acquisitions

# MPM_QSM main computational steps:

 1) phase unwrapping and B0 map calculation using ROMEO
 2) masking based on ROMEO quality map
 3) rotation to scanner space for oblique acquisitions using SPM
 4) PDF background field removal within SEPIA toolbox
 5) star QSM for dipole inversion as default (optional: non-linear dipole inversion) within SEPIA toolbox
 6) rotation back of QSM results to image space (for comparisons with PD, R2*, R1 and MT maps) using SPM

# MPM_SWI main computational steps:
 1) rotation of the magnitude and phase into scanner space
   (in the complex domain so that phase wraps are not blurred)
 2) CLEAR-SWI
   - SNR-optimal magnitude combination,
   - phase mask calculation and multiplication with combined magnitude
   - minimum intensity projection (MIP)
 3) rotation  of CLEAR-SWI and MIP back into image space

# Installation steps:

1. Download Zip with all the files from:

	https://github.com/fil-physics/MPM_QSM
	or clone it to your GitHub repository

2. Download compiled version of [ROMEO](https://github.com/korbinian90/CompileMRI.jl/releases/tag/v3.6.4) within MRItools either for windows or linux and unzip it in chosen destination

3. Download sepia toolbox:
	
	https://github.com/kschan0214/sepia.git
	and unzip it in chosen destination

4. Download SPM12:
	https://www.fil.ion.ucl.ac.uk/spm/software/download/

5. Set you local paths to MEDI and STI toolboxes downloaded in step 1 in file:
 	/your_path/sepia-master/SpecifyToolboxesDirectory.m

	as following:
	MEDI_dir = '/your/MEDI/path/';
	STISuite_dir = '/your/STI/path/';
	FANSI_dir = [];
	SEGUE_dir = [];

6. Add to your matlab path: SEPIA toolbox, MPM_QSM folder and SPM12
7. Edit MPM_QSM_caller.m user parameters, where you specify folders to you nifti files


# Publications:

Please remember to give credit to the authors of the methods used:

## For MPM_QSM:
1. MORSE-CODE used for image reconstruction:
   Oliver Josephs, Barbara Dymerska, et al. ???
2. SEPIA toolbox:
Chan KS, Marques JP. "SEPIA—Susceptibility mapping pipeline tool for phase images." Neuroimage 227 (2021), 117611.

3. SPM12 - rigid body registration:
Friston KJ, et al. "Movement-related effect in fMRI time-series." Magnetic Resonance in Medicine 35 (1995):346-355

4. complex fit of the phase:
Liu T, et al. "Nonlinear formulation of the magnetic field to source relationship for robust quantitative susceptibility mapping." Magnetic resonance in medicine 69.2 (2013): 467-476.

5. ROMEO phase unwrapping:
Dymerska B, and Eckstein K et al. "Phase unwrapping with a rapid opensource minimum spanning tree algorithm (ROMEO)." Magnetic Resonance in Medicine (2020).

6. PDF background field removal:
Liu T, et al. "A novel background field removal method for MRI using projection onto dipole fields." NMR in Biomedicine 24.9 (2011): 1129-1136.

7. starQSM:

Wei H, et al. "Streaking artifact reduction for quantitative susceptibility mapping of sources with large dynamic range." NMR in Biomedicine 28.10 (2015): 1294-1303.

## For MPM_SWI:
1. MORSE-CODE used for image reconstruction:
   Oliver Josephs, Barbara Dymerska, et al. ???
3. CLEARSWI used for SWI calculation:
   Eckstein, Korbinian, et al. "Improved susceptibility weighted imaging at ultra-high field using bipolar multi-echo acquisition and optimized image processing: CLEAR-SWI." Neuroimage 237 (2021): 118175.

