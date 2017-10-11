%HFR_DP_SETUP_README.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_SETUP_README.m
%
%  This package contains the master programs to process CODAR-type HF Radar
%  system data from the estimated spectra (the *.css or *.cs4) files to radial
%  radial velocities averaged over an azimuthal window for a particular 
%  COS radar site. 
%
%  This README file describes:
%  (1) what HFR_DP is, 
%  (2) what is needed to run HFR_DP, and
%  (3) what HFR_DP is not.
%
% This package distributed as is, with no additional warranties, as an open
% source community resource to enable HFR methodological development, and 
% is not meant to supplant existing COS (Codar Ocean Systems) software
% distributions.
%
% Similar to HFR_Progs, no part of this distrubution can be used for 
% commerical gain.
%
% v1        March 2017
%
% Anthony Kirincich
% Woods Hole Oceanographic Institution
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What HFR_DP does                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  The HFR_DP package is intended to be used as a platform to reprocess 
%  HFR radar data for direct scientific research and/or advancing HFR meth-
%  ological research.  The program and its structure
%  are based on the HFR_Progs matlab toolbox (by Kaplan and Cook) and the 
%  package contains versions of key HFR_Progs scripts that have been adapted
%  to provide additional output or functionality.  
%
%  The program follows many of the processing choices used by COS, as 
%  described by Lipa et al. (2006), for processing CODAR-type HF radar 
%  data sets. HFR_DP does not provide an exact replica of COS's processing 
%  methods (see below) but does provide additional features for how the radial 
%  metrics output can be used to apply additional quality control features
%  See Kirincich et al. (2012) or de Paola et al. (2015) for details.
%
%  HFR_DP can use either the COS-determined First Order Limits predetermined
%  by the COS Radial Suite or an alternative implementation of FOLs using image
%  processing techniques (Kirincich, 2017). 
%
%  HFR_DP performs the MUSIC direction finding algorithm on the spectral data
%  using the same parameter settings used by COS and found in the Site's 
%  header.txt file.
%
%  HFR_DP performs NO temporal averaging.  Thus, the end result of HFR_DP
%  is equivalent to COS-processed 'Radial Short' files.
%
%  There are numerous types of output products available.  All are available for 
%  each individual spectral estimate submitted for processsing.  Outputs
%  include:
%
%  - Radial Metrics output   (as a Matlab structure within a *.mat file)
%
%  - Radial Average output   (the spatial averaged result of the Radial 
%                             Metrics files, as a matlab structure within
%                             a *.mat file) 
%
%  - LLUV output             (the spatial averaged result of the Radial 
%                             Metrics files, saved in an ASCII text file 
%                             equivalent to COS lluv format.)
%
%  - *.jpg figures of the Radial metrics and Radial averaged results are
%                             also available.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What HFR_DP Needs                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  HFR_DP requires the following programs to be loaded and operational 
%  on the target machine:
%
%
%  (1) A MAC machine, though LINUX would likely work as well.
%
%  (2) MATLAB (tested on post-2013 releases only)
%
%  (3) The HFR_Progs and M_map packages, available at ROWG WIKI or github
%      site  (must be available within your path)
%
%  (4) The Matlab image processing toolbox is required to run ImageFOLs.  If
%      not using ImageFOLs, the input spectral files must have been previously 
%      used by COS Radial Suite to compute radial velocities, leaving the 
%      COS-determined FOLs within the resource fork of the spectral files.
%      To be clear, HFR_DP can only make NEW estimates of the FOLs using 
%      ImageFOL. 
%
%  (5) A subdirectory file structure for each site=XXXX is required, and 
%      must be formulated as such  (where ~ denotes some base directory):
%
%      ~/SITE_XXXX_css          %where the spectral files will be found
%      ~/SITE_XXXX_config       %where the header and pattern data will be found
%      ~/SITE_XXXX              %where the resulting 'Radial Metrics' output will go
%      ~/SITE_XXXX_radave       %where the resulting 'Radial Short' will go
%      ~/SITE_XXXX_radave_lluv  %where the resulting ascii file of the 'Radial Short'
%                               %  will go...suitable for transmission to the national archive
%
%      Any diversion from this structure and the user must alter the initial 
%      lines of HFR_DP_radial_prep_v?.m mat to compensate.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% What HFR_DP doesn't do                           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   HFR_DP does not provide an exact replica of COS's processing methods.
%   Specifically, HFR_DP does not:
%
%    (1)  Provide pre-processing noise cancelation.  All included noise
%         cancellation effects are included the imageFOL-based delineation
%         of what data is processed, but no formal de-striping, etc. is
%         done.  
%
%    (2)  Look for or account for ships, etc. within the spectral data, 
%
%    (3)  Perform any temporal averaging of the resulting radial velocities.
%         Thus, all resulting data, apart from the radial metrics themselves
%         is equivalent to COS's 'radial shorts' files.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% How to run HFR_DP                                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    (1) Add the HFR_Progs, M_map, and HFR_DP packages to your
%        matlab path.
%
%    (2) Create a working directory structure as shown above for each site
%        you intend to reprocess, or download the test data sets also
%        distributed by Kirincich (https://github.com/akirincich/HFR_DP_Data.git)
%
%    (3) Edit or resave HFR_DR_master_SITE.m with the information specific
%        to your SITE (i.e. LPWR), directories, and processing choices
%   
%    (4) Iterate as you need.
%
