%HFR_DP_master_XXXX.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_master_XXXX.m
%
%  This script contains the master program that runs the HFR_DP program of
%  spectra to radial average processing for CODAR-type HF radar datasets.
%
%  SETUP:  
%   To setup this template file or a particular COS radar site. 
%  The user must:
%  (1) define the site, array, directory locations, and processing 
%      choices in the first portion of the script, 
%  (2) save a new version of this template file with the XXXX replaced by
%      the site name
%  (3) Run the script within matlab, from the command line, or via a cron job to process
%      spectra CSS files to radial averaged files.
%
%  See HFR_DP_SETUP_README.m for directions on how to set up a machine to
%  perform the processing
%
% v1        March 2017
%
% Anthony Kirincich
% Woods Hole Oceanographic Institution
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     begin user choices       %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set site/array information

%%% the site name to process
site_name='LPWR';
%site_name='NWTP';

%%% set the common constants for the HFR array this site is a part of:
ARRAY=[];
ARRAY.name='MVCO';  
%%% array name, does not have to be four letters like radial sites, can be the
%%% same name for multiple radar sites that contribute to coverage within an area

%%% define the bounding box of the array, needed for any plotting operations.
ARRAY.min_lon=-(71+20/60); 
ARRAY.max_lon=-(69+40/60); 
ARRAY.min_lat=(40+40/60); 
ARRAY.max_lat=(41+30/60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set directory choices

%%% What is the working directory where the HFR_DP files are located
working_dir='/Users/anthony/Matlab/working/HFR_DP1';

%%% Where are the scripts
scripts_dir=[working_dir '/scripts/'];

%%% Where are the data  (see HFR_DP_SETUP_README.m for instructions
base_dir='/Users/anthony/Matlab/working/Data/RadialSites';
%base_dir='/Codar/SeaSonde/Data/RadialSites'; %if wish to place in codar directory

%%% Where do the DP logs go
log_dir=['/Users/anthony/Matlab/working/LOGS'];
%log_dir='/Codar/SeaSonde/Logs/codar_DP_Logs'; %if wish to place in codar directory


eval(['addpath ' scripts_dir]);
eval(['addpath ' working_dir]);
%%% move to inital base directory
eval(['cd ' working_dir])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set constants for the HFR
CONST=[];
%%% constants that won't change %%%
CONST.c=299792458;  %m/s speed of light in a vaccum
    %    c=299705000;   % in air   
    %    c=299644310;   % in air based on getting v=0 at ibragg (difference of 2 mm/s from first value....)
%number of antenna elements
CONST.N=3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set parameters for spectral processing

%%% the final azimuthal resolution of the radial velocities
CONST.bearing_width=5;

%%% the minimum Signal to Noise ratio allowed (used in Radial averaging and Image     
CONST.snr_min=5;  %dB

%%% for imageFOL, set user parameters
%%%  parameters are: [velD_change max_vel snr_min];  with vels defined in cm/s
%%%  see imageFOL for more details
CONST.imageFOL_user_param=[20 100 CONST.snr_min];  %for a 25MHz site
%CONST.imageFOL_user_param=[40 300 CONST.snr_min];  %for a 5MHz site, viewing the Gulf Stream
%CONST.imageFOL_user_param=[30 100 CONST.snr_min];  %for a 13MHz site

%%% define radave RM thresholds... to be used to weed out bad radials before the radial average is made 
%CONST.radave_thresholds=[snr_thresh angpeak_thresh angwidth_thresh(1) angwidth_thresh(2)];
CONST.radave_thresholds=[CONST.snr_min 5 0 50];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set processing choices

%%% Which first order line method to use
CONST.FOL_type='imageFOL';    %uses Image FOL methods
%CONST.FOL_type='codarFOL';    %imports the COS FOLs from the spectral file

%%% Which radial averaging method to output, both are estimated in the
%%% radave step, but only one can be sent out within the lluv file format.
%CONST.radave_type='regular';    %follows COS arthmetic averaging of the radials
                                %  within the CONST.bearing_width area and error calcuation
CONST.radave_type='metricQC';   %follows Kirincich et al 2012 to use thresholding
                                %  and power-weighted spatial means for rad ave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set plotting and printing switches
%plotting:
CONST.goplot=[1 0];   %where:
%  switch 1 controls the final plots of radial results, etc.  (1=plot, 0=do not plot)
%  switch 2 controls the intermediary plots within functions  (1=plot, 0=do not plot)
CONST.goprint=[1 0];   
%  switch 1 prints the final plots of radial results, etc.  (1=print, 0=do not print)
%  switch 2 prints the intermediary plots within functions  (1=print, 0=do not print)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set how 'new' files to process are found

%%% Method 1: choose files based on the file time
%CONST.files_to_process_method=1;
%%% Method 2: compare css files to RM files, and process new files
CONST.files_to_process_method=2;

%%% set the time frame to examine, only used by method 1
start_time=datenum(2011,1,1,0,0,0);
end_time=datenum(2012,1,1,0,0,0);
CONST.files_to_process_dates=[start_time end_time];


%%%%%%%%%%% end radial processing setup %%%%%%%%%%    
%%    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running the processing scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% setup logfile name for diary
logfn=['LOG_v7_' mfilename '_'  datestr(now,30) '.log'];
diary([log_dir '/' logfn]);
tic

%%% for this site, prepare for processing
HFR_DP_spectra2radial_prepwork_v3  
%%% sets up processing directories, loads site info and cal, and 
%%% identifies files that need processing


%  %limit the number of files  on in any 1 pass of this processing.  
%  % this can be helpful if running in realtime
%  if length(fnames)>5;
%      fnames=fnames(1:5);
%  end


%%
%%% now loop over each file in fnames to process data
for jjj=1:length(fnames)
%%
%%% display the file to be processed
disp(fnames(jjj))

%%%% start processing steps list now that HEAD has been established.
HEAD.ProcessingSteps={mfilename};
    
%%%% load data and run radial metric processing 
HFR_DP_spectra2radialmetric_process_v5  

%again, the spectra2radialmetric script returns a structure RM where the 
% RM.data array is same as COS RSv7 radial metric files with a few small 
% changes to units and field content:
%
% cols.   fields
% 1-2    lat lon
% 3-4    u v   (here nan as will not be used)
% 5      flag    (here nan, used by COS but not here)
% 6      range
% 7      bearing
% 8      rad vel
% 9      direction
% 10     rangecell
% 11     dopcell
% 12     angselect
% 13-15  musicSnglang musicDuel1ang musicDuel2ang
% 16-18  musicEigenRatio musicpowerRatio musicoffRatio
% 19-21  musicSnglpow(v) musicDuel1pow(v) musicDuel2pow(v) (in voltage (v) not power (db) )
% 22-24  musicSngl_pkwidth musicDuel1_pkwidth musicDuel2_pkwidth
% 25-27  musicDOASnglpeak musicDOADuel1peak musicDOADuel2peak
% 28-30  spectraA1(snr) spectraA2(snr) spectraA3(snr)
% 31-33  Eigenval1 Eigenval2 Eigenval3
% 34     Dual_reject  (still different, see the MUSIC function for details )
        
        
%%% save this RM ouput...you don't have choice %%%
%move to saving directory
eval(['cd ' outgoing_radialmat_file_dir])
f=char(fnames(jjj));
d=datevec(fdates(jjj));
i=find(f=='_'); i2=find(f=='.');
s=['RM_' site_name '_' num2str(d(1)) f(i(3):i2(1)) 'mat'];
eval(['save ' s ' RM;'])
%move to inital base directory
eval(['cd ' base_dir])  

%%  
%%% plot the RM result on a map of the array domain %%%
if CONST.goplot(1)==1
    figure(3); clf;
    HFR_DP_quickmap(ARRAY,working_dir);
    title(['RM for Site:' RM.Site_name ', Time:' datestr(RM.time)])
    %m_plot(RM.data(:,1),RM.data(:,2),'k.')
    m_quiver(RM.data(:,1),RM.data(:,2),RM.data(:,3),RM.data(:,4),10,'r');
    
    %%% if saving figures.
    if CONST.goprint(1)==1;
        eval(['cd ' outgoing_pics_file_dir])
        
        %save pictures of spectra
        figure(7); hold on;
        i=find(s=='.');   s2=s(1:i(1)-1);
        t=text(-100,5,[s2 '      raw spectra w/ridge lines']); set(t,'interpreter','none','rotation',90);
        set(gcf,'paperposition',[.25 .25 10 6])
        print('-djpeg','-r300',[s2 '_sfols_' date '.jpg'])
        
        %save pictures of radials
        figure(3); hold on;
        t=title([s '      RM locations']); set(t,'interpreter','none');
        set(gcf,'paperposition',[.25 .25 5 5])
        print('-djpeg','-r300',[s2 '_rmlocs_' date '.jpg'])
        %move back to inital base directory
        eval(['cd ' working_dir])
    end% if go print
end% if goplot
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Radial Average (equivalent to C.O.S. radial short result)

%%%% if RM is big enough continue to use RM to process radial metrics to radave files.
if length(RM.data)>500;
    
%%% compute radial averages   
[RADIAL,HEAD]=HFR_DP_RadialAverage_v5(RM,patt,HEAD,CONST);

%%% Clean large radials using ImageFOL user parameter for max velocity found
for k = 1:numel(RADIAL)
    [Rclean,I]=cleanRadials_v7(RADIAL(k),CONST.imageFOL_user_param(2));
end

%%% save this Rclean ouput...you don't have choice %%%
%move to saving directory
eval(['cd ' outgoing_radialavemat_file_dir])
f=char(fnames(jjj));
d=datevec(fdates(jjj));
i=find(f=='_'); i2=find(f=='.');
s=['RadAve_' site_name '_' num2str(d(1)) f(i(3):i2(1)) 'mat'];
eval(['save ' s ' Rclean;'])
%move to inital base directory
eval(['cd ' working_dir])  
      
%%% plot the result on a map of the array domain %%%
if CONST.goplot(1)==1
    figure(4); clf;
    HFR_DP_quickmap(ARRAY,working_dir);
    title(['Radave for Site:' RM.Site_name ', Time:' datestr(RM.time)])
    i=find(isnan(RADIAL.U(:,2))==0);
    m_quiver(RADIAL.LonLat(i,1),RADIAL.LonLat(i,2),RADIAL.U(i,2),RADIAL.V(i,2),10,'r')
    
    %%% if saving figures.
    if CONST.goprint(1)==1;
        eval(['cd ' outgoing_pics_file_dir])
        
        %save pictures of radials
        figure(4); hold on;
        t=title([s '      radave locations']); set(t,'interpreter','none');
        set(gcf,'paperposition',[.25 .25 5 5])
        i=find(s=='.');
        print('-djpeg','-r300',[s(1:i(1)-1) '_radave_' date '.jpg'])
        %move to inital base directory
        eval(['cd ' working_dir])
    end % if goprint
end% if goplot
%% 

%%%% convert existing, cleaned Rclean file to an ascii output in the
%%%% COS lluv format,

%move to saving directory
eval(['cd ' outgoing_radialavelluv_file_dir])

[Rclean,HEAD]=HFR_DP_RadAve2LLUV_v5(Rclean,CONST,PATT,HEAD,CSS_Head);

%move back to inital base directory
eval(['cd ' working_dir])  


end %if RM is long enough to compute radial averages...

end  % end jjj over fnames files


toc
diary off;
clear all; close all;

%%
return

