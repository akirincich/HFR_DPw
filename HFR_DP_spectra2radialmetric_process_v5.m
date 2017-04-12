%HFR_DP_spectra2radialmetric_process_v?.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_spectra2radialmetric_process_v?.m
%
%  This script runs the bulk of the processing steps used to derive radial 
%  velocity estimates from the spectra observed from CODAR-type HF radars
%  
%  In short, this script performs the following operations using the
%  user-defined choices and constants described by HFR_DP_master_XXXX.m,
%  including:
%
%  -Set up the outgoing radial metrics file for this spectral estimate
%  -Setup the file-spc
%  -Load the spectra file, and form the file-specific meta data critical
%      to estimating radial velocities.
%  -Establish or load the first order limits (FOLs)
%  -Use MUSIC to perform direction finding on the FO data following the CODAR-
%      type processing framework defined and described by Lipa et al 2006.
%  -Finallize the output MATLAB Structure, that contains the 'radial
%      metric' results 
%
%
% INPUT:  As this is a script, not a function, everything in the workspace
%           is passed to the script for its use
%
% OUTPUT: As this is a script, not a function, everything is passed back,
%          but the key addition is:
%
%   RM  --  A matlab structure that matrix of radial metric output,
%           including the following fields:
% 
% RM.Site_name            The 4 character site name
% RM.Site_loc             [Lon Lat] of the radar
% RM.fname                The file name processed
% RM.time                 Matlab time of the observed spectra/radials
% RM.rang                 Ranges (km)
% RM.c_vel                Doppler velocity of the spectral estimate (m/s)           
% RM.Stats                the fraction of single and duel angle solutions reported
% RM.data                 *** See Below ***
% RM.FOregi               the spectral indices that are processed
% RM.ProcessingSteps      A cell array list of all the scripts/functions
%                           used to get to this point in the processing
% RM.patt_UUID            The unique measured pattern identifier taken from
%                           the meas_patt.txt file for this radar
%
%
% Output format of RM.data is same as COS RSv7 radial metric files
% with a few small changes to units and field contents:
% cols   fields
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
%
%
% Version: 
% v5        March 2017
%
% Anthony Kirincich
% Woods Hole Oceanographic Institution
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%% mark that we are using this file to process...
HEAD.ProcessingSteps{end+1}=mfilename;

    %create output structure
    RM=[];
    RM.Site_name=HEAD.Site_name;
    RM.Site_loc=HEAD.Site_loc;
    RM.fname=[];
    RM.time=[];
    RM.rang=[];
    RM.c_vel=[];
    RM.Stats=[];
    RM.data=nan.*ones(1,34);
    RM.FOregi=[];
    
%%%% load CSS file
         [CSS_Head,Data]=HFR_spect_load_v2([incoming_spectra_file_dir '/' char(fnames(jjj))]);
        
        %%%% calculated constants from the CSS file %%%
        SpecHead=[];
        SpecHead.rad_res=CSS_Head.fRangeCellDistKm;
        SpecHead.Fc=CSS_Head.fStartFreqMHz.*1e6+CSS_Head.fBandwidthKHz.*1e3/2;
        %establish frequency, period
        SpecHead.Tr=1/CSS_Head.fRepFreqHz; %s
        SpecHead.fr=CSS_Head.fRepFreqHz;  %Hz
        SpecHead.fDmax=.5*SpecHead.fr;
        
        %array size
        [n m]=size(Data.a3);
        SpecHead.rangecell=1:n;  SpecHead.rang=SpecHead.rangecell.*SpecHead.rad_res;
        SpecHead.spectra_length=m;
        
        % more frequencies
        SpecHead.delta_f=2*SpecHead.fDmax/(SpecHead.spectra_length);     % correct .....
        SpecHead.doppler_freq=(-SpecHead.fDmax+SpecHead.delta_f:SpecHead.delta_f:SpecHead.fDmax);  %correct......
        
        %%% find what the bragg frequency is
        SpecHead.L=(CONST.c/SpecHead.Fc)/2;
        SpecHead.v_p=sqrt(9.81*SpecHead.L/2/pi);
        SpecHead.FBragg=2*SpecHead.v_p*SpecHead.Fc/CONST.c;
        [s,i]=sort(abs(abs(SpecHead.doppler_freq)-SpecHead.FBragg));
        SpecHead.iFBragg=i(1:2);
        
        %get doppler_vel
        SpecHead.doppler_vel=SpecHead.doppler_freq*CONST.c/SpecHead.Fc/2;
        %c_vel=doppler_vel-[-v_p*ones(1,(spectra_length-2)/2) 0 v_p*ones(1,(spectra_length-2)/2)];
        SpecHead.c_vel=SpecHead.doppler_vel-[-SpecHead.v_p*ones(1,(SpecHead.spectra_length)/2) SpecHead.v_p*ones(1,(SpecHead.spectra_length)/2)];

        %%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %occasionally, there will be an Inf in the data, not sure why. 
        % fix by changing all indices with data.a3==Inf to 1e-17
        % works out to be a amp of -204 dB
         i=find(Data.a3==Inf);
        %  set all of bad range bin to background for this bin
        %  and the bin before it
        Data.a3(i)=nan; Data.a1(i)=nan; Data.a2(i)=nan;
        Data.a12(i)=nan; Data.a13(i)=nan; Data.a23(i)=nan;
     
        j=find(isnan(mean(Data.a3,2))==1);
        i=find(abs(SpecHead.c_vel)<0.75);
        for jj=1:length(j)
            a=j(jj)-1:j(jj)+1;
            if a(end)>n; a(end)=n; end               
            Data.a1(j(jj),:)=(nanmean(Data.a1(a,:)));
            Data.a2(j(jj),:)=(nanmean(Data.a2(a,:)));
            Data.a3(j(jj),:)=(nanmean(Data.a3(a,:)));
            Data.a23(j(jj),:)=(nanmean(Data.a23(a,:)));
            Data.a13(j(jj),:)=(nanmean(Data.a13(a,:)));
            Data.a12(j(jj),:)=(nanmean(Data.a12(a,:)));
        %now change all data from bad rows inside bragg region to min
        %values
       Data.a3(j(jj),i)=1e-17; Data.a1(j(jj),i)=1e-17; Data.a2(j(jj),i)=1e-17;
         Data.a12(j(jj),i)=1e-17; Data.a13(j(jj),i)=1e-17; Data.a23(j(jj),i)=1e-17;   
        end
 %%               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% identify the first order region %%
FOreg=[]; FOregi=[];

if strcmp(CONST.FOL_type,'codarFOL')==1;  
%%%  Use COS estimated FOLs (aka 'ALims') available in the resource fork
%%%   of the css file AFTER THE FILE AS BEEN PROCESSED BY CSPRO.
%%%    Note that if the CSS has been 'improperly transferred', where by
%%%    the resource fork of the file has been stripped out, this will
%%%    not be available and the script will error out and you are stuck

%%% two ways to load the COS FOLS exist...
%%%  (1) call a perl script to strip and save as a text file, then load and
%%%   process
%[FOreg,FOregi,HEAD]=HFR_spectrsrc_load_v2([incoming_spectra_file_dir '/' char(fnames(jjj))],scripts_dir,CSS_Head,HEAD);

%%% (2)  use a shell command to copy the resource fork to a binary file, and
%%%   then read in following the directions of the perl script and COS, but
%%%   directly into matlab
%[FOreg,FOregi,Alims,HEAD]=HFR_spectrsrc_load_v3([incoming_spectra_file_dir '/' char(fnames(jjj))],scripts_dir,CSS_Head,HEAD);

[FOreg,FOregi,Alims,HEAD]=HFR_spectrsrc_load_v4([incoming_spectra_file_dir '/' char(fnames(jjj))],scripts_dir,CSS_Head,HEAD,SpecHead);

%%
%%% plot the result if you wish
if CONST.goplot(1)==1;
    gain3=10*log10(abs(Data.a3)) + (-40. + 5.8);  %only works with plus...typo in COS manual
    [n m]=size(gain3);

    figure(7); clf;
    subplot(211); pcolor(gain3); shading flat; colorbar; caxis([-160 -80]); title('Raw Ant3 Spectral Power with segments (red) and sFOLs (white) shown')
    hold on;
    [cs,ho]=contour(FOreg,[1 1],'w');
    plot(Alims,(1:n)'*[1 1 1 1],'k','linewidth',2)
end



elseif strcmp(CONST.FOL_type,'imageFOL')==1;    
%%% uses Image FOL methods (see Kirincich JOAT, 2016) for details, or the 
%%%  script to estimate the dynamic FOLs

%%% set last constant
CONST.v_incr=median(diff(SpecHead.c_vel));
%%% run imageFOLS
[FOreg,FOregi,Alims,HEAD,DN_out]=imageFOLs_v4(Data.a3,SpecHead.iFBragg,CONST.v_incr,CONST.imageFOL_user_param,CONST.goplot,HEAD);
                
end

if isempty(FOreg)==1;
    disp(['error...incorrect FOL_type: ' char(CONST.FOL_type) ' not identified'])
    asdfasdfsdf
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% process FO data with music if the file has enough FO data to be
%%% interesting
if length(FOregi)>100
    tic
    [R, HEAD]=HFR_DPmusic_v7(Data,FOregi,HEAD,SpecHead,patt,CONST);
    toc
else
    R=nan.*ones(1,29);
end
        
%%
%only move forward if have a significant number of returns...i.e.
%length(R) > 500
if length(R)>500;
    
    %%% Finish radial metric file creation,
    %%%  calc lon lat u v flag for each line
    r=nan*ones(length(R),5);
    
    x=R(:,1)*[1 1].*[sind(R(:,2)) cosd(R(:,2))];   %reverse sin cos to account for true vs math coord
    [lon,lat]=km2lonlat([HEAD.Site_loc(2).*ones(length(R),1)],[HEAD.Site_loc(1).*ones(length(R),1)],x(:,1),x(:,2));
    r(:,1:2)=[lon lat];
    %find U,V
    r(:,3:4) = [R(:,3).*cosd(true2math(R(:,4)))    R(:,3).*sind(true2math(R(:,4)))];
    
    %         %%%% test range, bearing estimates of x,y in km
    %         d=[45 135 225 315];
    %         dm=true2math(d)
    %         x=R(1,1).*[sind(d)' cosd(d)']  %reverse sin cos to account for true vs math coord
    %         xm=R(1,1).*[cosd(dm)' sind(dm)']
    
    %%
    
    %place results into metrics file
    RM.fname=char(fnames(jjj));
    RM.time=CSS_Head.nDateTime;
    RM.rang=SpecHead.rang;
    RM.c_vel=SpecHead.c_vel;
    RM.FOregi=FOregi;
    i=find(R(:,7)==1);
    
    % stats on the fraction of single and duel angle solutions being
    % reported
    RM.Stats=[length(i)./length(R) (length(R)-length(i))./length(R)];
    RM.data=[r R];
    
    RM(1).ProcessingSteps=HEAD.ProcessingSteps;
    RM.patt_UUID=patt.UUID;
    
    %again, output format of RM(jjj).data is same as COS RSv7 radial metric files
    %with a few small changes to units and field contentx
    % 1-2 lat lon
    % 3-4 u v   (here nan as will not be used)
    % 5 flag    (here nan, used by COS but not here)
    % 6 range
    % 7 bearing
    % 8 vel
    % 9 direction
    % 10     rangecell
    % 11     dopcell
    % 12     angselect
    % 13-15   musicSnglang musicDuel1ang musicDuel2ang
    % 16-18   musicEigenRatio musicpowerRatio musicoffRatio
    % 19-21 musicSnglpow(v) musicDuel1pow(v) musicDuel2pow(v)      (in voltage (v) not power (db) )
    % 22-24 musicSngl_pkwidth musicDuel1_pkwidth musicDuel2_pkwidth
    % 25-27 musicDOASnglpeak musicDOADuel1peak musicDOADuel2peak
    % 28-30 spectraA1(snr) spectraA2(snr) spectraA3(snr)
    % 31-33 Eigenval1 Eigenval2 Eigenval3
    % 34    Dual_reject  (still a... )
    
else  %if there wasn't enough data to make a real file
    RM(1).fname=char(fnames(jjj));
    RM(1).time=CSS_Head.nDateTime;
    RM(1).rang=SpecHead.rang;
    RM(1).c_vel=SpecHead.c_vel;
    %   i=find(R(:,7)==1);
    %   RM(jjj).Stats=[length(i)./length(R) (length(R)-length(i))./length(R)];
    %   RM(jjj).data=[r R];
    RM(1).Stats=nan;
    RM(1).data=nan.*ones(1,34);
    RM(1).ProcessingSteps=HEAD.ProcessingSteps;
    RM(1).patt_UUID=patt.UUID;
    
    
end  %length(R) >500

%%% go back main script to save RM and possibly the figures


