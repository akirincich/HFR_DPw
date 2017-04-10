function  [Rclean,HEAD]=HFR_DP_RadAve2LLUV_v5(Rclean,CONST,PATT,HEAD,CSS_Head);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_radave2lluv.m
%
% Script uses PATT, CSS_Head, HEAD, and Rclean to fill out the required 
%  elements of the COS lluv file format, writes a new file to the current 
%  directory with the outfilename:
%
%   ['RDLm' '_' char(Rclean.SiteName) '_' tt '.ruv'];
%
%  where tt is the time of the data
%
%  Version:
%  -v3  changed 
%    fprintf(fid,'%%RangeResolutionMeters: %1.7f\n',CSS_Head.fRangeCellDistKm);
%  to 
%    fprintf(fid,'%%RangeResolutionKMeters: %1.7f\n',CSS_Head.fRangeCellDistKm);
%  for consistence with HFRnet requirements
%
%  -v4 fixed transmit center frequency to allow sweep direction to 
%
%  -v5  3/2016
%       fixed to correct Range Resolution tag
%       added the spatial error estimate
%       prepped for package, functionized
%
%
%  Anthony Kirincich
%  WHOI-PO
%  akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
Rclean.ProcessingSteps{end+1}=mfilename;

%new constants for output to lluv
PROC_CONST=[];
PROC_CONST.MergedCount=1;
PROC_CONST.RadialMinimumMergePoints=1;
PROC_CONST.Temp_qual=999;
PROC_CONST.VelocityMaximum=999;
PROC_CONST.VelocityMimimum=999 ;
PROC_CONST.TemporalCount=1;
PROC_CONST.COS_Flag=999;

if strcmp(CONST.radave_type,'regular')==1;
    PROC_CONST.iwhich_ave=1;
elseif strcmp(CONST.radave_type,'metricQC')==1;
    PROC_CONST.iwhich_ave=2;
end

%%% use Rclean.U to figure out which rows are good.
ikeep=find(isnan(Rclean.RadComp(:,PROC_CONST.iwhich_ave))==0);


% %%%%%%%% setup proper outputfile name  and directory %%%%%%%%%%%
%make proper file name 
f_type='RDLm';
%set proper time stamp    
tt=datestr(Rclean.TimeStamp,31); tt=tt(1:16);
i=find(tt=='-' | tt==' '); tt(i)='_';
i=find(tt~=':'); tt=tt(i);
outfilename=[f_type '_' char(Rclean.SiteName) '_' tt '.ruv'];

%%% lets go!
fid=fopen(outfilename,'w');

fprintf(fid,'%%CTF: 1.00\n%%FileType: LLUV rdls "RadialMap"\n%%LLUVSpec: 1.12 2009 11 05\n%%Manufacturer: CODAR Ocean Sensors. SeaSonde\n');
fprintf(fid,'%%Site: %s ""\n',char(Rclean.SiteName));
tt=datestr(Rclean.TimeStamp,31); i=find(tt=='-' | tt==':'); tt(i)=' ';
fprintf(fid,'%%TimeStamp: %s\n',tt);

%fprintf(fid,'%%TimeZone: %s\n', HEAD.Time_zone);  %it is always GMT.
i=find(HEAD.Time_zone==' ');
tzone=HEAD.Time_zone(1:i(1)-1);
toff=HEAD.Time_zone(i(1)+1:i(2)-1);
%there is never daylight savings time
fprintf(fid,'%%TimeZone: "%s" +%s 0\n',tzone,toff);  %it is always GMT.
fprintf(fid,'%%TimeCoverage: %2.4f Minutes\n',CSS_Head.nCoverageMinutes);
fprintf(fid,'%%Origin: %3.7f %3.7f\n', HEAD.Site_loc);
fprintf(fid,'%%GreatCircle: "WGS84" 6378137.000  298.257223562997''\n');
fprintf(fid,'%%GeodVersion: "CGEO" 1.57  2009 03 10\n');
fprintf(fid,'%%LLUVTrustData: all %%%% all lluv xyuv rbvd\n');
fprintf(fid,'%%RangeStart: 1\n');
fprintf(fid,'%%RangeEnd: %2.0f\n',CSS_Head.nRangeCells);

%%% COS format for the LLUV file is to report the results in meters if less
%%% then 500 m, but KMeters if over
if CSS_Head.fRangeCellDistKm < 0.5
    fprintf(fid,'%%RangeResolutionMeters: %1.7f\n',CSS_Head.fRangeCellDistKm.*1000);
elseif CSS_Head.fRangeCellDistKm > 0.5
    fprintf(fid,'%%RangeResolutionKMeters: %1.7f\n',CSS_Head.fRangeCellDistKm);
end

fprintf(fid,'%%AntennaBearing: %3.1f True\n',PATT.Antbear);
fprintf(fid,'%%ReferenceBearing: 0 True\n');
fprintf(fid,'%%AngularResolution: 5 Deg\n');
fprintf(fid,'%%SpatialResolution: 5 Deg\n');
fprintf(fid,'%%PatternType: Measured\n');

tt=datestr(PATT.Patt_date,31); i=find(tt=='-' | tt==':'); tt(i)=' ';
fprintf(fid,'%%PatternDate: %s\n',tt);

fprintf(fid,'%%PatternResolution: %1.1f deg\n',PATT.DegRes);
fprintf(fid,'%%PatternSmoothing: %2.1f deg\n',PATT.DegSmooth);
fprintf(fid,'%%PatternUUID: %s\n',PATT.UUID);
if CSS_Head.bSweepUp==1
      fprintf(fid,'%%TransmitCenterFreqMHz: %2.6f\n',CSS_Head.fStartFreqMHz + CSS_Head.fBandwidthKHz/1000/2);
else 
    fprintf(fid,'%%TransmitCenterFreqMHz: %2.6f\n',CSS_Head.fStartFreqMHz + -1.*CSS_Head.fBandwidthKHz/1000/2);
end
fprintf(fid,'%%DopplerResolutionHzPerBin: %1.9f\n',HEAD.ALT_freq_doppt_dopres(3));

%fprintf(fid,'%%BraggSmoothingPoints: %2.0f\n',PROC_CONST.BraggSmoothingPoints);

fprintf(fid,'%%CurrentVelocityLimit: %3.1f\n',CONST.imageFOL_user_param(2));

%fprintf(fid,'%%BraggHasSecondOrder: %2.2f\n',PROC_CONST.BraggHasSecondOrder);
%fprintf(fid,'%%RadialBraggPeakDropOff: %2.2f\n',PROC_CONST.RadialBraggPeakDropOff);
%fprintf(fid,'%%RadialBraggPeakNull: %2.2f\n',PROC_CONST.RadialBraggPeakNull);
%fprintf(fid,'%%RadialBraggNoiseThreshold: %2.2f\n',PROC_CONST.RadialBraggNoiseThreshold);

fprintf(fid,'%%PatternAmplitudeCorrections: %1.4f %1.4f\n',PATT.ampfactors);

%%% new versions of CrossLoopPatterner (in RSS v7 plus) put phase calculations into
%%% the pattern file, older versions do not...write if present.
if isfield(PATT,'PhaseCorrects')==1
    fprintf(fid,'%%PatternPhaseCalculations: %3.4f %3.4f\n',PATT.PhaseCorrects);
end

fprintf(fid,'%%PatternAmplitudeCalculations: %1.4f %1.4f\n',HEAD.Ampfacs_ideal1_2_meas1_2(3:4));
%fprintf(fid,'%%PatternPhaseCalculations: 999 999\n');
fprintf(fid,'%%RadialMusicParameters: %2.3f %2.3f %2.3f\n',HEAD.Musicparams123_gsmoothwidth_deg_smearwidth_velthresh(1:3));

fprintf(fid,'%%MergedCount: %2.0f\n',PROC_CONST.MergedCount);
fprintf(fid,'%%RadialMinimumMergePoints: %2.0f\n',PROC_CONST.RadialMinimumMergePoints);
fprintf(fid,'%%MergeMethod: 1 MedianVectors\n');

%specific to the HFR_DP methods.
fprintf(fid,'%%RadialSpatialAveragingMethod_HFR_DP: %s\n',CONST.radave_type);
fprintf(fid,'%%FirstOrderCalc_HFR_DP: %s\n',CONST.FOL_type);

fprintf(fid,'%%PatternMethod: 1 PatternVectors\n');
fprintf(fid,'%%TransmitSweepRateHz: %2.4f\n',CSS_Head.fRepFreqHz);
fprintf(fid,'%%TransmitBandwidthKHz: %3.6f\n',CSS_Head.fBandwidthKHz);
fprintf(fid,'%%SpectraRangeCells: %4.0f\n',CSS_Head.nRangeCells);
fprintf(fid,'%%SpectraDopplerCells: %4.0f\n',CSS_Head.nDopplerCells);
fprintf(fid,'%%TableType: LLUV RDL_HFR_DP_v1\n');   %started 3/2017

%%%% original version with all field following COS methods
%fprintf(fid,'%%TableColumns: 18\n');
%fprintf(fid,'%%TableColumnTypes: LOND LATD VELU VELV VFLG ESPC ETMP MAXV MINV ERSC ERTC XDST YDST RNGE BEAR VELO HEAD SPRC\n'); 

%%% cut columns 5 6 7 8 9 from original, as per Otero instructions
%fprintf(fid,'%%TableColumns: 13\n');
%fprintf(fid,'%%TableColumnTypes: LOND LATD VELU VELV ERSC ERTC XDST YDST RNGE BEAR VELO HEAD SPRC\n'); 

%%% added back column 6, the spatial error estimate, to conform to standards
fprintf(fid,'%%TableColumns: 14\n');
fprintf(fid,'%%TableColumnTypes: LOND LATD VELU VELV ESPC ERSC ERTC XDST YDST RNGE BEAR VELO HEAD SPRC\n'); 

%fprintf(fid,'%%TableRows: %5.0f\n',length(find(isnan(Rclean.RadComp(:,PROC_CONST.iwhich_ave))==0)));
fprintf(fid,'%%TableRows: %5.0f\n',length(ikeep));

fprintf(fid,'%%TableStart:\n');
%%% for 18 columns
%fprintf(fid,'%%%%   Longitude   Latitude    U comp   V comp  VectorFlag    Spatial    Temporal     Velocity    Velocity  Spatial  Temporal X Distance  Y Distance   Range   Bearing   Velocity  Direction   Spectra\n');
%fprintf(fid,'%%%%     (deg)       (deg)     (cm/s)   (cm/s)  (GridCode)    Quality     Quality     Maximum     Minimum    Count    Count      (km)        (km)       (km)    (True)    (cm/s)     (True)    RngCell\n');

%%% for 13 columns
% fprintf(fid,'%%%%   Longitude   Latitude    U comp   V comp   Spatial  Temporal X Distance  Y Distance   Range       Bearing   Velocity  Direction   Spectra\n');
% fprintf(fid,'%%%%     (deg)       (deg)     (cm/s)   (cm/s)    Count    Count      (km)        (km)       (km)       (True)    (cm/s)     (True)    RngCell\n');

%%% for 14 columns
fprintf(fid,'%%%%   Longitude   Latitude    U comp   V comp     Spatial Error  Spatial  Temporal X Distance  Y Distance   Range       Bearing   Velocity  Direction   Spectra\n');
fprintf(fid,'%%%%     (deg)       (deg)     (cm/s)   (cm/s)       (cm/s)        Count    Count      (km)        (km)       (km)       (True)    (cm/s)     (True)    RngCell\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%% file format
%1 Longitude (deg) 
%2 Latitude   (deg)
%3 U comp (cm/s)   
%4 V comp (cm/s)
% cut %5 VectorFlag   (GridCode)  
%  %6 Spatial Quality   
% cut %7 Temporal   Quality     
% cut %8 Velocity    Maximum   
% cut %9 Velocity    Minimum  
%10  Spatial Count   
%11 Temporal Count     
%12 X Distance  (km)     
%13 Y Distance (km)     
%14    Range  (km)   
%15    Bearing    (True)   
%16    Velocity   (cm/s)    
%17  Direction     (True)   
%18  Spectra RngCell

%%%
for j=1:length(ikeep);
    s=Rclean.RangeBearHead(ikeep(j),:);    
   R=[Rclean.LonLat(ikeep(j),:)'   ;
       Rclean.U(ikeep(j),PROC_CONST.iwhich_ave)   ;
       Rclean.V(ikeep(j),PROC_CONST.iwhich_ave)    ;
      % PROC_CONST.COS_Flag;
       %%%Rclean.Flag(ikeep(j),iwhich_ave)    ;
       Rclean.Error(ikeep(j),PROC_CONST.iwhich_ave)  ;
      % PROC_CONST.Temp_qual    ;
      % PROC_CONST.VelocityMaximum   ; 
      % PROC_CONST.VelocityMimimum   ; 
         
       Rclean.Flag(ikeep(j),2)-floor(Rclean.Flag(ikeep(j),2)./100)*100 ;

       PROC_CONST.TemporalCount    ;
%       s(1)*cosd(true2math(s(3)))   ;   
%       s(1)*sind(true2math(s(3)))   ; 
       s(1)*cosd(s(2))   ;   
       s(1)*sind(s(2))   ; 
   
       Rclean.RangeBearHead(ikeep(j),1)
       math2true(Rclean.RangeBearHead(ikeep(j),2))   ;
       Rclean.RadComp(ikeep(j),PROC_CONST.iwhich_ave)    ;
       math2true(Rclean.RangeBearHead(ikeep(j),3)   ) ;
       round(Rclean.RangeBearHead(ikeep(j),1)./(CSS_Head.fRangeCellDistKm))] ;
       
    if isnan(R(3))==0  %proceed only if good data exists in this entry
i=find(isnan(R)==1);
R(i)=999;

%%%with length(R)==18
%str='   %3.7f  %3.7f  %3.3f\t%3.3f\t %3.0f\t   %3.3f\t%3.3f\t   %3.3f\t %3.3f\t%2.0f\t%2.0f   %2.4f\t%2.4f  \t%2.4f\t %3.1f\t%2.3f\t   %3.1f  \t%2.0f\n' ; 

%%%with length(R)==13
%%%cut columns 5 6 7 8 9 from original
%str='   %3.7f  %3.7f  %3.3f\t%3.3f\t %2.0f\t%2.0f      %2.4f    %2.4f  \t%2.4f  \t %3.1f\t%2.3f\t   %3.1f  \t%2.0f\n' ; 

%with length(R)==14
%%% add back column 6 (spatial ave) to  length(R)==13 estimate
str='   %3.7f  %3.7f  %3.3f\t%3.3f\t %3.3f\t %2.0f\t%2.0f      %2.4f    %2.4f  \t%2.4f  \t %3.1f\t%2.3f\t   %3.1f  \t%2.0f\n' ; 

fprintf(fid,str,R);
%    pause
    end       
end
%%% done with data TABLE
fprintf(fid,'%%TableEnd:\n');
fprintf(fid,'%%%%\n');
t=datestr(now,30);
fprintf(fid,'%%ProcessedTimeStamp: %s %s %s  %s %s %s \n',t(1:4),t(5:6),t(7:8),t(10:11),t(12:13),t(14:15));

%%
%%% finish with the processing steps used...
for ii=1:length(Rclean.ProcessingSteps)
    s=Rclean.ProcessingSteps{ii};
    i=find(s=='_');
    if  isempty(i)==0 & s(i(end)+1)=='v'
        s1=s(1:i(end)-1);
        s2=s(i(end)+2:end);
        fprintf(fid,'%%ProcessingTool: "%s" %s.0\n', s1,s2 );
    else
    end
end
  %%%% this will output the results in a format like:
%ProcessingTool: "RadialMerger" 10.7.1
%ProcessingTool: "SpectraToRadial" 10.9.0
%ProcessingTool: "RadialSlider" 11.2.0
%ProcessingTool: "RadialArchiver" 11.2.2
%ProcessingTool: "AnalyzeSpectra" 10.7.4

%%% all done, get out.
fprintf(fid,'%%End:\n');
fclose(fid);

