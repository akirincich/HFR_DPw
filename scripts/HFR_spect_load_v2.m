function [H,D]=HFR_spect_load_v2(fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [H,D]=HFR_spect_load_v2(fname)
%
%  Function reads in a CSS or a CSQ file from the COS binary format into matlab.
%  following COS documentation.  This has been tested on a number of file
%  versions (up to version 4, I think) with mostly success.  Kudos to Bill 
%  Rector for keeping the file formats normal and well documented
%  throughout the many interations.  Go Bill!
%
%   Versions
%  v2,  04/28/2014
%    added capability to not load D.D, the quality indicator, if a CSQ
%       file
%
%   3/2017   polished to package version,
%
%  Anthony Kirincich
%  WHOI PO
%  akirincich@whoi.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%% if as a function
fid=fopen([fname],'r','ieee-be');

%%%% get header information
H=[];
H.nCsaFileVersion = fread(fid,1,'int16');
H.nDateTime = datenum(1904,1,1)+fread(fid,1,'uint32')./3600/24;  %now in matlab time
H.nV1Extent = fread(fid,1,'int32');
H.nCsKind = fread(fid,1,'int16');
H.nV2Extent = fread(fid,1,'int32');
H.nSiteCodeName = char(fread(fid,4,'char'))';
H.nV3Extent = fread(fid,1,'int32');
H.nCoverageMinutes = fread(fid,1,'int32');

H.bDeletedSource = fread(fid,1,'int32');
H.bOverrideSourceInfo = fread(fid,1,'int32');

H.fStartFreqMHz = fread(fid,1,'float');
H.fRepFreqHz = fread(fid,1,'float');
H.fBandwidthKHz = fread(fid,1,'float');
H.bSweepUp = fread(fid,1,'int32');

H.nDopplerCells = fread(fid,1,'int32');
H.nRangeCells = fread(fid,1,'int32');
H.nFirstRangecell = fread(fid,1,'int32');

H.fRangeCellDistKm = fread(fid,1,'float');
H.nV4Extent = fread(fid,1,'int32');

%%% this portion is needed if running a version 5 file,
%%% would be noted above with H.nCsaFileVersion and 
%%%  H.nV4Extent would ~=0
%%%  but has not been well tested.
%
% H.nOutputInterval = fread(fid,1,'int32');
% H.nCreatorTypeCode = fread(fid,4,'char');
% H.nCreatorVerison = fread(fid,4,'char');
% H.nActiveChannels = fread(fid,1,'int32');
% H.nSpectraChannels = fread(fid,1,'int32');
% H.nActiveChannelBits = fread(fid,1,'int32');
% H.nV5Extent = fread(fid,1,'int32');
%%%


%make arrays for ant amp and cross amps
a=nan.*ones(H.nRangeCells,H.nDopplerCells);
D.a1=a; D.a2=a; D.a3=a;
D.a12=a+sqrt(-1).*a; D.a13=D.a12; D.a23=D.a12;
D.q=a;

%%

for jj=1:H.nRangeCells

D.a1(jj,:) = fread(fid,H.nDopplerCells,'float');
D.a2(jj,:) = fread(fid,H.nDopplerCells,'float');
D.a3(jj,:) = fread(fid,H.nDopplerCells,'float');

a=reshape(fread(fid,2*H.nDopplerCells,'float'),2,H.nDopplerCells);
D.a12(jj,:)=a(1,:)+a(2,:).*sqrt(-1);  
a=reshape(fread(fid,2*H.nDopplerCells,'float'),2,H.nDopplerCells);
D.a13(jj,:)=a(1,:)+a(2,:).*sqrt(-1);  
a=reshape(fread(fid,2*H.nDopplerCells,'float'),2,H.nDopplerCells);
D.a23(jj,:)=a(1,:)+a(2,:).*sqrt(-1);  

%add in quality data, if css file,   Quality factor (0-1) is a function of
%the state of the CSpro averaging process.  High if most of the data from the
% averaging period is being used, and low if only the current (end of
%averaging period) data is being used (see CSPro.pdf for details..
    if H.nCsKind==2;
        D.q(jj,:) = fread(fid,H.nDopplerCells,'float');
    end
end


fclose(fid);
