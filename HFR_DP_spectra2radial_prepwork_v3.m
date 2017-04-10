%HFR_DP_spectra2radial_prepwork_v?
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HFR_DP_spectra2radial_prepwork_v?
%
% this script sets up the proper working directories
% and file directories to process the site's spectra
% files and isolates the new files in the css folder that 
% need processing.
%
% version:
%  04/2/2014 updated to include for times of 0 files in RM folder.
%
%  7/11/2016  to version 2, fix issues
%
%  3/14/2017 only look for files starting with 'CSS' or 'RM'
%            transition to package with terminology edits
%
% Anthony Kirincich
% WHOI PO
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set the directories of files
%processing_file_dir=['/Codar/SeaSonde/Data/Processings/HFR_DP/Site_' site_name];
incoming_spectra_file_dir=[base_dir '/Site_' site_name '_css'];
incoming_config_file_dir=[base_dir '/Site_' site_name '_config'];
outgoing_radialmat_file_dir=[base_dir '/Site_' site_name ''];
outgoing_radialavemat_file_dir=[base_dir '/Site_' site_name '_radave'];
outgoing_radialavelluv_file_dir=[base_dir '/Site_' site_name '_radave_lluv'];
outgoing_pics_file_dir=[base_dir '/Site_' site_name '_pics'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load header and patt files

%load header file
s=[incoming_config_file_dir '/RadialConfigs/Header.txt'];
HEAD=open_header_v05012014(s);

%load meas pattern file
s=[incoming_config_file_dir '/RadialConfigs/MeasPattern.txt'];
[PATT]=open_measpatt_v04282014(s);

%%% using measured patterns, unrap data into usable form for music
patt=[];
%bring along pattern info
patt.Site_name=PATT.Site_name; patt.Site_loc=PATT.Site_loc;
patt.Patt_date=PATT.Patt_date; patt.UUID=PATT.UUID;
patt.angles=PATT.mpatt(1,:);
%%% The pattern bearings are CCW (counter-clockwise) degrees referenced 
%%% from the antenna bearing. The antenna bearing is found in Header.txt
%%% and is (CW) clockwise degrees from true North. See the File_RadialSetups guide.
%%% for details.
%adjust angles to be in 'true'
patt.angles=PATT.Antbear-patt.angles;
i=find(patt.angles>360); patt.angles(i)=patt.angles(i)-360;
i=find(patt.angles<0); patt.angles(i)=patt.angles(i)+360;

%Define antenna responses at those angles
patt.a1=PATT.mpatt(2,:)+PATT.mpatt(4,:).*sqrt(-1);
patt.a2=PATT.mpatt(6,:)+PATT.mpatt(8,:).*sqrt(-1);
patt.a3=ones(size(patt.angles));

%%
% make bearing limits for radave files...
% this is harder than it looks to do it a systematic way for all types of patterns
%  WHY? patt.angles starts on one side of the pattern and progresses to the other
%  regardless if there is a discontinuity in the true angles (i.e. crossing 0/360 
%  with jumps).  
% 
% The goal is to make a radave angle array that matches the 
%  direction and sense of rotation of the bearing as patt.angle, but has 
% standardized angles.  This will limit the radial concat issues one might have
% later if the pattern or Antbear is changed.

% %%%%%% old way...this doesn't work if the site wraps around True North
% %%%%%% but is not full coverage
% bearing_limits=[min(patt.angles)+2 max(patt.angles)-2];
% %establish 5deg bearing angle centers and ranges
% BearT=[bearing_limits(1):CONST.bearing_width:bearing_limits(2)];
% BearT_ends=[BearT-CONST.bearing_width/2 BearT(end)+CONST.bearing_width/2];

%%%%%% alternative...define as having centers such that no width spans 0/360
stock_Bear=CONST.bearing_width/2:CONST.bearing_width:360-CONST.bearing_width/2;
db=abs(nanmedian(diff(patt.angles)));
BearT=nan*ones(size(patt.angles));
for ii=1:length(patt.angles)
    [s,i]=sort(abs(stock_Bear-patt.angles(ii))); %find the nearest ave bearing
    if s(1)<= CONST.bearing_width/2 ; %if it is less than CONST.bearing_width/2 away
        BearT(ii)=stock_Bear(i(1));
    end
end
%%% now march through BearT and save only the first value, if there are
%%% more values than CONST.bearing_width/2, so, 3 or more for a 
u=unique(BearT);  
% can't just use unique(BearT,'stable') as also wish to kick out the ones
% that have less than CONST.bearing_width/db/2 values.
for ii=1:length(u)
i=find(BearT==u(ii));
if length(i)>CONST.bearing_width/db/2 % /db is included to account for APMs that don't have increments of 1
    BearT(i(2:end))=nan;
else %cut them all
    BearT(i(1:end))=nan;
end
end
%%% condense
BearT=BearT(~isnan(BearT));
%%% make ends.
BearT_ends=[BearT-CONST.bearing_width/2; BearT+CONST.bearing_width/2];
%%% fix to correct the 0 end if have a Bear? on the high side.
i=find(min(BearT_ends)==0 & [diff(BearT_ends)> CONST.bearing_width] );
if isempty(i)==0 %need to change a 0 to 360
    j=find(BearT_ends(:,i)==0); BearT_ends(j,i)=360;
end
%%% convert to math
Bear=true2math(BearT);
Bear_ends=true2math(BearT_ends);
i=find(min(Bear_ends)==0 & [abs(diff(Bear_ends))> CONST.bearing_width] );
if isempty(i)==0 %need to change a 0 to 360
    j=find(Bear_ends(:,i)==0); Bear_ends(j,i)=360;
end

patt.Bear=Bear;
patt.Bear_ends=Bear_ends;
patt.BearT=BearT;
patt.BearT_ends=BearT_ends;

%%% make a plot of the pattern, if wanted
if CONST.goplot(2)==1
    figure(2); clf
    p=polar(true2math(patt.angles)*pi./180,abs(patt.a1),'r'); hg
    set(p,'linewidth',2)
    p=polar(true2math(patt.angles)*pi./180,abs(patt.a2),'b');
    set(p,'linewidth',2)
    p=polar([true2math(HEAD.Bearing)*pi/180]*[1 1],[0 1],'k'); set(p,'linewidth',2);
    title([patt.Site_name '  Meas pattern (converted from DegT to math)']);
    %%% add the azimuths of the radave product, for show
    p=polar(Bear*pi./180,ones(size(Bear)),'kx');
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get list of spectra in css folder
f=dir([incoming_spectra_file_dir '/CSS*' site_name '*cs4']);
fnames={}; fdates=[];
for i=1:length(f)
    fnames(i)={f(i).name};
    aa=char(fnames(i));
    j=find(aa=='_');
    %year is between 2 and 3, month 3 and 4, day 4 and 5, time 5 and 6
    fdates(i)=datenum(2000+str2num(aa(j(2)+1:j(2)+2)),str2num(aa(j(3)+1:j(3)+2)),str2num(aa(j(4)+1:j(4)+2)),str2num(aa(j(5)+1:j(5)+2)),str2num(aa(j(5)+3:j(5)+4)),0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get list of RM files in RM  folder
f=dir([outgoing_radialmat_file_dir '/RM*' site_name '*mat']);
fmnames={}; fmdates=[];
for i=1:length(f)
    fmnames(i)={f(i).name};
    aa=char(fmnames(i));
    j=find(aa=='_');
    %year is between 2 and 3, month 3 and 4, day 4 and 5, time 5 and 6
    fmdates(i)=datenum(str2num(aa(j(2)+1:j(2)+4)),str2num(aa(j(3)+1:j(3)+2)),str2num(aa(j(4)+1:j(4)+2)),str2num(aa(j(5)+1:j(5)+2)),str2num(aa(j(5)+3:j(5)+4)),0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% identify which of the files in fnames to process

%%% Method 1: choose files based on the file time
if CONST.files_to_process_method==1;
    iwhich=find(fdates>=CONST.files_to_process_dates(1) & fdates<CONST.files_to_process_dates(2) );
    
%%% Method 2: compare css files to RM files, and process new files
elseif CONST.files_to_process_method==2;    
    if isempty(fmdates)==1
        iwhich=1:length(fdates);
    elseif isempty(fmdates)==0
        iwhich=find(fdates>fmdates(end));
    end
end


%%% return the files to process to running script in the fname and fdate
%%% arrays
fnames=fnames(iwhich);
fdates=fdates(iwhich);


