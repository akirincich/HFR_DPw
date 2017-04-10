function HEAD=open_header(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function HEAD=open_header(filename);
%
%  function opens COS Seasonde Header.txt file for instrument in question
%  and loads all pertinant site info for computing radials from CSS file
%
%  returns structure HEAD with all info inside
%
%  Should be good for most header files for Radial suite R6 and R7
%
%   03/23/2012
%   Anthony Kirincich
%   WHOI
%   akirincich@whoi.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(filename);
d=textscan(fid,'%s%s','Delimiter','!');
d1=d{1};


%extract parameters of interest, place into structure HEAD
%HEAD=[];

%LINE 1  get site name
a=char(d1(1));
%first and second spaces bound the site name, but to be sure
% find the first 4 chars after the first continuous group of spaces
i=find(a==' ');   id=diff(i);
id1=find(id>1);
HEAD.Site_name=a(i(id1(1))+1:i(id1(1))+4);


%LINE 2   get site location...text order should not change
a=char(d1(2)); llsign=[0 0];
%clean text string
i=find(a=='¡' | a==',' | a==''''); a(i)=' ';
i=find(a=='S'); if isempty(i)==0; llsign(1)=-1; end
i=find(a=='W'); if isempty(i)==0; llsign(2)=-1; end
i=find(a=='N' | a=='S' | a=='E' | a=='W'); a(i)=' ';
HEAD.Site_loc=str2num(a);
HEAD.Site_loc=[HEAD.Site_loc(1)+HEAD.Site_loc(2)/60 HEAD.Site_loc(3)+HEAD.Site_loc(4)/60];
if llsign(1)==-1; HEAD.Site_loc(1)=-HEAD.Site_loc(1); end
if llsign(2)==-1; HEAD.Site_loc(2)=-HEAD.Site_loc(2); end


%LINE3      antenna bearing
HEAD.Bearing=str2num(char(d1(3)));

%LINE4       radial grid info 
HEAD.Radial_lastrangecell_1st_step=str2num(char(d1(4)));

%LINE7-9 ALT Freq, doppler specs
HEAD.ALT_freq_doppt_dopres=[str2num(char(d1(7))) str2num(char(d1(8))) str2num(char(d1(9)))];

%LINE10       time zone
HEAD.Time_zone=char(d1(10));

%LINE11     FOL max vel, numpts
HEAD.FOL_maxvel_numpts=str2num(char(d1(11)));

%LINE12    Factor Down Peak limit 1st order Radials , 0 = no 2nd order, Factor down 1rst order Waves
HEAD.FOL_Fdpeaklim_2ndord_Fdpeakwaves=str2num(char(d1(12)));

%LINE15    Radials:Factor down peak nulls, Noise factor, Waves:Factor down peak nulls, Noise factor'
HEAD.FOL_Fdpeaknulls_noisefac_Fdpeaknullswaves_noisewaves=str2num(char(d1(15)));

%LINE17    Ampl Factors Ideal:Amp1,Amp2  Patt:Amp1,Amp2'
HEAD.Ampfacs_ideal1_2_meas1_2=str2num(char(d1(17)));

%LINE18   coastline angles
HEAD.Coastangles=str2num(char(d1(18)));

%LINE19  Music params eigrat,sigprat,diagrat;  Gaussian Smooth Med Width deg, Smear Width, Vel Thresh
HEAD.Musicparams123_gsmoothwidth_deg_smearwidth_velthresh=str2num(char(d1(19)));

%LINE21   Rads:  Coverage min., Output Interval min., Interval Offset min., 0=Watch Timespan
HEAD.Rads_cov_out_offset_watch=str2num(char(d1(21)));

%LINE22   Bearing resolution
HEAD.Bearres=str2num(char(d1(22)));

%LINE24   Ionosphere Noise Removal Factor
HEAD.ionremove=str2num(char(d1(24)));

%LINE25  Doppler Cell Noise Limit Adjustment
HEAD.Dopnoiselim=str2num(char(d1(25)));

%LINE27   Radial First Range Cell, WaveModel First Range Cell.
HEAD.rad1strc_wave1strc=str2num(char(d1(27)));

%LINE28   Amplitude Factor Averaging Period min.
HEAD.Ampfac_aveper=str2num(char(d1(28)));


if length(d1)>29
     %LINE30  RadialFiller RCLim AngLim CurLim AngGap RngGap
HEAD.RadialFiller=str2num(char(d1(30)));
end


fclose(fid);


