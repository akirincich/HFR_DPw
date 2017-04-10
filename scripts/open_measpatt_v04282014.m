function [PATT]=open_measpatt(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function PATT=open_measpatt(filename);
%
%  function opens COS Seasonde MeasPattern.txt file for instrument in question
%  and loads all pertinant site info for computing radials from CSS file
%
%  returns structure PATT with all info inside
%
%   03/23/2012
%   Anthony Kirincich
%   WHOI
%   akirincich@whoi.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fid=fopen(filename);
d=textscan(fid,'%s','Delimiter','\n');
d1=d{1};
n=length(d1);
%%
PATT=[];
%Line 1
PATT.numbear=str2num(char(d1(1)));

%%%% get info on pattern
% first find out if the file format is from SSSv7 or SSSv6
%Line N
a=char(d1(n)); i=find(a=='!');
if strcmp(a(i+1:end),' Phase Corrections')==1; %this is a v7 file
    PATT.file_vers=7;
    
    %Line n-9
    a=char(d1(n-9)); i=find(a=='!');     PATT.ampfactors=str2num(a(1:i(1)-1));
    %Line n-8
    a=char(d1(n-8)); i=find(a=='!');     PATT.Antbear=str2num(a(1:i(1)-1));
    %Line n-7
    a=char(d1(n-7));                     PATT.Site_name=a(1:4);
    %Line n-6
    a=char(d1(n-6)); i=find(a=='!');     PATT.Site_loc=str2num(a(1:i(1)-1));
    %Line n-5
    a=char(d1(n-5)); i=find(a=='!');     PATT.DegRes=str2num(a(1:i(1)-1));
    %Line n-4
    a=char(d1(n-4)); i=find(a=='!');     PATT.DegSmooth=str2num(a(1:i(1)-1));
    %Line n-3
    a=char(d1(n-3)); i=find(a=='!');     PATT.Patt_date=datenum(str2num(a(1:i(1)-1)));
    %Line n-2
    PATT.Notes=char(d1(n-2));
    %Line n-1
    a=char(d1(n-1)); i=find(a=='!');     PATT.UUID=(a(1:i(1)-1));
    %Line n
    a=char(d1(n)); i=find(a=='!');       PATT.PhaseCorrects=str2num(a(1:i(1)-1));
    
elseif strcmp(a(i+1:end),' UUID')==1;  %This is a v6 file
    PATT.file_vers=6;
    
    %Line n-8
    a=char(d1(n-8)); i=find(a=='!');     PATT.ampfactors=str2num(a(1:i(1)-1));
    %Line n-7
    a=char(d1(n-7)); i=find(a=='!');     PATT.Antbear=str2num(a(1:i(1)-1));
    %Line n-6
    a=char(d1(n-6));                     PATT.Site_name=a(1:4);
    %Line n-5
    a=char(d1(n-5)); i=find(a=='!');     PATT.Site_loc=str2num(a(1:i(1)-1));
    %Line n-4
    a=char(d1(n-4)); i=find(a=='!');     PATT.DegRes=str2num(a(1:i(1)-1));
    %Line n-3
    a=char(d1(n-3)); i=find(a=='!');     PATT.DegSmooth=str2num(a(1:i(1)-1));
    %Line n-2
    a=char(d1(n-2)); i=find(a=='!');     PATT.Patt_date=datenum(str2num(a(1:i(1)-1)));
    %Line n-1
    PATT.Notes=char(d1(n-1));
    %Line n
    a=char(d1(n)); i=find(a=='!');     PATT.UUID=(a(1:i(1)-1));
    
end


%%

%%%% get actual pattern by reading all the lines in the 
PATT.mpatt=nan*ones(9,PATT.numbear);
lcount=2;
for ii=1:9
    a=[];
    while length(a)<PATT.numbear
        a=[a str2num(char(d1(lcount)))];     
        lcount=lcount+1;
    end
    PATT.mpatt(ii,:)=a;
end

fclose(fid);
%%

