function [FOreg,FOregi,HEAD]=HFR_spectrsrc_load(fname,scripts_dir,CSS_Head,HEAD)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [FOreg,FOregi,HEAD]=HFR_spectrsrc_load_v?(fname,scripts_dir,CSS_Head,HEAD)
%
%  This function reads the resource fork of the CSS file in question in
%  order to find the COS-estimated (via CSPro) FOLs (or ALims).
%
%  This function calls a perl script written by Tom Cook (SIO) to strip the
%  Resource fork from the file as this cannot be done in matlab.  In theory
%  this might work on a linux distribution, but has only been tested on a
%  mac.  Kudos to Tom Cook for creating and sharing the perl script.
%  However, Tom asks that you not contact him with any questions regarding
%  out to run it with other distributions.  If you are running this on a
%  Linux, you should be able to figure it out.  
%
% Version:
% v1   2014    created
% v2   3/2017  edited for the distribution package
%
%  Anthony Kirincich
%  WHOI PO
%  akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% mark that we are using this file to process...
HEAD.ProcessingSteps{end+1}=mfilename;


%%%%% continue and get resource fork information.
%call perl script by Tom Cook (SIO) 
s=['!perl ' scripts_dir 'getAlims-MAC.pl ' fname ' > temp.txt'];
eval(s);
FOLS=load('temp.txt','-ascii');
eval(['!rm temp.txt'])

%change FOLS to FOregi...
FOregi=[]; FOreg=zeros(CSS_Head.nRangeCells,CSS_Head.nDopplerCells);
for ii=1:length(FOLS);
    a=[FOLS(ii,2):FOLS(ii,3)    FOLS(ii,4):FOLS(ii,5)]';
    aa=[ones(length(a),1).*FOLS(ii,1)+1 a];    
    FOregi=[FOregi; aa];
    
    FOreg(ii,FOLS(ii,2):FOLS(ii,3))=1;
    FOreg(ii,FOLS(ii,4):FOLS(ii,5))=1;
end




