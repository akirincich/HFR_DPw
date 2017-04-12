function [FOreg,FOregi, Alims, HEAD]=HFR_spectrsrc_load(fname,scripts_dir,CSS_Head,HEAD,SpecHead)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [FOreg,FOregi,HEAD]=HFR_spectrsrc_load_v?(fname,scripts_dir,CSS_Head,HEAD)
%
%  This function reads the resource fork of the CSS file in question in
%  order to find the COS-estimated (via CSPro) FOLs (or ALims). This
%  function is based on a perl script written by Tom Cook (SIO) to strip the
%  Resource fork from the file.  Kudos to Tom Cook for creating and sharing
%  the original perl script.
%
%  The function first uses a shell command to strip the resource fork from
%  the spectral file and save it as a temporary file in the current
%  directory. It then opens and reads in the binary data in to get the
%  'Alims'  and converts the result to FOreg, etc.
%
%  In theory this might work on a linux distribution, but has only been
%  tested on a mac. Tom asks that you not contact him with any questions
%  regarding how to run it with other distributions.  If you are running
%  this on a Linux, you should be able to figure it out.
%
% Version:
% v1   2014    created
% v2   07012016
%       generalized for older .cs files where the range cells are defined
%       differently
% v3   3/2017  edited for the distribution package
%
% v4  4/2017 corrected for issue with corrupted Alims
%
%  Anthony Kirincich
%  WHOI PO
%  akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% mark that we are using this file to process...
HEAD.ProcessingSteps{end+1}=mfilename;


%%
% break resource fork out into a temp file
eval(['!cp ' fname '/..namedfork/rsrc temp']);

fid=fopen('temp','r','b');

%open(FILE,"<:raw",$rsrcfile);
fseek(fid,260,0);

%read(FILE,$buff,8);
%my ($tnam,$num_rc)=unpack "a4N",$buff;
tnam=char(fread(fid,4,'char'))';
num_rc = fread(fid,1,'uint');

%read(FILE,$buff,4);
%my $first_rc=unpack "f", $buff;
first_rc=fread(fid,1,'float');

%read(FILE,$buff,4);
%my $bearing=unpack "f", $buff;
bearing=fread(fid,1,'float');

%read(FILE,$buff,16);
%my ($start_rc,$r1,$r2,$r3)=unpack "NNNN",$buff;
start_rc = fread(fid,1,'uint');
r1 = fread(fid,1,'uint');
r2 = fread(fid,1,'uint');
r3 = fread(fid,1,'uint');


% 	for ($i=$start_rc;$i<=$num_rc;$i++){
% 	    read(FILE,$buff,16);
% 	    ($Ldlb[$j],$Rdlb[$j],$Ldrb[$j],$Rdrb[$j])=unpack "NNNN", $buff;
% 	    $Ldlb[$j]--;
% 	    $Rdlb[$j]--;
% 	    $Ldrb[$j]=$Ldrb[$j]+$nDopplerCells/2;
% 	    $Rdrb[$j]=$Rdrb[$j]+$nDopplerCells/2;
% 	    $j++;
% 	}

rangecells=start_rc:1:num_rc-1;

if start_rc(1)==0;  %older files
    rangecells=start_rc:1:num_rc-1;
    rangecells=rangecells+1;
elseif start_rc(1)==1;  %newer files...
    rangecells=start_rc:1:num_rc;
end

for i=1:length(rangecells)
    Ldlb(i)= fread(fid,1,'uint');
    Rdlb(i)= fread(fid,1,'uint');
    Ldrb(i)= fread(fid,1,'uint')+CSS_Head.nDopplerCells/2;
    Rdrb(i)= fread(fid,1,'uint')+CSS_Head.nDopplerCells/2;
end

%clean up
fclose(fid);
eval(['!rm temp']);

%%

%screen for ALims greater than a max value... default to Bragg

i=find(abs(Ldlb)>CSS_Head.nDopplerCells); Ldlb(i)=SpecHead.iFBragg(1);
i=find(abs(Rdlb)>CSS_Head.nDopplerCells); Rdlb(i)=SpecHead.iFBragg(1);

i=find(abs(Ldrb)>CSS_Head.nDopplerCells); Ldrb(i)=SpecHead.iFBragg(2);
i=find(abs(Rdrb)>CSS_Head.nDopplerCells); Rdrb(i)=SpecHead.iFBragg(2);


%%% Finish by exporting to spectral points
if isempty(num_rc)==0
    FOLS=[rangecells' Ldlb' Rdlb' Ldrb' Rdrb'];
    Alims=FOLS(:,2:5);
    
    FOregi=[]; FOreg=zeros(CSS_Head.nRangeCells,CSS_Head.nDopplerCells);
    for ii=1:length(FOLS);
        a=[FOLS(ii,2):FOLS(ii,3)    FOLS(ii,4):FOLS(ii,5)]';
        aa=[ones(length(a),1).*FOLS(ii,1) a];
        FOregi=[FOregi; aa];
        
        FOreg(ii,FOLS(ii,2):FOLS(ii,3))=1;
        FOreg(ii,FOLS(ii,4):FOLS(ii,5))=1;
    end
    
else  %    Alims=nan.*ones(1,5);
    %return error because you cant proceed with out the Alims.
    disp('error...no Alims were found in the spectral file');
    asdfasdfsdf
end


