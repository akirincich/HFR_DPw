function [R, HEAD]=HFR_DPmusic(Data,FOregi,HEAD,SpecHead,patt,CONST)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [R, HEAD]=HFR_DPmusic(Data,FOregi,HEAD,SpecHead,patt,CONST)
%
%  follows MUSIC algorithm (Schmidt,1986; Lipa et al,2006) to compute and 
%    evaluate the DOA function for each spectral point in FOregi for both 
%    single and Duel angle solutions.  For each MUSIC result, the DOA 
%      function is evalulated to find the significant peaks using
%      pksfinder and a super-complex set of conditional statments to
%      account for a variety of measured pattern types and azimuthal extents
%      (see below)
%
% INPUT:  Should be self-explanatory given the prepatory work by the
% spectra2radialmetric script that calls this function.
%
%outputs: 
%
%   R  --  A matrix of radial metric output.  This output
%          follows COS documentation for column/field order for lack 
%          of a better option. Outgoing data starts with the range 
%          (column 6) of COS radial metrics file. The rest (lon lat 
%          u v flag) are added later.
%
%         The number rows of R is variable and depends on the number of
%         FOregi and the percent of duel angle solutions returned.
%
%  columns of R are as follows
% 1     range
% 2     bearing
% 3     vel
% 4     direction
% 5     rangecell
% 6     dopcell
% 7     angselect
% 8-10  musicSnglang musicDuel1ang musicDuel2ang
% 11-13 musicEigenRatio musicpowerRatio musicoffRatio
% 14-16 musicSnglpow(v) musicDuel1pow(v) musicDuel2pow(v)
% 17-19 musicSngl_pkwidth musicDuel1_pkwidth musicDuel2_pkwidth
% 20-22 musicDOASnglpeak musicDOADuel1peak musicDOADuel2peak
% 23-25 spectraA1(snr) spectraA2(snr) spectraA3(snr)
% 26-28 Eigenval1 Eigenval2 Eigenval3
% 29    Dual_reject***
%
%  see COS documentation (or de Poala and Terril (2009?) for more information 
%   on the radial metrics output can be used for and what these 
%   identifiers mean.
%
%
% ***Note that the flag 'Dual_reject' is formulated differently here than
%    in the COS file format and combines four flags, each 2 digits in length,
%    into one, 8 digit number, as:
%
%    M.total_flag=MS_flag(1)*1e4 + MS_flag(2)*1e2 + (7-(M.crit*[1 2 4]'));
%
% where,
% MS_flag(1) -single angle solution result flag
%     if a 2 can be subtracted from the flag--the peak is on the right boundary
%     if a 1 can be subtracted from the flag--the peak is on the left boundary
%    
%MS_flag(2) -duel angle solution result flag
%     if a 8 can be subtracted from the flag--no midspan peaks
%     if a 4 can be subtracted from the flag--too many peaks found
%     if a 2 can be subtracted from the flag--a peak is on the right boundary
%     if a 1 can be subtracted from the flag--a peak is on the left boundary
%    
%
%  M.crit  -- the result of the application of the 3 MUSIC parameters for
%  defining wether to accept the single or duel angle solution (see Lipa et
%  al 2006 for details of the parameters and below for their application.
%
%  Key points:
% (1) If either of the MS_flag are ~=0, you should ignore the results, the
%     flags just tell you why.
%
% (2) only select the duel angle solution if the single angle solution is 
%     okay, all three criteria are met, AND none of the duel angles are 
%     at the boundary...
%
%
% Versions:
%  3/25ish?/2012
%   created
%
%   3/28/2012
%   upgraded to allow running with parfor
%
%  4/4/2012
%   -corrected snr calculation error
%   -found negative signal values...assume do to spectral creation
%   averaging,  non-physical...right?  data is eliminated in thresholding stage 
%     as ant3snr must be >5 
%   -still need to understand what a complex signal result actually means... 
%   
%
%  6/8/2016   eliminated parfor loop in favor of just for as some versions
%              of Matlab tend to have issues with the pool config.
%
%  v7  3/2017   cleaned up for package, 
%               fixed issue with patterns that wrap around 0 but are not full
%               fixed issue when the total range of the DOA function is less that 3 db
%                       to return the DOA peak witdh at max ([1 range/2]) down instead
%                       of three
%
%    Anthony Kirincich
%    WHOI-PO
%    akirincich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%start with array for Sing or Duel1 (R1) for all FOregi, and a second for Duel2
% condense and combine after finished
R1=nan.*ones(length(FOregi),29);   %output array for sing or Duel1
R2=nan.*ones(length(FOregi),29);   %output for Duel2


%prep for calculations
noise123=[median(Data.a1(length(SpecHead.rang),:)) median(Data.a2(length(SpecHead.rang),:)) median(Data.a3(length(SpecHead.rang),:))];
mP=HEAD.Musicparams123_gsmoothwidth_deg_smearwidth_velthresh(1:3);
pksfinder_factor=.1;
dangles=nanmedian(abs(diff(patt.angles)));


%%
%%% Run music for each spectral point within FOLs, output to RM format

for i4=1:length(FOregi); 
    
%%% this loop can also be run as a parallel job (if multiple workers exist)
%%% if running with parfor NO plotting is allowed with parfor
%%% CONST.goplot(1:2)=0;
%parfor i4=1:length(FOregi); %run music for each spectral point within FOLs

    %isolate point of interest for this loop
    bp=FOregi(i4,:);
    
    %set up output variables
    M=[];     MS=[];
    M.angles=[0 0 0];
    M.crit=[0 0 0];
    MS_flag=[0 0];
    
    %so data from css files give v1,v2,v3 (all squared) as well as v1v2, v1v3, and v2v3...
    % thus covariance matrix can be filled out right from css file info
    S=[Data.a1(bp(1),bp(2)) Data.a12(bp(1),bp(2)) Data.a13(bp(1),bp(2));
        conj(Data.a12(bp(1),bp(2))) Data.a2(bp(1),bp(2)) Data.a23(bp(1),bp(2))
        conj(Data.a13(bp(1),bp(2))) conj(Data.a23(bp(1),bp(2))) Data.a3(bp(1),bp(2))];

    if isempty(find(isnan(S(:))==1))==1  %fix to ensure that data in S is good
    
    %%%% equation to solve is c = a*s*a*T
    [Q ,D]=eig(S); %Compute eigendecomposition of covariance matrix
    [D,I]=sort(diag(D),1,'descend'); %Find r largest eigenvalues
    Q=Q (:,I); %Sort the eigenvectors to put signal eigenvectors first
    
    
    %compute music assuming both 1 and 2 signal-directions (DOAs) are present
    for r=1:2
        Qs=Q (:,1:r); %Get the signal eigenvectors
        Qn=Q(:,r+1:CONST.N); %Get the noise eigenvectors       
        music_spectrum=[];
        for k=1:length(patt.angles)
            %set antenna array...different for each num of signals
            a=[patt.a1(k); patt.a2(k); patt.a3(k)];
            %find the value for each direction.
            music_spectrum(k)=diag( [ (a'*a)/(a'*Qn*Qn'*a) ]);
        end
        MS(r,:)=real(music_spectrum);
    end
    
    %get arrival angles for both types.
    %identify global maxima
    itheta1=find(MS(1,:)==max(MS(1,:))); itheta1=itheta1(1);
    itheta2=find(MS(2,:)==max(MS(2,:))); itheta2=itheta2(1);
    
    %%%%%% find the DOA maximum   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% While this is easy for the single angle solution (the max value)
    %%% use pksfinder to find all local maxima and minima of MS(2,:) for
    %%% the duel angle solution
    %%%
    %%%  pksfinder also marks beginnng and end of array as peaks
    %%%  so if there is only 1 midspan peak, length(pks2)=3 
    %%%
    %%%  Use MS(1,:) range to set pksfinder blanking distance
    %%%  for MS(2,:)...can't have pksfactor be based on MS(2,:) b/c the
    %%%  range is often too high and masks 2nd peak
    %%%  this appears to make a hugh difference in radial coverage
    %%%   may want to consider doing more work on this
    [pks2,dzdt2]=pksfinder(MS(2,:),range(MS(1,:))*pksfinder_factor);
    if sum(diff(pks2)==1)>0   %midspan peaks are right next to each other
        %run pksfinder with larger range limit
        [pks2,dzdt2]=pksfinder(MS(2,:),range(MS(2,:))*pksfinder_factor);    
    end
    
    
    %%% Run through flags to see if the peak is legit  %%%
    %%% first, watch out for peak issues in MS1
    if (itheta1==1) %peak on left boundary
        MS_flag(1)=MS_flag(1)+1;
    end
    if (itheta1==length(patt.angles)) %peak on right boundary
        MS_flag(1)=MS_flag(1)+2;
    end

    %%% define how long the pattern is, if almost full coverage, treat
    %%% peaks near 0 or 360 differently.
    patt_coverage=length(patt.angles)./median(abs(diff(patt.angles)));  
    
    %%% second pick up extra peaks for MS2, and fill out error flag with info
    if patt_coverage < 345      %non full-coverage antenna pattern
        if length(pks2)>3 % multiple peaks exist
            %find where dzdt changes sign from pos to neg
            itheta2=(find(diff(dzdt2)==-2)+1)';  % dzdt sign changing from pos (increasing) to neg (decreasing)
            % length(itheta2)=2 here if 2 midspan peaks were found
            %if itheta still is <2...also find pks on the boundaries
            if length(itheta2)<2;
                if dzdt2(1)==-1   %mag decreases away from left boundary... boundary hides peak...make boundary a peak
                    itheta2=[itheta2 1];
                    MS_flag(2)=1;  %lhbd hides peak
                end
                if dzdt2(end)==1   %mag increases towards right boundary... boundary hides peak...make boundary a peak
                    itheta2=[itheta2 length(patt.angles) ];
                    MS_flag(2)=MS_flag(2)+2;  %rhbd hides peak
                end
            end
            if length(itheta2)>2  % too many peaks for this antenna array...take the top two in mag
                [s,i]=sort(MS(2,itheta2),'descend');  %sort by magnitude
                itheta2=itheta2(i(1:2));
                MS_flag(2)=MS_flag(2)+4;   %peaks >3 took top 2,, bad cut
            end
        elseif  length(pks2)<=3  %no midspan peak...
            MS_flag(2)=MS_flag(2)+8;
            itheta2=[itheta2(1) itheta2(1)];  %set both peaks to be the same angle, just to preserve output dimensions
        end
        
    elseif patt_coverage > 345  %almost full coverage, 
          %%% as wrap-arounds can exist, allow boundary peaks to exist without error
        
        if length(pks2)>3 % multiple peaks exist
            %find where dzdt changes sign from pos to neg
            itheta2=(find(diff(dzdt2)==-2)+1)';  % dzdt sign changing from pos (increasing) to neg (decreasing)
            %  length(itheta2)=2 here if 2 midspan peaks were found
            %if itheta still is <2...also find pks on the boundaries
            if length(itheta2)<2;
                if dzdt2(1)==-1   %mag decreases away from left boundary... boundary hides peak...make boundary a peak
                    itheta2=[itheta2 1];
                end
                if dzdt2(end)==1   %mag increases towards right boundary... boundary hides peak...make boundary a peak
                    itheta2=[itheta2 length(patt.angles) ];
                end
            end
            
            if length(itheta2)>2  % too many peaks for this antenna array...take the top two in mag
                %%%%% BUT, if antenna pattern is 360 coverage, check for wrap-around
                %%%%%  (real peak could be on boundary and counted twice)
                %%%%%     pick the highest of boundary peaks only
                ia2=find(itheta2==1 | itheta2==length(patt.angles));  %bdry peaks
                ia1=find(itheta2~=1 & itheta2~=length(patt.angles));  %all others
                if length(ia2)==2;  %both boundaries are here;
                    %find the lower of the two boundary peaks and delete it.
                    [s,i]=sort(MS(2,itheta2(ia2)),'descend');
                    itheta2=itheta2([ia2(i(1)) ia1]);
                    %resort to insure that top 2 mag peaks are 1:2
                    [s,i]=sort(MS(2,itheta2),'descend');
                    itheta2=itheta2(i);
                elseif length(ia1)==2;    %two midspan peaks exist with other trash,
                    % save the midspans only
                    itheta2=itheta2(ia1);
                end
                %confirm that itheta2 is only 2 long
                itheta2=itheta2(1:2);
            end
        elseif  length(pks2)<=3  %no midspan peaks...
            MS_flag(2)=MS_flag(2)+8;
            itheta2=[itheta2(1) itheta2(1)];  %set both peaks to be the same angle, just to preserve output dimensions
        end
        
    end %if range(patt.angles)...
    
    %%% Again, MS_flag notation is %%%
    %for  MS_flag(1) -single angle solution
    %     if a 2 can be subtracted from the flag--the peak is on the right boundary
    %     if a 1 can be subtracted from the flag--the peak is on the left boundary
    
    %for  MS_flag(2) -duel angle solution
    %     if a 8 can be subtracted from the flag--no midspan peaks
    %     if a 4 can be subtracted from the flag--too many peaks found
    %     if a 2 can be subtracted from the flag--a peak is on the right boundary
    %     if a 1 can be subtracted from the flag--a peak is on the left boundary
    
    %if either of the flags is ~=0, you should ignore the results, the
    %flags just tell you why.
    
    %%%%% plot DOA results for this doppler/range point and pause, because
    %%%%% clearly you'd like to look at every one of the O(1000) DOA results
    if CONST.goplot(2)==1;
        figure(6); clf
        subplot(211);
        plot(patt.angles,real(MS(1,:)),'b.'); hg;  plot(patt.angles,real(MS(2,:)),'g.'); hg
        
        plot(patt.angles(itheta1),abs(MS(1,itheta1)),'ko','markerfacecolor','k');
        plot(patt.angles(pks2),abs(MS(2,pks2)),'bo','markerfacecolor','b');
    %    axis([min(patt.angles)-5 max(patt.angles)+5 0 abs(MS(2,itheta2(2)))*2])
        
        subplot(212)
        plot(patt.angles,dzdt2,'b.'); hg
        axis([min(patt.angles) max(patt.angles) -2 5])
        
        MS_flag
        pause      
    end
  
    %%
    
    %%%%%% continue to calculate power, and output variables %%%%
    %%%% in schmidt notation %%%
    %  S is the 3x3 covariance matrix
    %   So is the noise covariance matrix (also 3x3) but unknown...assume it is the identity matrix
    %    lambda_min is the noise eigenvalue
    
    %use found peak directions to calculate A and P (in schmidt) or S (in lipa et al)
    A1=nan*ones(3,1); A2=A1;
    A1=[patt.a1(itheta1(1)); patt.a2(itheta1(1)); patt.a3(itheta1(1))];
    %do same for MS2 for both 2 peaks and 1 peak
    if length(pks2)>3 %2 peaks
        A2=[patt.a1(itheta2(1)) patt.a1(itheta2(2)); patt.a2(itheta2(1)) patt.a2(itheta2(2)); patt.a3(itheta2(1)) patt.a3(itheta2(2))];
    elseif length(pks2)==3   %1 peak...this result cannot be accepted, but include to ensure output varable sizes are correct
        A2=[patt.a1(itheta2(1)); patt.a2(itheta2(1)); patt.a3(itheta2(1))];
    end
    
    So=[1 0 0; 0 1 0; 0 0 1]; %assume the noise covariance matrix is the identity matrix
    P1=((A1'*A1)^-1)*A1'*(S- D(3)*So)*A1*((A1'*A1)^-1);  %back out the signal power
    P2=((A2'*A2)^-1)*A2'*(S- D(3)*So)*A2*((A2'*A2)^-1);  %back out the signal powers
    
    %%% apply seasonde parameters to steer music to solution, returns
    %%% decision and criteria values. a value of 1 is a failure
    M.Mparam=[abs(D(1))./abs(D(2)) 0 0];
    %criteria 1
    if M.Mparam(1) < mP(1);    M.crit(1)=1;  end
    %criteria 2 and 3
    if length(P2)==2
        d=sort(abs(diag(P2)),'descend');
        M.Mparam(2:3)=[d(1)/d(2) abs(P2(1,1)*P2(2,2)) / abs(P2(1,2)*P2(2,1))];
        if M.Mparam(2) < mP(2);        M.crit(2)=1;    end
        if M.Mparam(3) > mP(3);        M.crit(3)=1;    end
    end
    
    %return the following data
    M.angles=[patt.angles(itheta1(1)) patt.angles(itheta2)];
    M.MS_flag=MS_flag;
    M.power=[P1 diag(P2)'];
    if length(diag(P2)')==1; M.power(3)=0; end
    M.DOApeak=abs([MS(1,itheta1(1)) MS(2,itheta2)]);
    M.Eigen=abs(D)';
    %find snr of each ant
    M.snrant123=[S(1,1)./noise123(1) S(2,2)./noise123(2) S(3,3)./noise123(3)];
    
    %find the 3dB width of each peak
    %  for the single angle
    i=find(patt.angles==M.angles(1));
    j=find(real(MS(1,:))>real(MS(1,i))-3);
    %%% sometimes, the total range of the DOA function for the single ang
    %%% can be small, even though a peak appears to be well formed.
    if range(real(MS(1,:)))<3   %sometimes the whole thing is quite small, adjust the 3db
         j=find(real(MS(1,:))>real(MS(1,i))-max([1 range(real(MS(1,:)))./2]));
    end
    M.halfpowwidth(1)=abs(length(j)./dangles);    
    %  for the duel angle
    for j3=1:length(M.angles)-1;
        i=find(patt.angles==M.angles(j3+1));
        ii=find(abs(patt.angles-M.angles(j3+1))<51);
        j=find(real(MS(2,ii))>real(MS(2,i))-3);
        if range(real(MS(2,:)))<3   %sometimes the whole thing is quite small, adjust the 3db
          j=find(real(MS(2,:))>real(MS(2,i))-max([1 range(real(MS(2,:)))./2]));
        end     
        M.halfpowwidth(j3+1)=abs(length(j)./dangles);
    end
       
    %if M.halfpowwidth(1)>200;  asdfasfasfd; end
    
    %combine all flags together (each get 2 digits)
    M.total_flag=MS_flag(1)*1e4 + MS_flag(2)*1e2 + (7-(M.crit*[1 2 4]'));
    % if the flag is >0 reject the duel angle solution,
    % if the flag is >1e4 reject the single angle solution as well.
    
    %only select the duel angle solution if the single angle solution is okay,
    % all three criteria are met, AND none of the duel angles are at the boundary...
    M.angselect=0;   %no angle is selected (and no data is returned)
    if MS_flag(1)==0;
        M.angselect=1;
    end
    if M.total_flag==0;
        M.angselect=[2 3];
    end
        
    %%
    
    %%%% put outgoing data in to same array as COS radial metrics, starting
    %%%% with range (column 6) of COS radial metrics file.
    % add the rest (lon lat u v flag) after all FOLs are processed.
    %COS radial metrics  (lon lat u v flag will be added to the front end later
    % 1 range
    % 2 bearing
    % 3 vel
    % 4 direction
    % 5     rangecell
    % 6     dopcell
    % 7     angselect
    % 8-10   musicSnglang musicDuel1ang musicDuel2ang
    % 11-13   musicEigenRatio musicpowerRatio musicoffRatio
    % 14-16 musicSnglpow(v) musicDuel1pow(v) musicDuel2pow(v)
    % 17-19 musicSngl_pkwidth musicDuel1_pkwidth musicDuel2_pkwidth
    % 20-22 musicDOASnglpeak musicDOADuel1peak musicDOADuel2peak
    % 23-25 spectraA1(snr) spectraA2(snr) spectraA3(snr)
    % 26-28 Eigenval1 Eigenval2 Eigenval3
    % 29    Dual_reject*** (see above for explaination
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     if M.angselect(1)==1 || M.angselect(1)==2;
        j=1;
            d=[SpecHead.rang(bp(1)) M.angles(M.angselect(j)) SpecHead.c_vel(bp(2)) M.angles(M.angselect(j))+180];    if d(4)>360; d(4)=d(4)-360; end
            R1(i4,:)=[d bp M.angselect(j) M.angles M.Mparam M.power M.halfpowwidth M.DOApeak M.snrant123 M.Eigen M.total_flag];
     end
     if length(M.angselect)==2 & M.angselect(2)==3;
         j=2;
            d=[SpecHead.rang(bp(1)) M.angles(M.angselect(j)) SpecHead.c_vel(bp(2)) M.angles(M.angselect(j))+180];    if d(4)>360; d(4)=d(4)-360; end
            R2(i4,:)=[d bp M.angselect(j) M.angles M.Mparam M.power M.halfpowwidth M.DOApeak M.snrant123 M.Eigen M.total_flag];
    end 

    
    end  %if isempty(find(isnan(S(:))==1))==0
    
end  %for i4=1:length(FORegi)
%%

%combine R1 and R2, while getting rid of empty lines (marked by nans)
i1=find(isnan(R1(:,1))==0); i2=find(isnan(R2(:,1))==0);
R=[R1(i1,:); R2(i2,:)];
clear R1 R2

%%

%%%%% add this script to the processing steps
HEAD.ProcessingSteps{end+1}=mfilename;

return
