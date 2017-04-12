function [RADIAL,HEAD]=HFR_DP_RadialAverage_v5(RM,patt,HEAD,CONST);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [RADIAL,HEAD]=HFR_DP_RadialAverage_v5(RM,HEAD,CONST);
%
%  This script takes the Radial Metric output from HFR_DP and calculates
%   spatial averages over the user-defined azimuthal limits.  Note that
%   there is no temporal averaging available in the present package.
%
%  The averaged result is returned in a Matlab structure, named RADIAL,
%   that is a modified version of the HFR_Progs 'radial' structure to account
%   for the addition information returned with the radial metrics.
%
%  Two versions of the spatially averaged radials are reported:
%
%  (1) Arithmetic spatial means that replicates the steps made
%      within the CODAR-type processing (Lipa et al 2006)
%
%  (2) Spatial means compute from QC'ed (threshold filtering) and MUSIC-derived
%      Power-weighted averaging following the work of Kirincich et al (2012)
%       and de Poalo et al (2015).   The threshold values are user-settable
%       and exist within the CONST.radave_thresholds set in the master
%       program.
%
%  NOTES:
%    (A) This script leans heavily on the HFR_Progs toolbox and recreates/modifies
%    some of the HFR_Progs scripts to handle the radial metrics output.  Thanks
%     Mike Cook and David Kaplan for their amazing initial efforts
%     developing the HFR_Progs toolbox and distributing it widely.
%
%    (B)  Background information on the file formats and the fields are
%    given in the script before the processing is done.
%
%  VERSIONS:
%     v2 ...does not save RAM for each, no plotting
%
%     v3 ...various fixes including: decrease to just two aves, correctly
%           get bearings around zero (for asit)
%
%     v4 ...(01/20/2014)
%       -corrected the returning raderr (sample std dev) to return std dev
%       and not variance as was done in v3 for the nvel_est=2 option
%       -corrected the weighted raderr std dev to be consistent with the
%       normalized weight, w1.
%
% v5    March 2017   made into a function for the package...
%                    corrected issues related to how patt.Bear_ends is used
%                    to find the averaging window, 
%
% Anthony Kirincich
% Woods Hole Oceanographic Institution
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%% background info on RM formats and what all the fields mean:  %%%%%%
%%% this is not critical for the script to run, but does give an interested
%%% user some background knowledge on what all the fields of the RM file
%%% are and how they might be able to be appled here.  Most of the info is
%%% found in the COS documention,
%
% Again, column output format of RM(jjj).data follows COS RSv7 radial metric
%  file format with a few small changes to units and field content
%   this is done for convienance and consistency
%
% cols.   Field
% 1-2     lat lon
% 3-4     u v   (here nan as will not be used)
% 5       flag    (here nan, used by COS but not here)
% 6       range
% 7       bearing
% 8       vel  ( in m/s, not cm/s )
% 9       direction
% 10      rangecell
% 11      dopcell
% 12      angselect
% 13-15   musicSnglang musicDuel1ang musicDuel2ang
% 16-18   musicEigenRatio musicpowerRatio musicoffRatio
% 19-21   musicSnglpow(v) musicDuel1pow(v) musicDuel2pow(v)      (in voltage (v) not power (db) )
% 22-24   musicSngl_pkwidth musicDuel1_pkwidth musicDuel2_pkwidth
% 25-27   musicDOASnglpeak musicDOADuel1peak musicDOADuel2peak
% 28-30   spectraA1(snr) spectraA2(snr) spectraA3(snr)
% 31-33   Eigenval1 Eigenval2 Eigenval3
% 34      Dual_reject  (still different from COS methods)
%
%
%
%%% Another way to describe these RM-specific fields follows the Radial metrics
%%%  documentation of COS by Bill Rector:
%
% %SPRC is the spectra range cell (from 0 which is at site)
% %SPDC is the spectra doppler index (from 0 left side/ maximum negative Bragg.)
%
% %MSEL is the MUSIC bearing selected (1=single, 2= dual angle1, 3=dual angle2)
% %MSA1 is MUSIC single angle bearing. (1440 if invalid)
% %MDA1 is MUSIC dual angle first bearing (1440 if invalid)
% %MDA2 is MUSIC dual angle second bearing (1440. if invalid)
%
% %MEGR is MUSIC eigein ratio of first to second.
% %MPKR is MUSIC dual angle power signal ratio of first to second.
% %MDFR is MUSIC dual angle off diagonal ratio.
%
% %MSP1  is MUSIC Single angle power in dB (10*log10)
% %MDP1 is MUSIC Dual angle first power in dB
% %MDP2 is MUSIC Dual angle second power in dB
%
% %MSW1 is MUSIC Single angle antenna response angular width of the 3dB down from response peak bearing.
% %MDW1 is is MUSIC Dual first angle antenna response angular width of the 3dB down from response peak bearing.
% %MDW2 is is MUSIC Dual second angle antenna response angular width of the 3dB down from response peak bearing.
%
% %MSR1 is the MUSIC Single angle antenna response peak value
% %MDR1 is the MUSIC Dual first angle antenna response peak value
% %MDR2 is the MUSIC Dual second angle antenna response peak value
%
% %MA1S is the spectra antenna Loop1 signal to noise in dB
% %MA2S is the spectra antenna Loop2 signal to noise in dB
% %MA3S is the spectra antenna monopole signal to noise in dB
%
% %MEI1 is the MUSIC single angle eigen value
% %MEI2 is the MUSIC dual angle first eigen value
% %MEI3 is the MUSIC dual angle second eigen value
%
% %MDRJ is the MUSIC dual angle rejection reason(s).This value is a composite sum of powers of 2 (think binary) where each power of 2 sum is a different indicator. Possible flags are
% %	If 0, then solution is OK to be dual angle.
% %	+1 failed eigen ratio check.
% %	+2 failed signal ratio check.
% %	+4 failed off diagonal ratio check.
% %	+8 failed MDA1 is too close to MDA2
% %	+16 only one dual angle bearing solution found.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Unpack processing thresholds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snr_thresh=CONST.radave_thresholds(1);
angpeak_thresh=CONST.radave_thresholds(2);
angwidth_thresh=CONST.radave_thresholds(3:4);

%%% define the list of ALL fields found within the data matrix coming from RM
vv={'LOND' 'LATD' 'VELU' 'VELV' 'VFLG' 'RNGE' 'BEAR' 'VELO' 'HEAD' 'SPRC' 'SPDC' 'MSEL' 'MSA1' 'MDA1' 'MDA2' 'MEGR' 'MPKR' 'MOFR' 'MSP1' 'MDP1' 'MDP2' 'MSW1' 'MDW1' 'MDW2' 'MSR1' 'MDR1' 'MDR2' 'MA1S' 'MA2S' 'MA3S' 'MEI1' 'MEI2' 'MEI3' 'MDRJ'};
%%% Basic List of 'radial' columns
cc = { 'LOND', 'LATD', 'RNGE', 'BEAR', 'HEAD', 'VELO' };
%%% full list of the non-radial fields that come out of RM
mm={'SPRC','SPDC','MSEL','MSA1','MDA1','MDA2','MEGR','MPKR','MOFR','MSP1','MDP1', 'MDP2','MSW1','MDW1','MDW2','MSR1','MDR1','MDR2','MA1S','MA2S','MA3S','MEI1','MEI2','MEI3','MDRJ'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Other parameters to set, unlikely to change
%%% Set the Number of different estimates of the means to make
%nvel_est=4;
nvel_est=2;  %only return the power weighted average and straight mean

%%% set the file index to process
i5=1;  % I'm only going to do this for one file at a time, although you could
%  recast this function to do it for an array of RM structures with
%  many files. See the way Kaplan treats some of the radial
%  functions in HFR_Progs for an example.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% place the RM data into a modified RADIAL structure
R = RADIALstructv7;   %get a blank structure.
R.FileName = RM(i5).fname;
%print to screen so we know what is happening
fprintf('Processing radial averages for RM of: %s\n',R.FileName)
% Give TimeStamp and matvars same size as FileName for consistency just
% in case we return early.
R.TimeStamp = RM(i5).time;   %time
R.TimeZone = 'GMT';          %timezone

% for m = RADIALmatvars
%     R.(m{:}) = zeros(0,1);
% end

% Add processing step to the structure
R.ProcessingSteps{end+1} = [RM(i5).fname ' to RM with HFR_spectra2radialmetrics_process_v?'];
R.SiteOrigin = RM(i5).Site_loc([2 1]);   % Origin Lat,Lon,
R.SiteName = RM(i5).Site_name;           % Site name
R.Type = 'Meas';                         % PatternType

%%% find where the data of interest is:
%%%   This might seem overly complex now, but it will allow for the RM file
%%%   to be flexible in format without breaking the script.

%%% Get the indices of critical fields for radial info only
II = [];
for k = 1:length(cc);
    ii = strmatch( cc{k}, vv, 'exact' );
    if isempty(ii)
        feval( warnfunc, [ mfilename ':MISSING_ESSENTIAL_DATA_COLUMN' ], ...
            [ 'One or more of the required data columns cannot be found.' ...
            ' Returning empty structure.' ] );
        %  return
    else;  II(k) = ii;
    end
end

% Get pieces of data we want if they are not empty.
R.LonLat = RM(i5).data( :, II(1:2) );
R.RangeBearHead = RM(i5).data( :, II(3:5) );
R.RadComp = RM(i5).data(:,II(6))*100;  %convert to cm/s for rest of processing stream;

% % Change 999 in any field into NaN - NaN's create massive problems for
% % using error field to calculate totals error.
% for m = RADIALmatvars
%     m = m{:};
%     R.(m)( R.(m) == 999 ) = NaN;
% end

% Change direction to cartesian convention.  Note, however, that bearing
% will still point away from radar and heading will still point towards radar.
R.RangeBearHead(:,2:3) = true2math( R.RangeBearHead(:,2:3) );

%%% Get the indices of critical fields for radial metrics
JJ = [];
for k = 1:length(mm);
    ii = strmatch( mm{k}, vv, 'exact' );
    if isempty(ii)
        feval( warnfunc, [ mfilename ':MISSING_ESSENTIAL_DATA_COLUMN' ], ...
            [ 'One or more of the required data columns cannot be found.' ...
            ' Returning empty structure.' ] );
        %  return
    else; JJ(k) = ii;
    end
end

%output the whole thing as a large matrix
R.Metric=RM(i5).data(:,JJ);
% with the metric, adjust BEARINGS to be in math coords
R.Metric(:,4)= true2math( R.Metric(:,4));
R.Metric(:,5)= true2math( R.Metric(:,5));
R.Metric(:,6)= true2math( R.Metric(:,6));

% Add currently unused Flag (1 is for original radials).
R.Flag = ones( size(R.RadComp) );

%
% % % U and V are optional but generally around - compute if absent
% % II = strmatch( 'VELU', vv, 'exact' );
% % JJ = strmatch( 'VELV', vv, 'exact' );
% % if isempty(II) || isempty(JJ)
% %     feval( warnfunc, [ mfilename ':MISSING_UV_DATA_COLUMN' ], ...
% %         'U or V component missing - calculating from RadComp and Heading');
% %     [R.U,R.V] = deal( R.RadComp .* cosd(R.RangeBearHead(:,3)), ...
% %         R.RadComp .* sind(R.RangeBearHead(:,3)) );
% % else
% %     R.U = RM(i5).data(:,II);
% %     R.V = RM(i5).data(:,JJ);
% % end
% % %calculate if they were nans
% % if isempty(II)==0 && sum(~isnan(RM(i5).data(:,II)))==0
% %     feval( warnfunc, [ mfilename ':UV_DATA_COLUMN_NANS' ], ...
% %         'U or V component are nans - calculating from RadComp and Heading');
% %     [R.U,R.V] = deal( R.RadComp .* cosd(R.RangeBearHead(:,3)), ...
% %         R.RadComp .* sind(R.RangeBearHead(:,3))  );
% % end
%
% % %
% % % Deal with two possible names for error column.
% % I1 = strmatch( 'ETMP', vv, 'exact' );  %a normal radial file would have this.
% % I2 = strmatch( 'STDV', vv, 'exact' );  %old name
% % I3 = strmatch( 'MSA1', vv, 'exact');   %a radial metric file would have this
% %
% % if ~isempty(I1)==1   %if error column is named ETMP
% %     R.Error = RM(i5).data(:,I2);
% % end
% % if ~isempty(I2)==1    %if error column is named STDV
% %     R.Error = RM(i5).data(:,I2);
% % end
% % if isempty(I1) & isempty(I2)
% %       feval( warnfunc, [ mfilename ':MISSING_ERROR_DATA_COLUMN' ], ...
% %              'Error column missing. Errors will all be NaN');
% %     R.Error = repmat(NaN,size(R.RadComp));
% % end
% % %%
% % if ~isempty(I3)==1    %radial metric file load the rest of the parameters.
% %
%
% % end  %%% ~isempty(I3)==1
%
% % % Change 999 in any field into NaN - NaN's create massive problems for
% % % using error field to calculate totals error.
% % for m = RADIALmatvars
% %     m = m{:};
% %     R.(m)( R.(m) == 999 ) = NaN;
% % end

%%

%%%%%%% new parts that convert the radial metric data to regular radials %%%%% %%%%%
if isempty(R.Metric)==0 & isempty(find(~isnan(R.Metric(:))==1))==0  %%% skip files that are empty
    %%
    %%%% doing this range binning separately for each file will allow variable
    %%%% ranges to be combined into the same timeseries...although later
    %%%% averaging will have to account for this as well.
    
    %%% make aver bin range points
    [range]=unique(R.RangeBearHead(:,1));
    %for the Bear and range...make to latlon using m_fdist
    FLonLat=[]; FRangeBearHead=[];
    a=1;
    for j=1:length(patt.BearT)
        for i=1:length(range)
            [FLonLat(a,1),FLonLat(a,2),FLonLat(a,3)]=m_fdist(R.SiteOrigin(1),R.SiteOrigin(2),patt.BearT(j),range(i).*1000);
            %also put out the current range and Bear
            FRangeBearHead(a,:)=[range(i) patt.BearT(j) FLonLat(a,3)];
            a=a+1;
        end
    end
    %%% sometimes the longitude comes out greater than +/180deg, im not
    %%% sure why.  and its not wrong mbut interfers with the plotting 
    %%%   fix
    i=find( FLonLat(:,1) >180);    FLonLat(i,1)=FLonLat(i,1)-360;
    FRangeBearHead(:,2)=true2math(FRangeBearHead(:,2));
    FRangeBearHead(:,3)=true2math(FRangeBearHead(:,3));
    
    %  %Test look at this file layout
    if CONST.goplot(2)==1
        figure(10); clf;
        plot(FLonLat(:,1),FLonLat(:,2),'r.');
        title([R.FileName ' Radave locations ';],'interpreter','none')
    end
    
    %%% do same for the bearing midpoints at the largest range.
    i=find(range==max(range));
    FmLonLat=[];
    a=1;
    for j=1:length(patt.BearT_ends)
%        [FmLonLat(a,1),FmLonLat(a,2),FmLonLat(a,3)]=m_fdist(R.SiteOrigin(1),R.SiteOrigin(2),patt.BearT_ends(j),range(i).*1000);
        [FmLonLat(a,1),FmLonLat(a,2),FmLonLat(a,3)]=m_fdist(R.SiteOrigin(1),R.SiteOrigin(2),patt.BearT(j),range(i).*1000);
        a=a+1;
    end
    %%% sometimes the longitude comes out greater than +/180deg, im not
    %%% sure why.  and its not wrong mbut interfers with the plotting 
    %%%   fix
    i=find( FmLonLat(:,1) >180);    FmLonLat(i,1)=FmLonLat(i,1)-360;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Start outgoing radial structure
    RADIAL(i5)=R;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% new section to process this data for the given patt file
    %
    % DOA Peak Value                                %msr1 mdr1 mdr2
    % DOA Function Width (at half power)            %msw1 mdw1 mdw2
    % Signal Power (a measure of the eigenvalue)    %msp1 mdp1 mdp2
    % Antenna 1,2,3 SNR for this doppler cell       %ma1s ma2s ma3s
    
    %Build output structure with gridded products....
    RADIAL(i5).LonLat=FLonLat(:,1:2);
    RADIAL(i5).RangeBearHead=FRangeBearHead;
    RADIAL(i5).RadComp=nan*ones(length(FLonLat),nvel_est); RADIAL(i5).Error=RADIAL(i5).RadComp;
    RADIAL(i5).U=RADIAL(i5).RadComp; RADIAL(i5).V=RADIAL(i5).RadComp;
    RADIAL(i5).Flag=nan*ones(length(FLonLat),2);
    RADIAL(i5).Metric=nan*ones(length(RADIAL(i5).LonLat),6);   %for carrying scripps-type average metrics forwards
    RADIAL(i5).ProcessingSteps(length(RADIAL(i5).ProcessingSteps)+1)={'weighted average in HFRRMetricLoad_process'};
    
    %%
    for jj=1:length(range);     %loop over each range cell,
        ii=find(R.Metric(:,1)==jj);   %find all data that falls in this range
        
        %find the bearing selected for each row of data...
        %   because the bearings of each of the angles in a duel angle solution are given...
        %   ...use Bearing_Select (column 3) to extract the MUSIC defined bearing
        % How?  make a mask for  M.Bearing that picks the defined
        % bearing only and gives mbearing, an array of the chosen
        % bearing for each data row matching the range cell...
        % all in three lines of code!
        a=R.Metric(ii,3)*[1 1 1];  b=ones(size(ii))*[1 2 3];
        c=a.*b-b.^2;  i=find(c~=0); c(i)=nan; c=c+1;
        mbearing=nanmean(R.Metric(ii,4:6).*c,2);
        
        RAM=[]; % clear a structure that would hold all the individal groups
        % of radial for this range circle
        %%% for each ave bearing window, in this range cell, place all the found velocities within a new matrix
        for kk=1:length(patt.Bear);
          %  if patt.Bear(kk)<2.5 | patt.Bear(kk)>357.5;    %%%allows the radial near mathdeg=0 to capture the 5deg swath needed
          %      i3=[find(mbearing<patt.Bear_ends(kk))' find(mbearing>=patt.Bear_ends(kk+1))'];
          %  else
               % i3=find(mbearing<patt.Bear_ends(kk) & mbearing>=patt.Bear_ends(kk+1));
               i3=find(mbearing<max(patt.Bear_ends(:,kk)) & mbearing>=min(patt.Bear_ends(:,kk)));  %account for new matrix form of 
          %  end
            
            if ~isempty(i3)==1;  %found music output for this Bear/Range
                RAM(kk).Range_cell=jj;
                RAM(kk).Bear=patt.Bear(kk);
                RAM(kk).vlist=length(i3);
                RAM(kk).LonLat=R.LonLat(ii(i3),:);
                RAM(kk).RangeBearHead=R.RangeBearHead(ii(i3),:);
                RAM(kk).RadComp=R.RadComp(ii(i3));
                RAM(kk).Metric=R.Metric(ii(i3),:);
                
                %%%% find the azimuthal deg grid point that matches this bin
                i2=find(RADIAL(i5).RangeBearHead(:,1)==range(jj));
                i22=find(RADIAL(i5).RangeBearHead(i2,2)==patt.Bear(kk));
                %so use i2(i22)
                
                %get power, volt, ang_width, and ang_peak for each velocity in list
                a=RAM(kk).Metric(:,3)*[1 1 1];  b=ones(RAM(kk).vlist,1)*[1 2 3];
                c=a.*b-b.^2;  i=find(c~=0); c(i)=nan; c=c+1;
                %%% in the RM files 10:12 are the signal voltages, not powers
                volt=nanmean(abs(RAM(kk).Metric(:,10:12)).*c,2);
                power=10*log10(abs(volt)) + (-40 + 5.8);
                ang_width=nanmean(RAM(kk).Metric(:,13:15).*c,2);
                ang_peak=nanmean(RAM(kk).Metric(:,16:18).*c,2); % what is the real unit of this number?
                ant_snr=RAM(kk).Metric(:,19:21);
                
                %%%%%%% time to cut data with bad...
                %%%snr
                isnr=find(ant_snr(:,3) < snr_thresh);     %only use snr in ant3 as others might be low b/c bearing is in a null region of the antenna (and then it would be a good thing)
                ant_snr(isnr,:)=nan;
                isnr_a=find(isnan(mean(ant_snr,2))==0);
                isnr_b=find(isnan(mean(ant_snr,2))==1);  %for flag
                
                %%%ang_width
                iaw_a=find(ang_width > angwidth_thresh(1) & ang_width < angwidth_thresh(2));
                iaw_b=find(ang_width < angwidth_thresh(1) | ang_width > angwidth_thresh(2)); %for flag
                
                %%% ang_peak
                iap_a=find(ang_peak > angpeak_thresh);
                iap_b=find(ang_peak < angpeak_thresh);  %for flags
                
                %%%%%% find the radials that pass all the tests
                igood=intersect(isnr_a,intersect(iaw_a,iap_a));
                %%% make flags for this
                %RADIAL(i4).Flag(i2(i22),2)=[length(igood)*1e6 + length(isnr_b)*1e4 + length(iaw_b)*1e2+length(iap_b)];
                %%%  [ #of_ang_peak_cut (2digits)  #of_ang_width_cut (2digits)    #of_snr3_cut  # of igood (2digits)  ];
                RADIAL(i5).Flag(i2(i22),2)=[ length(iap_b)*1e6 + length(iaw_b)*1e4 + length(isnr_b)*1e2 + length(igood) ];
                
                if length(igood)>=1   %continue to make estimates of the average radial velocity
                    
                    if nvel_est==4;  %produce 4 different types of averages
                        %%establish data and weights arrays
                        d=RAM(kk).RadComp(igood);
                        w1=volt(igood)./sum(volt(igood));
                        %w2=ang_peak(igood)./sum(ang_peak(igood));         % other potential weights to use
                        %w2=ang_peakvolt(igood)./sum(ang_peakvolt(igood)); % other  potential weights to use
                        w2=ant_snr(igood,3)./sum(ant_snr(igood,3));
                        w3=(w1+w2)./2;
                        %average
                        a=[nanmean(d) sum(d.*w1) sum(d.*w2) sum(d.*w3)];
                        
                        if length(igood)>=3
                            %do std devation   %%%% see http://en.wikipedia.org/wiki/Weighted_mean
                            b= sqrt( [sum( (d-a(1)).^2)./(length(igood)-1) ...
                                sum( w1.*((d-a(2)).^2))./(1-sum(w1.^2)) ...
                                sum( w2.*((d-a(3)).^2))./(1-sum(w2.^2)) ...
                                sum( w3.*((d-a(4)).^2))./(1-sum(w3.^2)) ] );
                        else
                            b=[nan nan nan nan];
                        end
                        
                    elseif nvel_est==2;  %produce power-weighted and arthimetic averages only.
                        %%establish data and weights arrays
                        d=RAM(kk).RadComp(igood);
                        w1=volt(igood)./sum(volt(igood));
                        %average
                        a=[nanmean(d) sum(d.*w1)];
                        
                        if length(igood)>=3
                            %do std devation   %%%% see http://en.wikipedia.org/wiki/Weighted_mean
                            %%% w1 is normalized to sum(w1)=1, thus an 'unbiased' weighted variance
                            %%%      is not possible and sig^2 = sum(w1.*(d-a(2).^2)./V1  where V1 =sum(w1)=1
                            b= sqrt( [sum( (d-a(1)).^2)./(length(igood)-1) ...
                                sum( w1.*((d-a(2)).^2))./1 ] );
                        else
                            b=[nan nan];
                        end
                    end   %if nn=
                    
                    RADIAL(i5).RadComp(i2(i22),:)=a;
                    RADIAL(i5).Error(i2(i22),:)=b;
                    %whats the average solution used here 1,2, or 3?
                    RADIAL(i5).Flag(i2(i22),1)=nanmean(RAM(kk).Metric(igood,3));
                    %return other average metrics for later QA/QC
                    % (1) DOA Peak Value  (2) DOA Function Width (at half power) (3) Signal Power (4-6) ant1-3 snr
                    RADIAL(i5).Metric(i2(i22),:)=[mean(ang_peak(igood)) mean(ang_width(igood)) mean(10*log10(nanmean(volt(igood))))  mean(ant_snr(igood,:),1) ];
                    
                    
                    %%%% plot results, if we are in the weeds with looking
                    %%%% at all steps of the processing
                    if length(i3)>5 & CONST.goplot(2)==1;
                        
                        angle=mbearing(i3);
                        sp=RAM(kk).RadComp;
                        ant3snr=ant_snr(:,3);
                        %regroup by angle
                        [s si]=sort(angle);
                         
                        figure(1); clf
                        subplot(311);
                        plot(1:length(angle),angle(si),'ko','markerfacecolor','k'); hold on;
                        title('angle'); ylabel('deg'); %xlabel('# radials')
                        axis([0 length(angle)+1 min(angle)-1 max(angle)+1])
                        
                        
                        s1=subplot(312);
                        set(s1,'box','off','xlim',[0 length(angle)+1],'xtick',[0:length(angle)],...
                            'ylim',[0 max(ant3snr(si))])
                        l1=line(1:length(angle),ant3snr(si),'color','k','parent',s1); hold on
                        l2=line(1:length(angle),ant3snr(si),'color','k','marker','.','parent',s1);
                        title('ant3snr and speed')
                        ylabel('ant3snr');
                                                
                        a1=axes('position',get(s1,'position')); %speed
                        set(a1,'yaxislocation','right','xaxislocation','top','box','off',...
                            'color','none','xlim',[0 length(angle)+1],'xtick',[0:length(angle)],...
                            'ylim',[-35 35],'xticklabel',[])
                        l2=line(1:length(angle),sp(si),'color',[.5 .5 .5],'parent',a1); hold on;
                        l2=line(1:length(angle),sp(si),'color',[.5 .5 .5],'marker','.','parent',a1); hold on;
                        pp=line([0 length(angle)+1],[a(1) a(1)],'color','k','parent',a1(1));
                        
                        pp=line([0 length(angle)+1],[a(2) a(2)],'color','k','linestyle','--','parent',a1);
                        ylabel('speed cm/s');
                                                
                        s1=subplot(313);
                        set(s1,'box','off','xlim',[0 length(angle)+1],'xtick',[0:length(angle)],...
                            'ylim',[-150 -100])
                        l1=line(1:length(angle),power(si),'color','k','parent',s1); hold on
                        l2=line(1:length(angle),power(si),'color','k','marker','.','parent',s1);
                        title('power and speed')
                        ylabel('power dB');
                        xlabel('# radials');
                        
                        a1=axes('position',get(s1,'position')); %speed
                        set(a1,'yaxislocation','right','xaxislocation','top','box','off',...
                            'color','none','xlim',[0 length(angle)+1],'xtick',[0:length(angle)],...
                            'ylim',[-35 35],'xticklabel',[])
                        l2=line(1:length(angle),sp(si),'color',[.5 .5 .5],'parent',a1); hold on;
                        l2=line(1:length(angle),sp(si),'color',[.5 .5 .5],'marker','.','parent',a1); hold on;
                        pp=line([0 length(angle)+1],[a(1) a(1)],'color','k','parent',a1(1));
                        
                        pp=line([0 length(angle)+1],[a(2) a(2)],'color','k','linestyle','--','parent',a1);
                        ylabel('speed cm/s');
                        
                        pause  %if your plotting this, you wish to look at every single result
                    end  % if length(i3)  
                end  %if length(igood)==                
            end %isempty(i)
        end %for kk
        
    end %jj length range
    
    %%
    
    %%%% make U,V velocities for plotting %%%%
    H=math2true(RADIAL(i5).RangeBearHead(:,3));
    RADIAL(i5).U(:,:)=RADIAL(i5).RadComp(:,:).*sin(H*ones(1,nvel_est).*pi/180);
    RADIAL(i5).V(:,:)=RADIAL(i5).RadComp(:,:).*cos(H*ones(1,nvel_est).*pi/180);
    
    HEAD.ProcessingSteps{end+1}=mfilename;
    RADIAL.ProcessingSteps=HEAD.ProcessingSteps;
    
end  %isempty(R.Metric)==0

