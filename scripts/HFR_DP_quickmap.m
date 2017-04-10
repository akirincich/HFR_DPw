function HFR_DP_quickmap(ARRAY,working_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function HFR_DP_quickmap(ARRAY,working_dir);
%
% HFR_DP_quickmap creates a basic background map to be used for plotting radials
% or totals estimated using the HFR_DP package.
%
% The function requires the m_map package to exist with the matlab path and 
% for the structure ARRAY to be populated as described below, or in the
% HFR_DP_master_SITE.m setup file.
%
% ARRAY=[];
% ARRAY.name='xxxxx';
% %define the bounding box, needed for any plotting operations.
% ARRAY.min_lon=-(71+0/60); 
% ARRAY.max_lon=-(69+0/60); 
% ARRAY.min_lat=(40+0/60); 
% ARRAY.max_lat=(42+0/60);
%
% The function looks for a coastline file called 'xxxxx_coast.mat' within
% the 'working_dir', also defined in HFR_DP_master_SITE.m, to speed
% plotting operations. If missing, or for a new ARRAY, the function saves 
% a coastline mat file in the working_dir.
%
% NOTE: if the bounding box of ARRAY is changed, but not ARRAY.name, the
% user will need to delete the 'xxxxx_coast.mat' file and re-run the
% function.
%
%
% March 2017
% Anthony Kirincich
% WHOI
% akirincich@whoi.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%set the land color
landcolor=[245 222 179]./255;
seacolor=[224 255 255]./255;
%Get the bounding box of the array using
marea= [ARRAY.min_lon ARRAY.max_lon ARRAY.min_lat ARRAY.max_lat];


%start the plot on the current figure
m_proj('Mercator','lon',marea(1:2),'lat',marea(3:4));

% look for a coastline file named [ARRAY.name '_coast'] in the working directory
f=dir([working_dir '/' ARRAY.name '_coast.mat']);
% if no files exist make one.
if isempty(f)==1
    m_gshhs_h('save',[working_dir '/' ARRAY.name '_coast.mat'])
end

%add background seacolor patch
m_patch(marea([1 1 2 2 1]),marea([3 4 4 3 3]),seacolor)
hold on;
% add the coastline from this file as a patch
m_usercoast([ARRAY.name '_coast'],'patch',landcolor);
%finish
hold on
m_grid;

