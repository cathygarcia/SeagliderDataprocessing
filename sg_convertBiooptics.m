% SG_CONVERTBIOOPTICS
% -------------------------------------------------------------------------
% Converts Sea-Bird/WET Lab ECO bbp, chl, cdom counts to values 
%
% INPUTS
% -------------------------------------------------------------------------
% REQUIRED INPUTS:
% missionname:  glider missionname
%
% OPTIONAL INPUTS:
%
% OUTPUTS SAVED
% -------------------------------------------------------------------------
% UP:           output for ascending data
% DWN:          output for descending data
% calib:        calibration coefficient structure
% QC_DWN:       qcflag table for DWN table
% QC_UP:        qcflag table for UP table
% outputs saved to qcdata folder
%
% USAGE
% -------------------------------------------------------------------------
% function sg_convertBiooptics(missionname)
%
% PROCESSING STEPS
% -------------------------------------------------------------------------
% 1. Load data file
% 2. Apply scale and dark count offset
% 3. Save tables to qcdata folder 'data/qcdata/missionname_L2.mat'

% -------------------------------------------------------------------------
% Author Catherine Garcia -- cathy.garcia@hawaii.edu -- Version 04 Aug 2023
%
% -------------------------------------------------------------------------

function sg_convertBiooptics(missionname)

filename = ['data/rawdata_L1/' missionname '_L1.mat'] ; %file path to mission
load(filename,'DWN0','UP0');
filename = ['data/qcdata_L2/' missionname '_L2.mat'] ; %file path to mission process data
load(filename,'DWN','UP','calib','QC_UP','QC_DWN');

% chl1 bug, signal all 695, messes with stuck signal QC flag
DWN0.wlbb2fl_chl_sig(DWN0.wlbb2fl_chl_sig == 695) = nan;
UP0.wlbb2fl_chl_sig(UP0.wlbb2fl_chl_sig == 695) = nan;

% -------------------------------------------------------------------------
% Convert ECO Puck counts (Sea-Bird, formerly WetLabs) 
opticVars = {'wlbb2fl_bb470_sig','wlbb2fl_bb470_ref',...
            'wlbb2fl_bb700_sig','wlbb2fl_bb700_ref',...
            'wlbb2fl_chl_sig','wlbb2fl_chl_ref',...
            'wlbbfl2_bb1_sig','wlbbfl2_bb1_ref',...
            'wlbbfl2_chl_sig','wlbbfl2_chl_ref',...
            'wlbbfl2_cdom_sig','wlbbfl2_cdom_ref'};
[TF,Loc] = ismember(opticVars,DWN0.Properties.VariableNames);
dT = cat(1,DWN0(:,Loc(TF)),UP0(:,Loc(TF)));
directionIdx = [DWN.direction;UP.direction];

% BB2F-VMG and BBFL2-VMT
if isfield(calib,'wlbb2f_470_dark')
    dT.beta470 = calib.wlbb2f_470_scale.*(dT.wlbb2fl_bb470_sig - calib.wlbb2f_470_dark);
    dT.beta700 = calib.wlbb2f_700_scale.*(dT.wlbb2fl_bb700_sig - calib.wlbb2f_700_dark);
    dT.chl1 = calib.wlbb2f_chl_scale.*(dT.wlbb2fl_chl_sig - calib.wlbb2f_chl_dark);
end
% BBFL2-VMT and BBFL2-IRB
if isfield(calib,'wlbbfl2_chl_dark')
    if isfield(calib,'wlbbfl2_650_dark') % BBFL2-IRB
        dT.beta650 = calib.wlbbfl2_650_scale.*(dT.wlbbfl2_bb1_sig - calib.wlbbfl2_650_dark); % replaced dT.wl600sig with dT.wlbbfl2_bb1_sig (and bb600 with wlbbfl2_bb1) *** Changed by Cathy Garcia 09/2020
    elseif isfield(calib,'wlbbfl2_660_dark') % BBFL2-VMT
        dT.beta660 = calib.wlbbfl2_660_scale.*(dT.wlbbfl2_bb1_sig - calib.wlbbfl2_660_dark); % replaced dT.wl600sig with dT.wlbbfl2_bb1_sig (and bb600 with wlbbfl2_bb1)*** Changed by Cathy Garcia 09/2020
    end
    dT.cdom = calib.wlbbfl2_cdom_scale.*(dT.wlbbfl2_cdom_sig - calib.wlbbfl2_cdom_dark); % replaced dT.cdom_sig with dT.wlbbfl2_cdom_sig *** Changed by Cathy Garcia 09/2020
    dT.chl2 = calib.wlbbfl2_chl_scale.*(dT.wlbbfl2_chl_sig - calib.wlbbfl2_chl_dark); % replaced dT.chl_sig with dT.wlbbfl2_chl_sig *** Changed by Cathy Garcia 09/2020
end

% new 66.11 version has all cal data in 'WETLabsCalData' struct

bb2fl_list = {'bb470','bb700','chl'};
bbfl2_list = {'bb1','chl','cdom'}; %bbfl2_list = {'bb1','cdom','chl'}; 28 Jul 2021 CAG if chlorophyll a magnitude too high (0 - 5 range) check that this is in the correct order
bb2fl_out = {'beta470','beta700','chl1'};
bbfl2_out = {'beta650','chl2','cdom'};

if isfield(calib,'WETLabsCalData')
    WETLabsCalData = calib.WETLabsCalData;
    if isfield(WETLabsCalData,'wlbb2f') % Ignores the difference between bb2fl and bb2f sensors
        WETLabsCalData.wlbb2fl = WETLabsCalData.wlbb2f;
        WETLabsCalData = rmfield(WETLabsCalData,'wlbb2f');
    end
    sensorList = fieldnames(WETLabsCalData);
    nS = length(sensorList);
    for ii = 1:nS
        channelList = fieldnames(WETLabsCalData.(sensorList{ii}));
        nC = length(channelList);
        for jj = 1:nC
            if strcmpi(sensorList{ii},'wlbb2fl')
                channelVar = bb2fl_list{jj};
                outputVar = bb2fl_out{jj};
            elseif strcmpi(sensorList{ii},'wlbbfl2')
                channelVar = bbfl2_list{jj};
                outputVar = bbfl2_out{jj};
            else
                channelVar = channelList{jj};
            end
            dT.(outputVar) = ...
                WETLabsCalData.(sensorList{ii}).(channelList{jj}).scaleFactor...
                .*(eval(['dT.' sensorList{ii} '_' channelVar '_sig'])...
                - WETLabsCalData.(sensorList{ii}).(channelList{jj}).darkCounts);
        end
    end
end


% Split DWN, UP & join with raw data
dT = removevars(dT,opticVars(TF));
newvars = dT.Properties.VariableNames;
for ii = 1:numel(newvars)
    DWN.(newvars{ii}) = dT.(newvars{ii})(directionIdx == "down");
    UP.(newvars{ii}) = dT.(newvars{ii})(directionIdx == "up");
end


% save processed data
save(['data/qcdata_L2/' missionname '_L2'],'UP','DWN','calib','QC_DWN','QC_UP');

end