% SG_PROCESSSCHLA
% -------------------------------------------------------------------------
% Process Sea-Bird/WET Lab ECO Puck Volume chla,
% apply quality control flags, and calibrate with HPLC or FCHLA
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
% function sg_processBbp(missionname)
%
% PROCESSING STEPS
% -------------------------------------------------------------------------
% 1. Load data file
% 2. Apply quality control flags
%   - visual qc (after all qc, remaining bad data outliers)
%   - stuck qc (argo test, 50% cutoff)
%   - range qc (glogal argo test, regional HOT test)
%   - biofoul qc (from oxygen, bbp code)
% 3. Briggs 2011 spike code applied
% 4. Chl calibration
% 5. choose one chl1,chl2 to save as chl_cal from rmsd
% 6. Save tables to qcdata folder 'data/qcdata/missionname_L2.mat'

% MATLAB CODE NEEDED
% -------------------------------------------------------------------------
% sg_chla_calibrationTable
% sg_chla_calibrationMethod
% sg_greylistQC
% sg_spikes
%
% -------------------------------------------------------------------------
% Author Catherine Garcia -- cathy.garcia@hawaii.edu -- Version 04 Aug 2023
%
% -------------------------------------------------------------------------


function sg_processChla(missionname)

filename = ['data/rawdata_L1/' missionname '_L1.mat'] ; %file path to mission process data
load(filename,'DWN0','UP0');
filename = ['data/qcdata_L2/' missionname '_L2.mat'] ; %file path to mission process data
load(filename,'DWN','UP','calib','QC_UP','QC_DWN');

%% -------------------------------------------------------------------------
% QC Analysis on raw data

% 1. visual qc
if any(strcmp(DWN.Properties.VariableNames,'chl1'))
    [qcup,qcdown] = sg_greylistQC(missionname,UP,DWN,'chl1');
    QC_UP.qcflag_visual_chl1 = qcup;
    QC_DWN.qcflag_visual_chl1 = qcdown;
end
if any(strcmp(DWN.Properties.VariableNames,'chl2'))
    [qcup,qcdown] = sg_greylistQC(missionname,UP,DWN,'chl2');
    QC_UP.qcflag_visual_chl2 = qcup;
    QC_DWN.qcflag_visual_chl2 = qcdown;
end

% 2. Optics Stuck Value test
if any(strcmp(DWN.Properties.VariableNames,'chl1'))
    % chl1
    [qcup,qcdown] = sg_opticsQC_stuckValue(DWN0,UP0,'wlbb2fl_chl_sig');
    QC_UP.('qcflag_stuck_chl1') = qcup;
    QC_DWN.('qcflag_stuck_chl1') = qcdown;
end
if any(strcmp(DWN.Properties.VariableNames,'chl2'))
    %chl2
    [qcup,qcdown] = sg_opticsQC_stuckValue(DWN0,UP0,'wlbbfl2_chl_sig');
    QC_UP.('qcflag_stuck_chl2') = qcup;
    QC_DWN.('qcflag_stuck_chl2') = qcdown;
end

%*****************Range flags moved to post calibration******************%
% % 3. Flag all variables (and related ones) as fail for out of acceptable
% % range
% % global range test from argo float QC/QC DAC
% % chl1
% 
% QC_UP.qcflag_range_chl1 = ones(height(UP),1);
% QC_DWN.qcflag_range_chl1 = ones(height(DWN),1);
% QC_UP.qcflag_range_chl1(UP.chl1 < -0.1 | UP.chl1 > 50) = 4; % chl (mg m-3_
% QC_DWN.qcflag_range_chl1(DWN.chl1 < -0.1 | DWN.chl1 > 50) = 4; % chl (mg m-3_
% % chl2
% QC_UP.qcflag_range_chl2 = ones(height(UP),1);
% QC_DWN.qcflag_range_chl2 = ones(height(DWN),1);
% QC_UP.qcflag_range_chl2(UP.chl2 < -0.1 | UP.chl2 > 50) = 4; % chl (mg m-3)
% QC_DWN.qcflag_range_chl2(DWN.chl2 < -0.1 | DWN.chl2 > 50) = 4; % chl (mg m-3)
%*****************Range flags moved to after calibration******************%

% 4. Assign variable qc flags based on above tests
if any(contains(DWN.Properties.VariableNames,'chl1'))
    % chl1
    flags = {'qcflag_visual_general','qcflag_visual_chl1','qcflag_biofoulOptics','qcflag_stuck_chl1'};
    var_flags = QC_UP{:,flags};
    QC_UP.chl1_qcflag = max(var_flags,[],2,'omitnan');
    var_flags = QC_DWN{:,flags};
    QC_DWN.chl1_qcflag = max(var_flags,[],2,'omitnan');
end
if any(contains(DWN.Properties.VariableNames,'chl2'))
    % chl2
    flags = {'qcflag_visual_general','qcflag_visual_chl2','qcflag_biofoulOptics','qcflag_stuck_chl2'};
    var_flags = QC_UP{:,flags};
    QC_UP.chl2_qcflag = max(var_flags,[],2,'omitnan');
    var_flags = QC_DWN{:,flags};
    QC_DWN.chl2_qcflag = max(var_flags,[],2,'omitnan');
end

%% -------------------------------------------------------------------------
% Get spikeflag for chl uncal data
opticVars = {'chl1','chl2'};
for ii = 1:numel(opticVars)
    if sum(strcmp(DWN.Properties.VariableNames,opticVars{ii})) == 1
        % DWN
        [sgspike] = sg_spikes(DWN,opticVars); %DWN despiked
        DWN.([opticVars{ii} '_base']) = sgspike.([opticVars{ii} '_base']);
        DWN.([opticVars{ii} '_spikes']) = sgspike.([opticVars{ii} '_spikes']);
        DWN.([opticVars{ii} '_despiked']) = sgspike.([opticVars{ii} '_nospike']);
        % UP
        [sgspike] = sg_spikes(UP,opticVars); %UP despiked
        UP.([opticVars{ii} '_base']) = sgspike.([opticVars{ii} '_base']);
        UP.([opticVars{ii} '_spikes']) = sgspike.([opticVars{ii} '_spikes']);
        UP.([opticVars{ii} '_despiked']) = sgspike.([opticVars{ii} '_nospike']);
    end
end



%% -------------------------------------------------------------------------
% calibrate chlorophyll

% get calibration points
inputs = struct();
inputs.distance = 10E3;
inputs.depthz = 10;
inputs.timelag = 2/24;
inputs.calibration_filepath = 'glider-kit/calibrationFiles/hplc_cal_19-Jul-2023.mat';
inputs.method = 'linear';
[settings] = sg_chla_calibrationSettings(missionname);

% chl2
sg_chlvar = 'chl2_despiked';
chl2_comp = nan;
if sum(strcmp(DWN.Properties.VariableNames,'chl2')) > 0
    % get calibration matches & calibrate glider chlorophyll
    [chl2_comp] = sg_chla_calibrationTable(DWN(QC_DWN.chl2_qcflag==1,:),UP(QC_UP.chl2_qcflag==1,:),sg_chlvar,settings);
    if istable(chl2_comp)
        imatches = ~isnan(chl2_comp.sg_chl);
        if nnz(~isnan(chl2_comp.chla(imatches))) > 0
            inputs.chlvar = 'chla';
            [slope,intercept,~] = sg_chla_calibrationMethod(chl2_comp,inputs); % hplc chl exists
        else
            inputs.chlvar = 'fchla85m';
            [slope,intercept,~] = sg_chla_calibrationMethod(chl2_comp,inputs); % only fchla, use corrected version
        end
        % save figure
        figdir = 'figures/chl2_cal/';
        if ~isfolder(figdir); mkdir(figdir); end
        exportgraphics(gcf,[figdir missionname '.png'],'Resolution',300); close all;
        % calibrate chlorophyll
        UP.chl2_cal = UP.chl2.*slope + intercept;
        DWN.chl2_cal = DWN.chl2.*slope + intercept;
        chl2_comp.sg_chl_cal = chl2_comp.sg_chl.*slope + intercept; % mg m-3
        % calib ouptputs to save
        calib.chl2_cal.slope = slope;
        calib.chl2_cal.intercept = intercept;
        calib.chl2_comp = chl2_comp;
        calib.method = inputs.method;
    end
end

% chl1
sg_chlvar = 'chl1_despiked';
chl1_comp = nan;
if sum(strcmp(DWN.Properties.VariableNames,'chl1')) > 0
    if nnz(~isnan(DWN.chl1)) > 0
    % get calibration matches & calibrate glider chlorophyll
    [chl1_comp] = sg_chla_calibrationTable(DWN(QC_DWN.chl1_qcflag==1,:),UP(QC_UP.chl1_qcflag==1,:),sg_chlvar,settings);
    if istable(chl1_comp)
        imatches = ~isnan(chl1_comp.sg_chl);
        if nnz(~isnan(chl1_comp.chla(imatches))) > 0
            inputs.chlvar = 'chla';
            [slope,intercept,~] = sg_chla_calibrationMethod(chl1_comp,inputs); % hplc chl exists
        else
            inputs.chlvar = 'fchla85m';
            [slope,intercept,~] = sg_chla_calibrationMethod(chl1_comp,inputs); % only fchla, use corrected version
        end
        % save figure
        figdir = 'figures/chl1_cal/';
        if ~isfolder(figdir); mkdir(figdir); end
        exportgraphics(gcf,[figdir missionname '.png'],'Resolution',300); close all;
        % calibrate chlorophyll
        UP.chl1_cal = UP.chl1.*slope + intercept;
        DWN.chl1_cal = DWN.chl1.*slope + intercept;
        chl1_comp.sg_chl_cal = chl1_comp.sg_chl.*slope + intercept; % mg m-3
        % calib ouptputs to save
        calib.chl1_cal.slope = slope;
        calib.chl1_cal.intercept = intercept;
        calib.chl1_comp = chl1_comp;
        calib.method = inputs.method;
    end
    end
end


%*****************Range flags moved to post calibration******************%
% 3. Flag all variables (and related ones) as fail for out of acceptable
% range

% 4. Assign variable qc flags based on above tests
if any(contains(DWN.Properties.VariableNames,'chl1_cal'))
    % global range test from argo float QC/QC DAC
    QC_UP.qcflag_range_chl1 = ones(height(UP),1);
    QC_DWN.qcflag_range_chl1 = ones(height(DWN),1);
    QC_UP.qcflag_range_chl1(UP.chl1_cal < -0.1 | UP.chl1_cal > 50) = 4; % chl (mg m-3)
    QC_DWN.qcflag_range_chl1(DWN.chl1_cal < -0.1 | DWN.chl1_cal > 50) = 4; % chl (mg m-3)
    % chl1 flags
    flags = {'chl1_qcflag','qcflag_range_chl1'};
    var_flags = QC_UP{:,flags};
    QC_UP.chl1_qcflag = max(var_flags,[],2,'omitnan');
    var_flags = QC_DWN{:,flags};
    QC_DWN.chl1_qcflag = max(var_flags,[],2,'omitnan');
end
if any(contains(DWN.Properties.VariableNames,'chl2_cal'))
    % global range test from argo float QC/QC DAC
    QC_UP.qcflag_range_chl2 = ones(height(UP),1);
    QC_DWN.qcflag_range_chl2 = ones(height(DWN),1);
    QC_UP.qcflag_range_chl2(UP.chl2_cal < -0.1 | UP.chl2_cal > 50) = 4; % chl (mg m-3)
    QC_DWN.qcflag_range_chl2(DWN.chl2_cal < -0.1 | DWN.chl2_cal > 50) = 4; % chl (mg m-3)
    % chl2 flags
    flags = {'chl2_qcflag','qcflag_range_chl2'};
    var_flags = QC_UP{:,flags};
    QC_UP.chl2_qcflag = max(var_flags,[],2,'omitnan');
    var_flags = QC_DWN{:,flags};
    QC_DWN.chl2_qcflag = max(var_flags,[],2,'omitnan');
end
%*****************Range flags moved to after calibration******************%

% choose best chlorophyll (TBD: merge if one is partially missing?)
rmsd = nan(1,2);
if istable(chl1_comp)
    predicted = chl1_comp.sg_chl_cal;
    if strcmp(missionname,'sg511_21') || strcmp(missionname,'sg626_4')
        reference = chl1_comp.fchla;
    else
        reference = chl1_comp.chla.*1E-3; % convert chl units
    end
    inan = ~isnan(predicted) & ~isnan(reference);
    rmsd(1) = rms_dev(predicted(inan),reference(inan));
end
if istable(chl2_comp)
    predicted = chl2_comp.sg_chl_cal;
    if strcmp(missionname,'sg511_21') || strcmp(missionname,'sg626_4')
        reference = chl2_comp.fchla;
    else
        reference = chl2_comp.chla.*1E-3; % convert chl units
    end
    inan = ~isnan(predicted) & ~isnan(reference);
    rmsd(2) = rms_dev(predicted(inan),reference(inan));
end

% choose to report sensor with best coverage & no drift; if both appear to
% be good choose best match to calibration points
chl1_missionlist = {'sg146_9','sg148_6','sg148_11'};
chl2_missionlist = {'sg148_12','sg148_16','sg511_20','sg511_21','sg512_3','sg512_4','sg624_4'};
if any(contains(chl1_missionlist,missionname))
    ix = 1;
elseif any(contains(chl2_missionlist,missionname))
    ix = 2;
else
    [~,ix] = min(rmsd);
end

if ix == 1
    UP.chl_cal = UP.chl1_cal;
    DWN.chl_cal = DWN.chl1_cal;
    QC_UP.chl_cal_qcflag = QC_UP.chl1_qcflag;
    QC_DWN.chl_cal_qcflag = QC_DWN.chl1_qcflag;
    calib.chl_cal_var = 'chl1';
elseif ix == 2
    UP.chl_cal = UP.chl2_cal;
    DWN.chl_cal = DWN.chl2_cal;
    QC_UP.chl_cal_qcflag = QC_UP.chl2_qcflag;
    QC_DWN.chl_cal_qcflag = QC_DWN.chl2_qcflag;
    calib.chl_cal_var = 'chl2';
end


%% -------------------------------------------------------------------------
% Get spikeflag for chl cal data
opticVars = {'chl1_cal','chl2_cal','chl_cal'};
for ii = 1:numel(opticVars)
    if sum(strcmp(DWN.Properties.VariableNames,opticVars{ii})) == 1
        % DWN
        [sgspike] = sg_spikes(DWN,opticVars); %DWN despiked
        DWN.([opticVars{ii} '_base']) = sgspike.([opticVars{ii} '_base']);
        DWN.([opticVars{ii} '_spikes']) = sgspike.([opticVars{ii} '_spikes']);
        DWN.([opticVars{ii} '_despiked']) = sgspike.([opticVars{ii} '_nospike']);
        % UP
        [sgspike] = sg_spikes(UP,opticVars); %UP despiked
        UP.([opticVars{ii} '_base']) = sgspike.([opticVars{ii} '_base']);
        UP.([opticVars{ii} '_spikes']) = sgspike.([opticVars{ii} '_spikes']);
        UP.([opticVars{ii} '_despiked']) = sgspike.([opticVars{ii} '_nospike']);
    end
end


% save processed data
save(['data/qcdata_L2/' missionname '_L2'],'UP','DWN','calib','QC_DWN','QC_UP');

%}
end