% SG_PROCESSSAAOXYGEN
% -------------------------------------------------------------------------
% Process Aanderaa oxygen phase values, convert to units of 
% dissolved oxygen concentration, and apply quality control flags
%
% INPUTS
% -------------------------------------------------------------------------
% REQUIRED INPUTS:
% missionname:  glider missionname
%
% OPTIONAL INPUTS:
% tauoptode:    time in seconds for bittig inverse filter (DEFAULT = 30)
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
% function sg_processAAoxygen(missionname)
%
% PROCESSING STEPS
% -------------------------------------------------------------------------
% 1. Load raw data file "DWN0" & "UP0:
% 2. Calculate oxygen concentration
%   - recover mission sg146_9 phase if needed
%   - determine phase (calphase or dphase)
%   - pressure effect on lumiphore (Bittig 2015)
%   - calculate DO from phase
% 3. Apply quality control flags
%   - visual qc (after all qc, remaining bad data outliers)
%   - range qc (glogal argo test, regional HOT test)
%   - biofoul qc (from oxygen, bbp code)
%   - temp,salin qc input flag (o2 conversion uses these inputs)
% 4. Winkler calibration
% 5. Bittig 2014 inverse filter applied
% 6. Save tables to qcdata folder 'data/qcdata/missionname_L2.mat'
%
% MATLAB CODE NEEDED
% -------------------------------------------------------------------------
% aa_recoverRawPhase
% sg_aaoptode
% sg_greylistQC
% sg_oxygenWinklerCalibration
% optodeBittigCorrection
% -------------------------------------------------------------------------
% Author Catherine Garcia -- cathy.garcia@hawaii.edu -- Version 04 Aug 2023
%
% -------------------------------------------------------------------------

function sg_processAAoxygen(missionname,varargin)

if ~exist("tauoptode",'var')
    tauoptode = 30;
end

filename = ['data/rawdata_L1/' missionname '_L1.mat'] ; %file path to mission
load(filename,'DWN0','UP0');
filename = ['data/qcdata_L2/' missionname '_L2.mat'] ; %file path to mission process data
load(filename,'DWN','UP','calib','QC_UP','QC_DWN');

% -------------------------------------------------------------------------
% Aanderraa raw oxygen processing
dT = [DWN(:,{'tempc','sp','sa','ct','sigmath0','press','direction'});...
    UP(:,{'tempc','sp','sa','ct','sigmath0','press','direction'})];
dT.optode_temp = [DWN0.optode_temp;UP0.optode_temp]; % save as backup temperature

%
%  patch to recover mission sg146_9 dphase, phaseCoeff bad on calibration
%  sheet, from Roo Nicholson
%
if strcmp(missionname,'sg146_9') %if isfield(calib,'run_sg146_9_patch') % set in sg_calib_constant.m file for this mission only
    %if run_sg146_9_patch == 1 % set to true, run patch
        dT.optode_dphase = [DWN0.optode_dphase;UP0.optode_dphase];
        rphase = aa_recoverRawPhase(dT.optode_dphase, [-88.3602,3.65857]);
        NewPhaseCoef =  [-2.78588, 1.16177]; % old coefficients
        dT.optode_phase = rphase.*(NewPhaseCoef(2)) + NewPhaseCoef(1);
    %end
else
    % determine phase
    if isfield(calib,'optode_C00Coef')
        phase_var = 'optode_dphase';
    elseif isfield(calib,'optode_FoilCoefA0') || isfield(calib,'optode_Coef')
        phase_var = 'optode_calphase';
    end
    dT.optode_phase = [DWN0.(phase_var);UP0.(phase_var)];
end

% Phase pressure adjustment accounts for O-2 independent pressure effect on
% lumiphore (Bittig 2015)
%
% convert optode dphase or calphase to raw phase
if exist('optode_PhaseCoef0','var')
    pcoef1 = 0.1; %set pcoef1 = 0 for no pressure adjustment
    rawphase = (dT.optode_phase - optode_PhaseCoef0)./optode_PhaseCoef1;
    % Bittig et al. (2015) correction
    rawphase = rawphase + pcoef1*dT.press./1000;
    dT.optode_phase = optode_PhaseCoef0 + optode_PhaseCoef1.*rawphase;
else
    pcoef1 = 0; % for no pressure adjustment
end

% Calculate oxygen concentration (Aanderaa Optode)
sg_optodeCoef;
if isfield(calib,'optode_C00Coef')
    optodemod = 3830;
    %     AA.CM = [optode_C00Coef optode_C01Coef optode_C02Coef optode_C03Coef;...
    %         optode_C10Coef optode_C11Coef optode_C12Coef optode_C13Coef;...
    %         optode_C20Coef optode_C21Coef optode_C22Coef optode_C23Coef;...
    %         optode_C30Coef optode_C31Coef optode_C32Coef optode_C33Coef;...
    %         optode_C40Coef optode_C41Coef optode_C42Coef optode_C43Coef];
    dT.optode_oxygen = sg_aaoptode(optodemod,AA,dT.optode_phase,dT.tempc,dT.sp,dT.press,pcoef1); % dphase
elseif isfield(calib,'optode_FoilCoefA0')
    optodemod = 4330;
    %     AA.optode_FoilCoefA = optode_FoilCoefA;
    %     AA.optode_FoilCoefB = optode_FoilCoefB;
    %     AA.FoilPolyDegT = optode_FoilPolyDegT;
    %     AA.FoilPolyDegO = optode_FoilPolyDegO;
    if isfield(calib,'optode_S0')
        AA.S0 = optode_S0;
    end
    dT.optode_oxygen = sg_aaoptode(optodemod,AA,dT.optode_phase,dT.tempc,dT.sp,dT.press,pcoef1); % calphase
elseif isfield(calib,'optode_SVUCoef')
    optodemod = 4330;
    AA.SVU = optode_SVUCoef;
    dT.optode_oxygen = sg_aaoptode(optodemod,AA,dT.optode_phase,dT.tempc,dT.sp,dT.press,pcoef1); % calphase
else
    error('unidentified optode model "optodemod"');
end

% save dwn, up as new table
DWN.optode_oxygen = dT.optode_oxygen(dT.direction == "down");
UP.optode_oxygen = dT.optode_oxygen(dT.direction == "up");
DWN.optode_phase = dT.optode_phase(dT.direction == "down");
UP.optode_phase = dT.optode_phase(dT.direction == "up");
DWN.optode_temp = dT.optode_temp(dT.direction == "down");
UP.optode_temp = dT.optode_temp(dT.direction == "up");

% -------------------------------------------------------------------------
% QC Analysis on raw data

% 1. flag sbe aandera specifically
[qcup,qcdown] = sg_greylistQC(missionname,UP,DWN,'aanderaa');
QC_UP.qcflag_visual_aanderaa = qcup;
QC_DWN.qcflag_visual_aanderaa = qcdown;

% 2. Flag all variables (and related ones) as fail for out of acceptable
% range
% global range test from argo float QC/QC DAC
% preallocate
QC_UP.qcflag_range_oxygenAA = ones(height(UP),1);
QC_DWN.qcflag_range_oxygenAA = ones(height(DWN),1);
QC_UP.qcflag_range_phaseAA = ones(height(UP),1);
QC_DWN.qcflag_range_phaseAA = ones(height(DWN),1);
QC_UP.qcflag_range_tempAA = ones(height(UP),1);
QC_DWN.qcflag_range_tempAA = ones(height(DWN),1);
% regional ranges (run first to not overwrite global range)
QC_UP.qcflag_range_oxygenAA(UP.optode_oxygen < 0 | UP.optode_oxygen > 600) = 3; % oxygen regional range (umol L-1)
QC_DWN.qcflag_range_oxygenAA(DWN.optode_oxygen < 0 | DWN.optode_oxygen > 600) = 3; % oxygen regional range (umol L-1)
QC_UP.qcflag_range_phaseAA(UP.optode_phase < 10 | UP.optode_phase > 70) = 3; % oxygen phase regional range (deg)
QC_DWN.qcflag_range_phaseAA(DWN.optode_phase < 10 | DWN.optode_phase > 70) = 3; % oxygen phase regional range (deg)
QC_UP.qcflag_range_tempAA(UP.optode_temp < -2.5 | UP.optode_temp > 40) = 3; % temperature global range (deg C)
QC_DWN.qcflag_range_tempAA(DWN.optode_temp < -2.5 | DWN.optode_temp > 40) = 3; % temperature global range (deg C)
% global ranges (overwrites regional range to for large outliers)
QC_UP.qcflag_range_oxygenAA(UP.optode_oxygen < 0 | UP.optode_oxygen > 600) = 4; % oxygen global range (umol L-1)
QC_DWN.qcflag_range_oxygenAA(DWN.optode_oxygen < 0 | DWN.optode_oxygen > 600) = 4; % oxygen global range (umol L-1)
QC_UP.qcflag_range_tempAA(UP.optode_temp < -2.5 | UP.optode_temp > 40) = 4; % temperature global range (deg C)
QC_DWN.qcflag_range_tempAA(DWN.optode_temp < -2.5 | DWN.optode_temp > 40) = 4; % temperature global range (deg C)

% --- Patch for sg148_12 ---------------
% Low o2 data is real, likely a "cuddy" off oxygen minimim zone in east
% tropical pacific that moved north and west to hawaii
% other wise there is no data out of range
if strcmp(missionname,'sg148_12')
    QC_UP.qcflag_range_oxygenAA(QC_UP.qcflag_range_oxygenAA > 1) = 1;
    QC_DWN.qcflag_range_oxygenAA(QC_DWN.qcflag_range_oxygenAA > 1) = 1;
end

% 3. Biofouling evidence 
QC_UP.qcflag_biofoul = ones(height(UP),1);
QC_DWN.qcflag_biofoul = ones(height(DWN),1);
% only include data within regional phase range and that is not "generally"
% a bad mission
iflagDWN = QC_DWN.qcflag_visual_aanderaa < 3 & QC_DWN.qcflag_visual_general < 3 & QC_DWN.qcflag_range_phaseAA < 3;
iflagUP = QC_UP.qcflag_visual_aanderaa < 3 & QC_UP.qcflag_visual_general < 3 & QC_UP.qcflag_range_phaseAA < 3;
[biofouldate,pvalue] = sg_BiofilmEvidence(DWN(iflagDWN,:),UP(iflagUP,:),missionname);
calib.biofoul.dateHST = biofouldate;
calib.biofoul.pvalue = pvalue;
% only certain missions
if ~isnan(biofouldate)
    
    % set flag
    ibad = find(UP.dateHST >= biofouldate,1);
    QC_UP.qcflag_biofoul(ibad:end) = 4;
    ibad = find(DWN.dateHST >= biofouldate,1);
    QC_DWN.qcflag_biofoul(ibad:end) = 4; % set qc flag to 4 for fail

    % save figure
    figdir = 'figures/biofilm/';
    if ~isfolder(figdir); mkdir(figdir); end
    exportgraphics(gcf,[figdir missionname '.png'],'Resolution',300);
    close

end

% Spike test
dives = unique(DWN.divenum); % DWN
nd = numel(dives);
qcdown = ones(size(DWN.lat));
for dd = 1:nd
    itest = DWN.divenum == dives(dd) & qcdown ~=4 & ~isnan(DWN.optode_oxygen);
    data = DWN.optode_oxygen(itest);
    press = DWN.press(itest);
    [qc] = argoFiltersTest('spike',data,press,'abovebelow','oxygen');
    qcdown(itest) = qc;
end
dives = unique(UP.divenum); % DWN
nd = numel(dives);
qcup = ones(size(UP.lat));
for dd = 1:nd
    itest = UP.divenum == dives(dd) & qcup ~=4 & ~isnan(UP.optode_oxygen);
    data = UP.optode_oxygen(itest);
    press = UP.press(itest);
    [qc] = argoFiltersTest('spike',data,press,'abovebelow','oxygen');
    qcup(itest) = qc;
end

% 4. Assign variable qc flags based on above tests
% optode oxygen
flags = {'qcflag_visual_general','qcflag_visual_aanderaa','qcflag_range_oxygenAA','qcflag_range_phaseAA','qcflag_biofoul','temp_qcflag','salin_qcflag'};
iset3 = strcmp(flags,'temp_qcflag') | strcmp(flags,'salin_qcflag');
var_flags = QC_UP{:,flags};
var_4mask = var_flags == 4; var_4mask(:,~iset3) = false;
var_flags(var_4mask) = 2;
QC_UP.optode_oxygen_qcflag = max(var_flags,[],2,'omitnan');
var_flags = QC_DWN{:,flags};
var_4mask = var_flags == 4; var_4mask(:,~iset3) = false;
var_flags(var_4mask) = 2;
QC_DWN.optode_oxygen_qcflag = max(var_flags,[],2,'omitnan');
% optode phase
flags = {'qcflag_visual_general','qcflag_visual_aanderaa','qcflag_range_phaseAA','qcflag_biofoul'};
var_flags = QC_UP{:,flags};
QC_UP.optode_phase_qcflag = max(var_flags,[],2,'omitnan');
var_flags = QC_DWN{:,flags};
QC_DWN.optode_phase_qcflag = max(var_flags,[],2,'omitnan');
% optode temperature
flags = {'qcflag_visual_general','qcflag_visual_aanderaa','qcflag_range_tempAA'};
var_flags = QC_UP{:,flags};
QC_UP.optode_temp_qcflag = max(var_flags,[],2,'omitnan');
var_flags = QC_DWN{:,flags};
QC_DWN.optode_temp_qcflag = max(var_flags,[],2,'omitnan');

% -------------------------------------------------------------------------
% Winkler Calibration
% settings
calfilepath = 'glider-kit/calibrationFiles/oxy_cal_19-Jul-2023.mat';
method = 'gain';
% patch to remove large surface spikes from calibration
if strcmp(missionname,'sg512_1') || strcmp(missionname,'sg626_1')
    upspike = UP.vmdepth < 15 & (UP.optode_oxygen < 166 | UP.optode_oxygen > 188);
    downspike = DWN.vmdepth < 15 & (DWN.optode_oxygen < 166 | DWN.optode_oxygen > 188);
elseif strcmp(missionname,'sg512_4') 
    upspike = UP.vmdepth < 15 & (UP.optode_oxygen < 196 | UP.optode_oxygen > 201);
    downspike = DWN.vmdepth < 15 & (DWN.optode_oxygen < 196 | DWN.optode_oxygen > 201);
else
    upspike = isnan(UP.optode_oxygen); % filler mask
    downspike = isnan(DWN.optode_oxygen);  % filler mask
end

gliderTable = [DWN(QC_DWN.optode_oxygen_qcflag <3 & ~downspike,{'optode_oxygen','divenum','dateUTC','vmdepth','lat','lon','sigmath0'});...
    UP(QC_UP.optode_oxygen_qcflag <3 & ~upspike,{'optode_oxygen','divenum','dateUTC','vmdepth','lat','lon','sigmath0'})]; % DWN & UP

[~,oxygencalInfo] = sg_oxygenWinklerCalibration(gliderTable,'optode_oxygen',calfilepath,method,missionname);
calib.optode_oxygen_cal = oxygencalInfo;
DWN.optode_oxygen_cal = DWN.optode_oxygen.*oxygencalInfo.calibrationCoefficient;
UP.optode_oxygen_cal = UP.optode_oxygen.*oxygencalInfo.calibrationCoefficient;
% save figure
figdir = 'figures/aaoxy_cal/';
if ~isfolder(figdir); mkdir(figdir); end
exportgraphics(gcf,[figdir missionname '.png'],'Resolution',300);
close

% -------------------------------------------------------------------------
% Bittig et al. (2014) inverse filter correction
% Optode instrument response correction
M = [DWN(:,{'optode_oxygen_cal','divenum','vmtime','press','direction'});...
    UP(:,{'optode_oxygen_cal','divenum','vmtime','press','direction'})]; % DWN & UP
M.optode_oxygen_qcflag = [QC_DWN.optode_oxygen_qcflag; QC_UP.optode_oxygen_qcflag]; % use flag for bittig correction
M.optode_oxygenc_cal = nan.*M.optode_oxygen_cal;

% settings
timestep = 1;
movingpoints = 60;

% optimize tauoptode from 1:100 seconds (DEFAULT 30 seconds)
%sg_optimizeTauoptode(M); % Add in future version


% loop through dives
dives = unique(M.divenum);
dives(isnan(dives)) = [];
ndives = numel(dives);
for dd = 1:ndives

    % flagged & unflagged
    igoodbad = M.divenum == dives(dd);
    dT = M(igoodbad,:);
    M.optode_oxygenc_cal(igoodbad) = optodeBittigCorrection(dT.vmtime,dT.optode_oxygen_cal,tauoptode,timestep,movingpoints,dT.press);

    % quality data ONLY to overwrite above
    igood = M.divenum == dives(dd) & M.optode_oxygen_qcflag<3;
    dT = M(igood,:);
    M.optode_oxygenc_cal(igood) = optodeBittigCorrection(dT.vmtime,dT.optode_oxygen_cal,tauoptode,timestep,movingpoints,dT.press);

end
DWN.optode_oxygenc_cal = M.optode_oxygenc_cal(M.direction == "down");
UP.optode_oxygenc_cal = M.optode_oxygenc_cal(M.direction == "up");

% save processed data
save(['data/qcdata_L2/' missionname '_L2'],'UP','DWN','calib','QC_DWN','QC_UP');


end