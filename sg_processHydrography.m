% SG_PROCESSHYDROGRAPHY
% -------------------------------------------------------------------------
% Qality control for Sea-Bird CT Sail temp, cond, salin, & sigma; data
% reprocessed with 1st order lag without flagged data
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
% function sg_processHydrography(missionname)
%
% PROCESSING STEPS
% -------------------------------------------------------------------------
% 1. Load raw data file "DWN0" & "UP0:
% 2. Apply quality control flags
%   - visual qc (after all qc, remaining bad data outliers)
%   - pitch & density inversion qc (pitch for large sections, dens inv for
%   single point inversions)
%   - range qc (glogal argo test, regional HOT test)
% 4. Reprocess tempc, condc, salin all, sigma
% 5. Save tables to qcdata folder 'data/qcdata/missionname_L2.mat'

% MATLAB CODE NEEDED
% -------------------------------------------------------------------------
% gsw toolbox
%
% -------------------------------------------------------------------------
% Author Catherine Garcia -- cathy.garcia@hawaii.edu -- Version 04 Aug 2023
%
% -------------------------------------------------------------------------

function sg_processHydrography(missionname)

filename = ['data/rawdata_L1/' missionname '_L1.mat'] ; %file path to mission
load(filename,'DWN0','UP0','calib');


% New variables to save for processed data
varlist = {'direction','index','divenum','lat','lon','dateUTC','dateHST','press','vmdepth','vmtime','pitch'};
DWN = DWN0(:,varlist);
UP = UP0(:,varlist);
QC_DWN = DWN0(:,'index');
QC_UP = UP0(:,'index');

% -------------------------------------------------------------------------
% QC Analysis

% 0. VISUAL INSPECTION: remove obviously bad dives where sensors are damaged or bad
% These visual flags are added to code after first round of processing, and
% then kept for all future runs

% flag entire dive for every sensor
[qcup,qcdown] = sg_greylistQC(missionname,UP0,DWN0,'general');
QC_UP.qcflag_visual_general = qcup;
QC_DWN.qcflag_visual_general = qcdown;

% flag sbe ctsail specifically
[qcup,qcdown] = sg_greylistQC(missionname,UP0,DWN0,'ctsail');
QC_UP.qcflag_visual_ctsail = qcup;
QC_DWN.qcflag_visual_ctsail = qcdown;

% 1. DENSITY INVERSIONS:
% flag data where pitch angles opposite expected
% flag argo density inversion test

% Sharp changes in pitch lead to density invervions
QC_UP.qcflag_pitch = ones(size(UP0.pitch));
QC_UP.qcflag_pitch(UP0.pitch < 5) = 4; % flag negative pitch angles
QC_DWN.qcflag_pitch = ones(size(DWN0.pitch));
QC_DWN.qcflag_pitch(DWN0.pitch > -5 & DWN0.vmdepth > 75) = 4; % flag positive pitch angles, below 75m

% Density inversions also useful on top of pitch changes for bottom of
% dives, uses argo density inverion test
[qcup,qcdown] = sg_qcGradientTest(UP0,DWN0);
QC_UP.qcflag_densinv = qcup;
QC_DWN.qcflag_densinv = qcdown;

% 2. Flag all variables (and related ones) as fail for out of acceptable
% range
% global range test from argo float QC/QC DAC
% preallocate
QC_UP.qcflag_range_temp = ones(height(UP),1);
QC_UP.qcflag_range_salin = ones(height(UP),1);
QC_UP.qcflag_range_cond = ones(height(UP),1);
QC_UP.qcflag_range_sigma = ones(height(UP),1);
QC_UP.qcflag_range_press = ones(height(UP),1);
QC_DWN.qcflag_range_temp = ones(height(DWN),1);
QC_DWN.qcflag_range_salin = ones(height(DWN),1);
QC_DWN.qcflag_range_cond = ones(height(DWN),1);
QC_DWN.qcflag_range_sigma = ones(height(DWN),1);
QC_DWN.qcflag_range_press = ones(height(DWN),1);
% regional ranges (run first to not overwrite global range)
QC_UP.qcflag_range_temp(UP0.tempc < -2.5 | UP0.tempc > 40) = 3; % temperature global range (deg C)
QC_UP.qcflag_range_salin(UP0.sp < 2 | UP0.sp > 41) = 4; % salinity global range (PSU)
QC_UP.qcflag_range_cond(UP0.condc < 0 | UP0.condc > 85) = 3; % conductivity global range (mS/cm)
QC_UP.qcflag_range_sigma(UP0.sigmath0 < 0 | UP0.sigmath0 > 60) = 3; % pot. dens. anomaly (kg m-3)
QC_UP.qcflag_range_press(UP0.press < -2.4 | UP0.press > 2000) = 3; % pressure (db)
QC_DWN.qcflag_range_temp(DWN0.tempc < -2.5 | DWN0.tempc > 40) = 3; % temperature global range (deg C)
QC_DWN.qcflag_range_salin(DWN0.sp < 2 | DWN0.sp > 41) = 3; % salinity global range (PSU)
QC_DWN.qcflag_range_cond(DWN0.condc < 0 | DWN0.condc > 85) = 3; % conductivity global range (mS/cm)
QC_DWN.qcflag_range_sigma(DWN0.sigmath0 < 0 | DWN0.sigmath0 > 60) = 3; % pot. dens. anomaly (kg m-3)
QC_DWN.qcflag_range_press(DWN0.press < -2.4 | DWN0.press > 2000) = 3; % pressure (db)
% global ranges (overwrites regional range to for large outliers)
QC_UP.qcflag_range_temp(UP0.tempc < -2.5 | UP0.tempc > 40) = 4; % temperature global range (deg C)
QC_UP.qcflag_range_salin(UP0.sp < 2 | UP0.sp > 41) = 4; % salinity global range (PSU)
QC_UP.qcflag_range_cond(UP0.condc < 0 | UP0.condc > 85) = 4; % conductivity global range (mS/cm)
QC_UP.qcflag_range_sigma(UP0.sigmath0 < 0 | UP0.sigmath0 > 60) = 4; % pot. dens. anomaly (kg m-3)
QC_UP.qcflag_range_press(UP0.press < -2.4 | UP0.press > 2000) = 4; % pressure (db)
QC_DWN.qcflag_range_temp(DWN0.tempc < -2.5 | DWN0.tempc > 40) = 4; % temperature global range (deg C)
QC_DWN.qcflag_range_salin(DWN0.sp < 2 | DWN0.sp > 41) = 4; % salinity global range (PSU)
QC_DWN.qcflag_range_cond(DWN0.condc < 0 | DWN0.condc > 85) = 4; % conductivity global range (mS/cm)
QC_DWN.qcflag_range_sigma(DWN0.sigmath0 < 0 | DWN0.sigmath0 > 60) = 4; % pot. dens. anomaly (kg m-3)
QC_DWN.qcflag_range_press(DWN0.press < -2.4 | DWN0.press > 2000) = 4; % pressure (db)

% 3. Assign variable qc flags based on above tests
% temperature
flags = {'qcflag_visual_general','qcflag_visual_ctsail','qcflag_pitch','qcflag_densinv','qcflag_range_temp'};
iset3 = strcmp(flags,'qcflag_pitch') | strcmp(flags,'qcflag_pitch');
var_flags = QC_UP{:,flags};
var_4mask = var_flags == 4; var_4mask(:,~iset3) = false;
var_flags(var_4mask) = 3;
QC_UP.temp_qcflag = max(var_flags,[],2,'omitnan');
var_flags = QC_DWN{:,flags};
var_4mask = var_flags == 4; var_4mask(:,~iset3) = false;
var_flags(var_4mask) = 3;
QC_DWN.temp_qcflag = max(var_flags,[],2,'omitnan');
% conductivity
flags = {'qcflag_visual_general','qcflag_visual_ctsail','qcflag_pitch','qcflag_densinv','qcflag_range_cond'};
iset3 = strcmp(flags,'qcflag_pitch') | strcmp(flags,'qcflag_pitch');
var_flags = QC_UP{:,flags};
var_4mask = var_flags == 4; var_4mask(:,~iset3) = false;
var_flags(var_4mask) = 3;
QC_UP.cond_qcflag = max(var_flags,[],2,'omitnan');
var_flags = QC_DWN{:,flags};
var_4mask = var_flags == 4; var_4mask(:,~iset3) = false;
var_flags(var_4mask) = 3;
QC_DWN.cond_qcflag = max(var_flags,[],2,'omitnan');
% salinity
flags = {'qcflag_visual_general','qcflag_visual_ctsail','qcflag_pitch','qcflag_densinv','qcflag_range_salin','temp_qcflag','cond_qcflag'};
var_flags = QC_UP{:,flags};
QC_UP.salin_qcflag = max(var_flags,[],2,'omitnan');
var_flags = QC_DWN{:,flags};
QC_DWN.salin_qcflag = max(var_flags,[],2,'omitnan');
% potential density anomaly
flags = {'qcflag_visual_general','qcflag_visual_ctsail','qcflag_pitch','qcflag_densinv','qcflag_range_sigma','salin_qcflag','temp_qcflag','cond_qcflag'};
var_flags = QC_UP{:,flags};
QC_UP.sigma_qcflag = max(var_flags,[],2,'omitnan');
var_flags = QC_DWN{:,flags};
QC_DWN.sigma_qcflag = max(var_flags,[],2,'omitnan');

% -------------------------------------------------------------------------
% Reprocess hydrography & temp,cond 1st order lag without flagged data
M = [DWN0;UP0]; % combine to run dives together
MQC = [QC_DWN;QC_UP];
dives = unique(M.divenum); dives(isnan(dives)) = [];
ndives = numel(dives);
for dd = 1:ndives
    % index dives and good data
    idive = M.divenum == dives(dd);
    % temp lag-----------------------------------------
    tempqc = idive & MQC.temp_qcflag < 4;
    tautemp = 0.6;
    M.tempc(tempqc) = M.temp(tempqc) + tautemp*gradient(M.temp(tempqc), M.vmtime(tempqc));
    % cond lag-----------------------------------------
    condqc = idive & MQC.cond_qcflag < 4;
    taucond = 1./( 0.01 + (M.spdg(condqc)/10).*(M.spdg(condqc)/10)/5.97 );
    M.condc(condqc) = M.cond(condqc) + taucond.*gradient(M.cond(condqc), M.vmtime(condqc));
end
M.sp = gsw_SP_from_C(M.condc,M.tempc,M.press);
M.sa = gsw_SA_from_SP(M.sp,M.press,M.lon,M.lat); % Absolute salinity
M.ct = gsw_CT_from_t(M.sa,M.tempc,M.press); % Conservative temperature
M.sigmath0 = gsw_sigma0(M.sa,M.ct); % replaced 1. salin with sa and 2. tempc with ct Cathy 08June2022

% vars to add to processed data tables
varsAppend = {'tempc','condc','sp','sa','ct','sigmath0'};
vdown = M.direction == "down";
vup = M.direction == "up";
DWN = [DWN, M(vdown,varsAppend)];
UP = [UP, M(vup,varsAppend)];


% save processed data
save(['data/qcdata_L2/' missionname '_L2'],'UP','DWN','calib','QC_DWN','QC_UP');


end








