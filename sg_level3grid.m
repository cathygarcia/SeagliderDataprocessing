% SG_LEVEL3GRID
% -------------------------------------------------------------------------
% Applies qc flags to processed data, and bins to 2 m grid
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
% function sg_level3grid(missionname)
%
% PROCESSING STEPS
% -------------------------------------------------------------------------
% 1. Load qc data file
% 2. Apply quality control flags
% 3. Grids to 2 m
% 4. Reshapes to [n depths, n dives]
% 5. Grabs metadata (dcm, mld)
% 6. Save tables to qcdata folder 'data/griddata_L3/missionname_L3.mat'
%
% -------------------------------------------------------------------------
% Author Catherine Garcia -- cathy.garcia@hawaii.edu -- Version 04 Aug 2023
%
% -------------------------------------------------------------------------

function sg_level3grid(missionname)

filename = ['data/qcdata_L2/' missionname '_L2.mat'] ; %file path to mission process data
load(filename,'DWN','UP','QC_UP','QC_DWN');

% Apply flags
[DWN,UP] = sg_removeFlagvars('temp_qcflag',{'tempc'},DWN,UP,QC_DWN,QC_UP);% temp flag
[DWN,UP] = sg_removeFlagvars('cond_qcflag',{'condc'},DWN,UP,QC_DWN,QC_UP); % condc flag
[DWN,UP] = sg_removeFlagvars('salin_qcflag',{'sp','sa'},DWN,UP,QC_DWN,QC_UP);% salin flag
[DWN,UP] = sg_removeFlagvars('sigma_qcflag',{'sigmath0'},DWN,UP,QC_DWN,QC_UP);% sigmath0 glag
[DWN,UP] = sg_removeFlagvars('chl1_qcflag',DWN.Properties.VariableNames(contains(DWN.Properties.VariableNames,'chl1')),DWN,UP,QC_DWN,QC_UP);% chl1 flag
[DWN,UP] = sg_removeFlagvars('chl2_qcflag',DWN.Properties.VariableNames(contains(DWN.Properties.VariableNames,'chl2')),DWN,UP,QC_DWN,QC_UP);% chl2 flag
[DWN,UP] = sg_removeFlagvars('chl_cal_qcflag',DWN.Properties.VariableNames(contains(DWN.Properties.VariableNames,'chl_cal')),DWN,UP,QC_DWN,QC_UP);% chl_cal flag
[DWN,UP] = sg_removeFlagvars('optode_oxygen_qcflag',DWN.Properties.VariableNames(contains(DWN.Properties.VariableNames,'optode_oxygen')),DWN,UP,QC_DWN,QC_UP);% optode oxy flag
[DWN,UP] = sg_removeFlagvars('optode_temp_qcflag',{'optode_temp'},DWN,UP,QC_DWN,QC_UP);% optode temp glag
[DWN,UP] = sg_removeFlagvars('oxygenSBE_qcflag',{'oxygenSBE','oxygenSBE_cal'},DWN,UP,QC_DWN,QC_UP);% sbe43 oxygen flag
[DWN,UP] = sg_removeFlagvars('bbp470_qcflag',DWN.Properties.VariableNames(contains(DWN.Properties.VariableNames,'bbp470') & ~contains(DWN.Properties.VariableNames,'spikeflag')),DWN,UP,QC_DWN,QC_UP);% bbp470 flag
[DWN,UP] = sg_removeFlagvars('bbp700_qcflag',DWN.Properties.VariableNames(contains(DWN.Properties.VariableNames,'bbp700') & ~contains(DWN.Properties.VariableNames,'spikeflag')),DWN,UP,QC_DWN,QC_UP);% bbp700 flag
[DWN,UP] = sg_removeFlagvars('bbp650_qcflag',DWN.Properties.VariableNames(contains(DWN.Properties.VariableNames,'bbp650') & ~contains(DWN.Properties.VariableNames,'spikeflag')),DWN,UP,QC_DWN,QC_UP);% bbp650 flag
[DWN,UP] = sg_removeFlagvars('bbp660_qcflag',DWN.Properties.VariableNames(contains(DWN.Properties.VariableNames,'bbp660') & ~contains(DWN.Properties.VariableNames,'spikeflag')),DWN,UP,QC_DWN,QC_UP);%
[DWN,UP] = sg_removeFlagvars('cdom_qcflag',{'cdom'},DWN,UP,QC_DWN,QC_UP);%

% -------------------------------------------------------------------------
% save qc'd data only
filename = ['data/griddata_L3/' missionname '_L3.mat'] ; 
save(filename,"UP","DWN");


% -------------------------------------------------------------------------
% Grid data every 2m *** add iso grid
depth = 2:2:600;
outVarNames = {'divenum','dateUTC','dateHST','vmtime','lat','lon','tempc','condc','sp','sa','ct','sigmath0','oxygenSBE','oxygenSBE_cal',...
    'optode_oxygen','optode_phase','optode_temp','optode_oxygen_cal','optode_oxygenc_cal',...
    'optode_oxygen_uM','optode_oxygen_cal_uM','optode_oxygenc_cal_uM',...
    'beta470','beta700','beta650','beta660',...
    'bbp470','bbp700','bbp650','bbp660',...,...
    'bbp470_despiked','bbp700_despiked','bbp650_despiked','bbp660_despiked',...,...
    'bbp470c','bbp700c','bbp650c','bbp660c',...
    'bbp470c_despiked','bbp700c_despiked','bbp650c_despiked','bbp660c_despiked',...
    'cdom','chl1','chl1_cal','chl2','chl2_cal','chl_cal',...,...
    'chl1_cal_despiked','chl1_despiked','chl2_cal_despiked','chl_cal_despiked',...
    'pcmDiff'};
[TF,~] = ismember(outVarNames,DWN.Properties.VariableNames); % only keep valid names
outVarNames = outVarNames(TF);
diveRange = [min(DWN.divenum), max(DWN.divenum)];


% old version uses naninterp1 to linearly interp on z_grid
[UG,DG] = sg_grid(filename,depth,'gridVar','vmdepth','diveRange',diveRange(1):diveRange(2),...
   'outVars',outVarNames);


% -------------------------------------------------------------------------
% Put data on a regular grid (dives,depth)
varNamesOut = DG.Properties.VariableNames;
varNamesIn = DG.Properties.VariableNames;
varUnits = {};
dives = unique([UP.divenum(~isnan(UP.divenum));DWN.divenum(~isnan(DWN.divenum))]);

[sgd,sgu,dived,diveu] = putOnGrid(DG,UG,varNamesIn,varNamesOut,varUnits,depth,dives);

%patch for nan dates (e.g. 148_11 some UP dives only have data below 600m
%and dates are missed by grid);
inan = find(isnan(dived.dateHST)); % down HST
for ii = 1:numel(inan)
    ix = DWN.divenum == dived.dive(inan(ii));
    dived.dateHST(inan(ii)) = mean(DWN.dateHST(ix),'omitnan');
    dived.dateUTC(inan(ii)) = mean(DWN.dateUTC(ix),'omitnan');
end
inan = find(isnan(diveu.dateHST)); % up HST
for ii = 1:numel(inan)
    ix = UP.divenum == diveu.dive(inan(ii));
    diveu.dateHST(inan(ii)) = mean(UP.dateHST(ix),'omitnan');
    diveu.dateUTC(inan(ii)) = mean(UP.dateUTC(ix),'omitnan');
end
% -------------------------------------------------------------------------
% Get additional metadata for dataanalysis (DCM, surface layer, mixed layer)
[sgd,dived] = sg_metadata(sgd,dived);
[sgu,diveu] = sg_metadata(sgu,diveu);
% fill in down
var = 'mld001'; dived.(var)(isnan(dived.(var))) = diveu.(var)(isnan(dived.(var)));
var = 'mld003'; dived.(var)(isnan(dived.(var))) = diveu.(var)(isnan(dived.(var)));
var = 'mld0125'; dived.(var)(isnan(dived.(var))) = diveu.(var)(isnan(dived.(var)));
var = 'fcm'; dived.(var)(isnan(dived.(var))) = diveu.(var)(isnan(dived.(var)));
var = 'zcm'; dived.(var)(isnan(dived.(var))) = diveu.(var)(isnan(dived.(var)));
var = 'sigcm'; dived.(var)(isnan(dived.(var))) = diveu.(var)(isnan(dived.(var)));
% fill in up
var = 'mld001'; diveu.(var)(isnan(diveu.(var))) = dived.(var)(isnan(diveu.(var)));
var = 'mld003'; diveu.(var)(isnan(diveu.(var))) = dived.(var)(isnan(diveu.(var)));
var = 'mld0125'; diveu.(var)(isnan(diveu.(var))) = dived.(var)(isnan(diveu.(var)));
var = 'fcm'; diveu.(var)(isnan(diveu.(var))) = dived.(var)(isnan(diveu.(var)));
var = 'zcm'; diveu.(var)(isnan(diveu.(var))) = dived.(var)(isnan(diveu.(var)));
var = 'sigcm'; diveu.(var)(isnan(diveu.(var))) = dived.(var)(isnan(diveu.(var)));

% -------------------------------------------------------------------------
% Compute characteristics on isopycnal levels (average sigma at depth bins)
% compute average sigma at depth bins
sig_grid_d = mean(sgd.sigmath0,2,'omitnan');
sig_grid_d = sort(sig_grid_d);
sig_grid_u = mean(sgu.sigmath0,2,'omitnan');
sig_grid_u = sort(sig_grid_u);
% initialize table of isopycnal data
isod = sgd;
isod.sig = sig_grid_d;
isou = sgu;
isou.sig = sig_grid_u;
varNames = sgd.Properties.VariableNames(2:end);
nd = numel(dived.dive);
for vv = 1:numel(varNames)
    isod.(varNames{vv}) = NaN*isod.(varNames{vv});
    isou.(varNames{vv}) = NaN*isou.(varNames{vv});
end

for ii = 1:nd
    ind_nan = isnan(sgd.sigmath0(:,ii)); % DWN
    if sum(~ind_nan) > 1
        for vv = 1:numel(varNames)
            isod.(varNames{vv})(:,ii) = interp1(sgd.sigmath0(~ind_nan,ii),sgd.(varNames{vv})(~ind_nan,ii),sig_grid_d);
        end
        clear ind_nan
    end
    ind_nan = isnan(sgu.sigmath0(:,ii)); % UP
    if sum(~ind_nan) > 1
        for vv = 1:numel(varNames)
            isou.(varNames{vv})(:,ii) = interp1(sgu.sigmath0(~ind_nan,ii),sgu.(varNames{vv})(~ind_nan,ii),sig_grid_u);
        end
        clear ind_nan
    end
end


% save processed data
save(['data/griddata_L3/' missionname '_L3'],'DWN','UP','UG','DG','sgd','sgu','isod','isou','dived','diveu');


% -------------------------------------------------------------------------
% function to remove flagged data & set to NaN
    function [DWN,UP] = sg_removeFlagvars(flag,flagvars,DWN,UP,QC_DWN,QC_UP)
        if sum(strcmp(QC_DWN.Properties.VariableNames,flag)) == 1 && ~isempty(flagvars)% is flag exists, flagvars should exist
            DWN(QC_DWN.(flag) > 2,flagvars) = table(nan);
            UP(QC_UP.(flag) > 2,flagvars) = table(nan);
        end
    end

end