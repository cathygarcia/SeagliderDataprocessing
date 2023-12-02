% MAKENetCDFFILE
% -------------------------------------------------------------------------
% Wrapper code to write finalized data to shareable .nc file for public use
%
% INPUTS
% -------------------------------------------------------------------------
% REQUIRED INPUTS:
% missionname:  glider mission file path to raw .eng files
%
% OPTIONAL INPUTS:
% divesToRemove:   optode time delay for phase correction (default=30seconds)
%
% OUTPUTS
% -------------------------------------------------------------------------
% UP:           output for ascending data in
% DWN:          output for descending data
% outputs also saved to rawdata folder
%
% USAGE
% -------------------------------------------------------------------------
% function makeNetCDFfile(missionname)
%
% PROCESSING STEPS
% -------------------------------------------------------------------------
%  Write variables to nc file
%
% -------------------------------------------------------------------------
% Author Catherine Garcia -- cathy.garcia@hawaii.edu -- Version 13 Jun 2023
%
% -------------------------------------------------------------------------

function makeNetCDFfile(missionname)


% Make NetCDF file
savepath = 'data/ncfiles/';


%% ------------ MISSION FILE ----------------------------------------------------
% create nc file
ncfile = [savepath missionname '.nc'];
if exist(ncfile,'file')==2; delete(ncfile); end % delete old file

% ************ Metadata ***********
datafile = ['data/qcdata_L2/' missionname '_L2.mat'] ; %file path to mission
load(datafile,'UP','DWN','QC_UP','QC_DWN','calib');
M = [DWN;UP];
MQC = [QC_DWN;QC_UP];
[~,idx] = sort(M.dateHST,'ascend'); % sort by date
M = M(idx,:);
MQC = MQC(idx,:);
% make new direction column for ncwrite, difficulty with categorical array
M.directionlogical = nan.*M.lat;
M.directionlogical(M.direction == "down") = 0;
M.directionlogical(M.direction == "up") = 1;
[rows,~] = size(M);

% load sensor file
sensordata = readtable('data/table_ofSensors_bymission.xlsx'); % sensor information
sensorrow = find(strcmp(sensordata.Mission,missionname)); % mission row

% Variable names on all missions
sgvar = {'lat',... % data variable name
    'lon',...
    'dateUTC',...
    'divenum',...
    'press',...
    'vmdepth',...
    'directionlogical'};
ncvars = {'lat',... % name for ncfile variable
    'lon',...
    'datetime',...
    'divenum',...
    'press',...
    'depth',...
    'direction'};
sgunits = {'decimal degrees N [-90, 90]',...% variable units
    'decimal degrees E [-180, 180]',...
    'yyyy-MM-ddTHH:mm:ss',...
    'dive number',...
    'decibars',...
    'meters',...
    '1 = ascent, 0 = descent/down'};
sglongname = {'Latitude',...% long descriptor name
    'Longitude',...
    'Datetime',...
    'Dive number',...
    'Pressure',...
    'Depth',...
    'Profile direction'};

for ii = 1:numel(ncvars)
    if sum(strcmp(M.Properties.VariableNames,sgvar{ii})) == 1 % case data var exists
        var = M.(sgvar{ii});
    else
        continue
    end
    varname = ncvars{ii};
    varunits = sgunits{ii};
    varlongname = sglongname{ii};
    nccreate(ncfile,varname,'Dimensions',{varname,rows});
    ncwrite(ncfile,varname,var);
    ncwriteatt(ncfile,varname,'long name',varlongname);
    ncwriteatt(ncfile,varname,'units',varunits);
    % add pressure sensor
    if strcmp(sgvar{ii},'press')
        varsensor = sensordata{sensorrow,sgvar{ii}};
        ncwriteatt(ncfile,varname,'sensor',varsensor{:});
    end
end

% ************ Main data ***********
% data variable name
sgvar = {'pitch',... % data variable name
    'tempc',...
    'condc',...
    'sa',...
    'sigmath0',...
    'optode_oxygen_uM',...
    'optode_oxygenc_cal_uM',...
    'chl',...
    'chl_cal',...
    'bbp470',...
    'bbp470c',...
    'bbp470_spikeflag',...
    'bbp650',...
    'bbp650c',...
    'bbp650_spikeflag',...
    'bbp660',...
    'bbp660c',...
    'bbp660_spikeflag',...
    'bbp700',...
    'bbp700c',...
    'bbp700_spikeflag',...
    'temp_qcflag',... % variable qc flag
    'cond_qcflag',...
    'salin_qcflag',...
    'sigma_qcflag',...
    'optode_oxygen_qcflag',...
    'chl_cal_qcflag',...
    'bbp470_qcflag',...
    'bbp650_qcflag',...
    'bbp660_qcflag',...
    'bbp700_qcflag'};
ncvars = {'pitch',...  % name for ncfile variable
    'tempc',...
    'condc',...
    'salin',...
    'sigma',...
    'oxygen_raw',...
    'oxygen',...
    'chla_raw',...
    'chla',...
    'bbp470',...
    'bbp470c',...
    'bbp470_spikeflag',...
    'bbp650',...
    'bbp650c',...
    'bbp650_spikeflag',...
    'bbp660',...
    'bbp660c',...
    'bbp660_spikeflag',...
    'bbp700',...
    'bbp700c',...
    'bbp700_spikeflag',...
    'temp_qcflag',... % variable qc flag
    'cond_qcflag',...
    'salin_qcflag',...
    'sigma_qcflag',...
    'oxygen_qcflag',...
    'chla_qcflag',...
    'bbp470_qcflag',...
    'bbp650_qcflag',...
    'bbp660_qcflag',...
    'bbp700_qcflag'};
sgunits = {'degrees',... % variable units
    '° Celcius',...
    'mS cm-1',...
    'g kg-1',...
    'kg m-3',...
    'µmol L-1',...
    'µmol L-1',...
    'mg m-3',...
    'mg m-3',...
    'm-1',...
    'm-1',...
    '0 for no spike, 1 for spike',...
    'm-1',...
    'm-1',...
    '0 for no spike, 1 for spike',...
    'm-1',...
    'm-1',...
    '0 for no spike, 1 for spike',...
    'm-1',...
    'm-1',...
    '0 for no spike, 1 for spike',...
    'QC pass = 1','QC probably bad = 3','QC fail = 4',...
    'QC pass = 1','QC probably bad = 3','QC fail = 4',...
    'QC pass = 1','QC probably bad = 3','QC fail = 4',...
    'QC pass = 1','QC probably bad = 3','QC fail = 4',...
    'QC pass = 1','QC probably bad = 3','QC fail = 4',...
    'QC pass = 1','QC probably bad = 3','QC fail = 4',...
    'QC pass = 1','QC probably bad = 3','QC fail = 4',...
    'QC pass = 1','QC probably bad = 3','QC fail = 4',...
    'QC pass = 1','QC probably bad = 3','QC fail = 4',...
    'QC pass = 1','QC probably bad = 3','QC fail = 4'};
sglongname = {'Pitch angle',... % long descriptor name
    'Temperature',...
    'Conductivity',...
    'Absolute salinity',...
    'Potential density anomaly at 0 dbar',...
    'Raw dissolved oxygen concentration',...
    'Dissolved oxygen concentration',...
    'Raw chlorophyll-a concentration',...
    'Chlorophyll-a concentration',...
    'Particulate backscattering coefficient (λ = 470 nm)',...
    'Particulate backscattering coefficient (λ = 470 nm) Corrected',...
    'bbp at 470 nm spike flag',...
    'Particulate backscattering coefficient (λ = 650 nm)',...
    'Particulate backscattering coefficient (λ = 650 nm) Corrected',...
    'bbp at 650 nm spike flag',...
    'Particulate backscattering coefficient (λ = 660 nm)',...
    'Particulate backscattering coefficient (λ = 660 nm) Corrected',...
    'bbp at 660 nm spike flag',...
    'Particulate backscattering coefficient (λ = 700 nm)',...
    'Particulate backscattering coefficient (λ = 700 nm) Corrected',...
    'bbp at 700 nm spike flag',...
    'temperature QC flag',...
    'conductivity QC flag',...
    'salinity QC flag',...
    'sigma QC flag',...
    'oxygen QC flag',...
    'chla QC flag',...
    'bbp470 QC flag',...
    'bbp650 QC flag',...
    'bbp660 QC flag',...
    'bbp700 QC flag'};
% loop over vars
for ii = 1:numel(ncvars)
    if sum(strcmp(M.Properties.VariableNames,sgvar{ii})) == 1 % case var exists
        var = M.(sgvar{ii});
    elseif sum(strcmp(MQC.Properties.VariableNames,sgvar{ii})) == 1 % case  qcflag var exists
        var = MQC.(sgvar{ii});
    elseif strcmp(sgvar{ii},'chl') % case chl : var depends on chl1 or chl2 reported
        if strcmp(calib.chl_cal_var,'chl1') % report chl1
               var = M.chl1;
        elseif strcmp(calib.chl_cal_var,'chl2') % report chl2
              var = M.chl2;
        end
    else
        continue
    end
    % set logical to double
    if islogical(var)
        var = double(var);
    end

    varname = ncvars{ii};
    varunits = sgunits{ii};
    varlongname = sglongname{ii};
    nccreate(ncfile,varname,'Dimensions',{varname,rows});
    ncwrite(ncfile,varname,var);
    ncwriteatt(ncfile,varname,'long name',varlongname);
    ncwriteatt(ncfile,varname,'units',varunits);

    % add sensor info if present
    if sum(strcmp(sensordata.Properties.VariableNames,sgvar{ii})) == 1
        varsensor = sensordata{sensorrow,sgvar{ii}};
        ncwriteatt(ncfile,varname,'sensor',varsensor{:});
        % chl sensor depends on chl1, or chl2 reported
        if strcmp(sgvar{ii},'chl')
        elseif strcmp(calib.chl_cal_var,'chl1') % chl1, report wlbb2fl (bbp470,bbp700,chl)
            varsensor = sensordata{sensorrow,'bbp470'};
            ncwriteatt(ncfile,varname,'sensor',varsensor{:});
            chl_var = 'chl1';
        elseif strcmp(calib.chl_cal_var,'chl2') % chl2, report wlbbfl2 (bbp650 or bbp 660,chl,cdom)
            varsensor = sensordata{sensorrow,'cdom'};
            ncwriteatt(ncfile,varname,'sensor',varsensor{:});
            chl_var = 'chl2';
        end
    end
end

%% ------------ MISSION EXPANDED QC FLAG FILE ----------------------------------------------------
% create nc file
ncqcfile = [savepath missionname '_QCflags.nc'];
if exist(ncqcfile,'file')==2; delete(ncqcfile); end % delete old file

% ************ QC flags ***********
% data variable name
% data variable name
sgvar = {'temp_qcflag',... % variable qc flag
    'cond_qcflag',...
    'salin_qcflag',...
    'sigma_qcflag',...
    'optode_oxygen_qcflag',...
    'chl_cal_qcflag',...
    'bbp470_qcflag',...
    'bbp650_qcflag',...
    'bbp660_qcflag',...
    'bbp700_qcflag',...
    'qcflag_pitch',...
    'qcflag_densinv',...
    'qcflag_biofoul',...
    'qcflag_biofoulOptics',...
    'qcflag_visual_general',...
    'qcflag_visual_CTflag',...
    'qcflag_visual_aanderaa',...
    ['qcflag_visual_' chl_var],...
    'qcflag_visual_bbp470',...
    'qcflag_visual_bbp650',...
    'qcflag_visual_bbp660',...
    'qcflag_visual_bbp700',...
    'qcflag_range_temp',...
    'qcflag_range_cond',...
    'qcflag_range_salin',...
    'qcflag_range_sigmath0',...
    'qcflag_range_oxygenAA',...
    ['qcflag_range_' chl_var],...
    'qcflag_range_bbp470',...
    'qcflag_range_bbp650',...
    'qcflag_range_bbp660',...
    'qcflag_range_bbp700',...
    ['qcflag_stuck_' chl_var],...
    'qcflag_stuck_bbp470',...
    'qcflag_stuck_bbp650',...
    'qcflag_stuck_bbp660',...
    'qcflag_stuck_bbp700'};
ncvars = {'temp_qcflag',... % variable qc flag for ncfile
    'cond_qcflag',...
    'salin_qcflag',...
    'sigma_qcflag',...
    'optode_oxygen_qcflag',...
    'chl_cal_qcflag',...
    'bbp470_qcflag',...
    'bbp650_qcflag',...
    'bbp660_qcflag',...
    'bbp700_qcflag',...
    'qcflag_pitch',...
    'qcflag_densinv',...
    'qcflag_biofoul',...
    'qcflag_biofoulOptics',...
    'qcflag_visual_general',...
    'qcflag_visual_CTflag',...
    'qcflag_visual_aanderaa',...
    'qcflag_visual_chl',...
    'qcflag_visual_bbp470',...
    'qcflag_visual_bbp650',...
    'qcflag_visual_bbp660',...
    'qcflag_visual_bbp700',...
    'qcflag_range_temp',...
    'qcflag_range_cond',...
    'qcflag_range_salin',...
    'qcflag_range_sigma',...
    'qcflag_range_oxygen',...
    'qcflag_range_chla',...
    'qcflag_range_bbp470',...
    'qcflag_range_bbp650',...
    'qcflag_range_bbp660',...
    'qcflag_range_bbp700',...
    'qcflag_stuck_chl',...
    'qcflag_stuck_bbp470',...
    'qcflag_stuck_bbp650',...
    'qcflag_stuck_bbp660',...
    'qcflag_stuck_bbp700'};
sgunits = {'QC pass = 1','QC probably bad = 3','QC fail = 4'};
sglongname = {'temperature QC flag',... % longer descriptor
    'conductivity QC flag',...
    'salinity QC flag',...
    'sigma QC flag',...
    'oxygen QC flag',...
    'chla QC flag',...
    'bbp470 QC flag',...
    'bbp650 QC flag',...
    'bbp660 QC flag',...
    'bbp700 QC flag',...
    'pitch QC flag',...
    'density inversion QC flag',...
    'biofoul QC flag',...
    'biofoul optics QC flag',...
    'Visual QC flag general',...
    'Visual QC flag CT sail',...
    'Visual QC flag aanderaa',...
    'Visual QC flag chla',...
    'Visual QC flag bbp470',...
    'Visual QC flag bbp650',...
    'Visual QC flag bbp660',...
    'Visual QC flag bbp700',...
    'Range QC flag temp',...
    'Range QC flag cond',...
    'Range QC flag salin',...
    'Range QC flag sigma',...
    'Range QC flag oxygen',...
    'Range QC flag chla',...
    'Range QC flag bbp470',...
    'Range QC flag bbp650',...
    'Range QC flag bbp660',...
    'Range QC flag bbp700',...
    'Stuck QC flag chla',...
    'Stuck QC flag bbp470',...
    'Stuck QC flag bbp650',...
    'Stuck QC flag bbp660',...
    'Stuck QC flag bbp700'};

% loop over vars
for ii = 1:numel(ncvars)
    if sum(strcmp(MQC.Properties.VariableNames,sgvar{ii})) == 1 % case var exists
        var = MQC.(sgvar{ii});
        if islogical(var)
            var = double(var);
        end
    else
        continue
    end
    varname = ncvars{ii};
    varunits = sgunits{1};
    varlongname = sglongname{ii};
    nccreate(ncqcfile,varname,'Dimensions',{varname,rows});
    ncwrite(ncqcfile,varname,var);
    ncwriteatt(ncqcfile,varname,'long name',varlongname);
    ncwriteatt(ncqcfile,varname,'units',varunits);

end
%save UP and DWN casts
ncdisp(ncqcfile)
