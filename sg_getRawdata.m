function sg_getRawdata(missionname)

addpath(genpath('garcia_glider-kit/')); % make sure path to processing code is valid

workingdir = pwd;
sgdir = [workingdir '/data/unprocessed_L0/' missionname '/'] ; %file path to dives

% get calib file
cd(sgdir);
sg_calib_constants
save([missionname '_calib.mat']);
calib = load([missionname '_calib.mat']);

% get glider number
glider = strsplit(missionname,'_');
glider = glider{1};
glider = strsplit(glider,'sg');
glider = glider{2};

% TOTAL DIVES
divefileInfo = dir([sgdir '*.eng']); % .eng files contain data
nd = size(divefileInfo,1);
diveList = NaN(nd,1);
for dd = 1:nd % get min and max dive from file names
    tempStr = strsplit(divefileInfo(dd).name,{['p' num2str(glider)],'.eng'});
    diveList(dd) = str2double(tempStr{2});
end
diveRange = [min(diveList) max(diveList)];

clear divefileInfo diveList nd dd tempStr

% outputMode
outputMode = 'table';

[UP0,DWN0] = sg_cat(diveRange,calib,'outputMode',outputMode);
UP0.index = (1:height(UP0))';
DWN0.index = (1:height(DWN0))';

% save UP and DWN casts
cd(workingdir)
save(['data/rawdata_L1/' missionname '_L1'],'UP0','DWN0','calib');

end