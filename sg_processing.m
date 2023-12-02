function sg_processing(missionname)

addpath(genpath('glider-kit'))

% Notes on toolbox dependencies
% 1. Gibbs Seawater Toolbox needed to run code
% (https://www.teos-10.org/software.htm)

% Level 1 data - includes basic tempc, condc, salin, sigmath0 processing
% combine all raw unprocessed dive files
sg_getRawdata(missionname);

% Level 2 data QC/QC & Calibrations/Corrections
% qc & reprocess hydrography
sg_processHydrography(missionname)

% process Aanderaa oxygen
% Note: tauoptode optimization in beta,
% 30s default, 35s for fast response optodess
if strcmp(missionname,'sg626_1') || strcmp(missionname,'sg512_1')
    tauoptode = 35;
else
    tauoptode = 30;
end
sg_processAAoxygen(missionname,tauoptode)

% process SBE oxygen (SBE43)
% Note: data needs furth quality control
%sg_processSBE43oxygen(missionname)

% process bio-optics (ECO Puck)
sg_convertBiooptics(missionname)

% bbp
sg_processBbp(missionname)

% chlorophyll a
sg_processChla(missionname)
% (NPQ correction in beta)

% CDOM
% Note: data not trustworthy due to fluoresence interference with
% instrument coating
%sg_processCDOM(missionname)

% Level 3 data Gridded on 2m bin
% DWN and UP data save with QC = 1 flags only
% sg_metadata also saved with dived, diveu structures
sg_level3grid(missionname)

% NC file for public use
makeNetCDFfile(missionname)

% Level 4 data - extras calculated from L3 data
% Future versions:
% cphyto
% size fractions
% POC from bbp
% GOP from o2 data, bbp?

end
