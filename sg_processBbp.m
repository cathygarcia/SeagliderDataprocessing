% SG_PROCESSSBBP
% -------------------------------------------------------------------------
% Convert Sea-Bird/WET Lab ECO Puck Volume scattering function (VSF, beta) 
% to particulate bckscattering coefficient, bbp, after beta seawater subtraction, 
% apply quality control flags, and subtract secondary dark bbp minimum
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
% 2. Set bbp wavelength (lambda) and qc range values
% 3. Convert beta to bbp
% 4. Apply quality control flags
%   - visual qc (after all qc, remaining bad data outliers)
%   - stuck qc (argo test, 50% cutoff)
%   - range qc (glogal argo test)
%   - biofoul qc (from oxygen, bbp code)
%   - temp,salin qc input flag (beta to bbp conversion uses these inputs)
% 4. Dark bbp subtraction
% 5. Briggs 2011 spike code applied
% 6. Save tables to qcdata folder 'data/qcdata/missionname_L2.mat'

% MATLAB CODE NEEDED
% -------------------------------------------------------------------------
% betasw_ZHH2009
% sg_greylistQC
% sg_spikes
% 
% -------------------------------------------------------------------------
% Author Catherine Garcia -- cathy.garcia@hawaii.edu -- Version 04 Aug 2023
%
% -------------------------------------------------------------------------


function sg_processBbp(missionname)

filename = ['data/rawdata_L1/' missionname '_L1.mat'] ; %file path to mission process data
load(filename,'DWN0','UP0');
filename = ['data/qcdata_L2/' missionname '_L2.mat'] ; %file path to mission process data
load(filename,'DWN','UP','calib','QC_UP','QC_DWN');

% inputs for for loop
lambdas = [470,700,650,660]; % bbp wavelength
rangeval = [.1 4E-2 4E-2 4E-2]; % for range test
bbpDarkCounts = struct();
signalvars = {'wlbb2fl_bb470_sig','wlbb2fl_bb700_sig','wlbbfl2_bb1_sig','wlbbfl2_bb1_sig'};

% figure for deep correction
close all
fg = figure();
fg.Position = [1224 659 1195 492];
t = tiledlayout(1,3);
tile = 0;
xlimits = [-4E-4 1.5E-3;-2E-4 1E-3;-2E-4 1E-3;-2E-4 1E-3];

for ll = 1:numel(lambdas)

    lambdavar = num2str(lambdas(ll));

    betavar = ['beta' lambdavar];
    if sum(strcmp(DWN.Properties.VariableNames,betavar)) == 1

        % merge table inputs
        directionIdx = [DWN.direction;UP.direction];
        dT.tempc = [DWN.tempc;UP.tempc];
        dT.sp = [DWN.sp;UP.sp];
        dT.beta = [DWN.(betavar);UP.(betavar)];

        % -------------------------------------------------------------------------
        % particle backscattering coefficients
        chi_p = 1.076;
        ldT = height(dT);
        betasw = NaN(ldT,1);

        % compute scatter at 117 deg due to seawater
        for jj = 1:height(dT)
            [betasw(jj),~,~]= betasw_ZHH2009(lambdas(ll),dT.tempc(jj),117,dT.sp(jj)); % unfortunately this function doesn't work with vectors in T or S
        end

        % calculation for particle backscattering coefficients
        dT.(['bbp' lambdavar]) = 2*pi*chi_p*(dT.beta-betasw);

        % move to dwn up
        DWN.(['bbp' lambdavar]) = dT.(['bbp' lambdavar])(directionIdx == "down");
        UP.(['bbp' lambdavar]) = dT.(['bbp' lambdavar])(directionIdx == "up");

        % -------------------------------------------------------------------------
        % QC Analysis on raw data

        % Visual qc
        [qcup,qcdown] = sg_greylistQC(missionname,UP,DWN,['bbp' lambdavar]);
        QC_UP.(['qcflag_visual_bbp' lambdavar]) = qcup;
        QC_DWN.(['qcflag_visual_bbp' lambdavar]) = qcdown;

        % Optics Stuck Value test
        [qcup,qcdown] = sg_opticsQC_stuckValue(DWN0,UP0,signalvars{ll});
        QC_UP.(['qcflag_stuck_bbp' lambdavar]) = qcup;
        QC_DWN.(['qcflag_stuck_bbp' lambdavar]) = qcdown;

        % Biofoul evidence following phase qc
        % not all bioptics appear to be effected
        % (this will loop over iteself with each lambda, but thats okay)
        QC_DWN.qcflag_biofoulOptics = QC_DWN.qcflag_biofoul;
        QC_UP.qcflag_biofoulOptics = QC_UP.qcflag_biofoul;
        % optics flags
        % mission sg146_5, no biofouling on other sensors
        % mission sg146_6, slight increase in spikes on all bbp's
        % mission sg148_8, clear biofouling on other sensors
        % mission sg148_9, clear biofouling on other sensors
        % mission sg511_20, clear biofouling on other sensors, wlbb2fl even earlier
        % mission sg512_3, no biofouling on other sensors
        % mission sg626_3, clear biofouling on other sensors, maybe earlier?
        if any(strcmp(missionname,{'sg146_5','sg512_3'}))
            QC_DWN.qcflag_biofoulOptics(QC_DWN.qcflag_biofoulOptics == 4) = 2;
            QC_UP.qcflag_biofoulOptics(QC_UP.qcflag_biofoulOptics == 4) = 2;
        end
        if any(strcmp(missionname,{'sg146_6','sg148_8','sg148_9','sg511_20','sg626_3'}))
            QC_DWN.qcflag_biofoulOptics(QC_DWN.qcflag_biofoulOptics == 4) = 3;
            QC_UP.qcflag_biofoulOptics(QC_UP.qcflag_biofoulOptics == 4) = 3 ;
        end

        % Flag all variables (and related ones) as fail for out of acceptable
        % range
        % global range test from argo float QC/QC DAC
        % range val set before loop
        QC_UP.(['qcflag_range_bbp' lambdavar]) = ones(height(UP),1);
        QC_DWN.(['qcflag_range_bbp' lambdavar]) = ones(height(DWN),1);
        QC_UP.(['qcflag_range_bbp' lambdavar])(UP.(['bbp' lambdavar]) < -(rangeval(ll)) | UP.(['bbp' lambdavar]) > rangeval(ll)) = 4; % bbp (m-1)
        QC_DWN.(['qcflag_range_bbp' lambdavar])(DWN.(['bbp' lambdavar]) < -(rangeval(ll))  | DWN.(['bbp' lambdavar]) > rangeval(ll)) = 4; % bbp (m-1)

        
        % 4. Assign variable qc flags based on above tests
        % bbp flag
        flags = {'qcflag_visual_general',['qcflag_visual_bbp' lambdavar],'qcflag_biofoulOptics',['qcflag_range_bbp' lambdavar],['qcflag_stuck_bbp' lambdavar],'temp_qcflag','salin_qcflag'};
        iset3 = strcmp(flags,'qcflag_temp') | strcmp(flags,'qcflag_salin');
        var_flags = QC_UP{:,flags};
        var_4mask = var_flags == 4; var_4mask(:,~iset3) = false;
        var_flags(var_4mask) = 2;
        QC_UP.(['bbp' lambdavar '_qcflag']) = max(var_flags,[],2,'omitnan');
        var_flags = QC_DWN{:,flags};
        var_4mask = var_flags == 4; var_4mask(:,~iset3) = false;
        var_flags(var_4mask) = 2;
        QC_DWN.(['bbp' lambdavar '_qcflag']) = max(var_flags,[],2,'omitnan');


        % -------------------------------------------------------------------------
        % Get spikeflag for bbp data
        [sgspikeDWN] = sg_spikes(DWN,{['bbp' lambdavar]});
        [sgspikeUP] = sg_spikes(UP,{['bbp' lambdavar]});
        % bbp at lambda
        DWN.(['bbp' lambdavar '_spikeflag']) = ~isnan(sgspikeDWN.(['bbp' lambdavar '_spikes']));
        UP.(['bbp' lambdavar '_spikeflag']) = ~isnan(sgspikeUP.(['bbp' lambdavar '_spikes']));

        % -------------------------------------------------------------------------
        % Deep dark correction (see Thomalla 2017, Bol 2018)

        % make 2nd percentile correction
        % loop through wavelengths
        spikeflagDWN = DWN.(['bbp' lambdavar '_spikeflag']);
        spikeflagUP = UP.(['bbp' lambdavar '_spikeflag']);
        deepbbpDWN = prctile(DWN.(['bbp' lambdavar])(DWN.vmdepth >= 190 & DWN.vmdepth <= 200 & ~spikeflagDWN & QC_DWN.(['bbp' lambdavar '_qcflag']) < 3),.2); % DWN
        deepbbpUP = prctile(UP.(['bbp' lambdavar])(UP.vmdepth >= 190 & UP.vmdepth <= 200 & ~spikeflagUP & QC_UP.(['bbp' lambdavar '_qcflag']) < 3),.2); % UP
        deepbbp = (deepbbpDWN + deepbbpUP)./2;
        DWN.(['bbp' lambdavar 'c']) = DWN.(['bbp' lambdavar]) - deepbbp;
        UP.(['bbp' lambdavar 'c']) = UP.(['bbp' lambdavar]) - deepbbp;
        bbpDarkCounts.(['bbp' lambdavar]) = deepbbp; % save for output


        % -------------------------------------------------------------------------
        % Get spikeflag for bbp data, but save it this time
        opticVars = {['bbp' lambdavar],['bbp' lambdavar 'c']};
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

        % figure of despiked data
        tile = tile + 1;
        nexttile(tile);
        x = [DWN.([opticVars{1} '_despiked']);UP.([opticVars{1} '_despiked'])];
        xc = [DWN.([opticVars{2} '_despiked']);UP.([opticVars{2} '_despiked'])];
        y = [DWN.vmdepth;UP.vmdepth];
        iqc = [QC_DWN.(['bbp' lambdavar '_qcflag']);QC_UP.(['bbp' lambdavar '_qcflag'])];
        iqc = iqc <3;
        scatter(x(iqc),y(iqc),'ko','filled');
        hold on;
        xline(deepbbp,'r--','LineWidth',1);
        xline(0,'g-');
        scatter(xc(iqc),y(iqc),'r.');
        set(gca,'YDir','reverse');
        title(['bbp ' lambdavar]);
        xlim(xlimits(ll,:));
    end

end

% save figure
if ~isfolder('figures/bbpCorrection'); mkdir('figures/bbpCorrection'); end
title(t,missionname,'Interpreter','none');
exportgraphics(gcf,['figures/bbpCorrection/' missionname '.png'],'Resolution',150);

% save processed data
calib.bbpDarkCounts = bbpDarkCounts;
save(['data/qcdata_L2/' missionname '_L2'],'UP','DWN','calib','QC_DWN','QC_UP');

end