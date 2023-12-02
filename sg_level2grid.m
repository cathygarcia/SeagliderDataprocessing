% SG_LEVEL2GRID
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

function [sgu,sgd,diveu,dived] = sg_level2grid(DWN,UP,QC_UP,QC_DWN)


% -------------------------------------------------------------------------
% Grid data every 2m *** add iso grid
depth = 2:2:600;
outVarNames = {'vmdepth','chl1','chl2','sigmath0','dateUTC','dateHST','divenum','lat','lon'};
[TF,~] = ismember(outVarNames,DWN.Properties.VariableNames); % only keep valid names
outVarNames = outVarNames(TF);
diveRange = [min(DWN.divenum) max(DWN.divenum)];

[UG, DG] = sg_grid_l2(DWN,UP,depth,'vmdepth',diveRange(1):diveRange(2),outVarNames);


% -------------------------------------------------------------------------
% Put data on a regular grid (dives,depth)
varNamesOut = DG.Properties.VariableNames;
varNamesIn = DG.Properties.VariableNames;
varUnits = {};
dives = unique([UP.divenum(~isnan(UP.divenum));DWN.divenum(~isnan(DWN.divenum))]);

[sgd,sgu,dived,diveu] = putOnGrid_l2(DG,UG,varNamesIn,varNamesOut,varUnits,depth,dives);


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


%% -------------------------------------------------------------------------
% function to place date on 2m grid

    function [UG, DG] = sg_grid_l2(DWN,UP,z_grid,gridVar,diveRange,outVars)

        % check dimensions
        % if ~ismember(gridVar,outVars)
        %     outVars = [gridVar outVars];
        % end

        outVars = unique([gridVar outVars]);

        alldives = [DWN.divenum;UP.divenum];
        if isempty(diveRange)
            dvs = unique(alldives);
        else
            dvs = intersect(diveRange,alldives);
        end
        ndvs = length(dvs);
        nz = length(z_grid);
        szD = size(DWN);
        szU = size(UP);

        % Interpolate and concatenate data
        % a very small random number is added so that there are no duplicate values
        % in interpolation
        DWN.(gridVar) = DWN.(gridVar) + (rand(szD(1),1).*1e-7);
        UP.(gridVar) = UP.(gridVar) + (rand(szU(1),1).*1e-7);

        DG = nan(nz.*ndvs,length(outVars));
        UG = DG;

        for ii = 1:ndvs

            % copy and order data from an individual dive
            D = sortrows(DWN(DWN.divenum == dvs(ii),outVars),gridVar);
            %     disp('sg_grid: Following Dives (if any listed) have NaNs present in gridvar');
            %     if sum(isnan(D.(gridVar))) > 0
            %         disp(num2str(ii))
            %     end
            D(isnan(D.(gridVar)),:) = [];
            U = sortrows(UP(UP.divenum == dvs(ii),outVars),gridVar);
            U(isnan(U.(gridVar)),:) = [];
            % inputs to naninterp1 function
            x = D.(gridVar);
            Y = nan(size(D{:,outVars}));
            for vv = 1:numel(outVars)
                Y(:,vv) = D{:,vv};
            end
            xi = z_grid;
            % gridded array for a single dive
            DG(1+(ii-1)*nz:ii*nz,:) = naninterp1(x, Y, xi);
            clear x Y xi vv
            try
                % inputs to naninterp1 function
                x = U.(gridVar);
                Y = nan(size(U{:,outVars}));
                for vv = 1:numel(outVars)
                    Y(:,vv) = U{:,vv};
                end
                xi = z_grid;
                UG(1+(ii-1)*nz:ii*nz,:) = naninterp1(x, Y, xi);
                clear x Y xi vv
                %     DG(1+(ii-1)*nz:ii*nz,1) = z_grid;
                %     UG(1+(ii-1)*nz:ii*nz,1) = z_grid;
            catch
                debug = 1;
            end
        end

        %convert full array back to table
        DG = array2table(DG,'VariableNames',outVars);
        UG = array2table(UG,'VariableNames',outVars);
        DG.([gridVar '_g']) = repmat(reshape(z_grid,length(z_grid),1),ndvs,1);
        UG.([gridVar '_g']) = repmat(reshape(z_grid,length(z_grid),1),ndvs,1);

    end

%% -------------------------------------------------------------------------
% function to interpolate over NaNs

    function [out] = naninterp1(x, Y, xi)
        MIN_DATA_POINTS = 5;
        szY = size(Y);
        out = nan(length(xi),szY(2));
        for ii = 1:szY(2)
            gd = ~isnan(Y(:,ii));
            % at least two points needed to interpolate
            if sum(gd) >= MIN_DATA_POINTS
                out(:,ii) = interp1(x(gd),Y(gd,ii),xi);
            else
                out(:,ii) = nan.*xi;
            end
        end
    end


% -------------------------------------------------------------------------
% function to reshape data
    function [sgd,sgu,dived,diveu] = putOnGrid_l2(DG,UG,varNamesIn,varNamesOut,varUnits,depth,dives)

        % RESHAPE dimensions
        ld = length(depth);
        nd = numel(dives);

        % Time & Position always included
        % DOWN CAST
        sgd = table(depth',reshape(DG.dateUTC,ld,nd)); sgd.Properties.VariableNames = ({'depth','date'});
        sgd.lon = reshape(DG.lon,ld,nd);
        sgd.lat = reshape(DG.lat,ld,nd);
        % UP CAST
        sgu = table(depth',reshape(UG.dateUTC,ld,nd)); sgu.Properties.VariableNames = ({'depth','date'});
        sgu.lon = reshape(UG.lon,ld,nd);
        sgu.lat = reshape(UG.lat,ld,nd);

        % Water properties
        for ii = 1:numel(varNamesIn)
            sgd.(varNamesOut{ii}) = reshape(DG.(varNamesIn{ii}),ld,nd);
            sgu.(varNamesOut{ii}) = reshape(UG.(varNamesIn{ii}),ld,nd);
        end

        % Variable Units
        if ~isempty(varUnits) == 1
            sgu.Properties.VariableUnits = varUnits;
        end

        % Dive information always included
        % DOWN CAST
        dived.lon = mean(sgd.lon,1,'omitnan');
        dived.lat = mean(sgd.lat,1,'omitnan');
        dived.dive = dives;
        dived.dateUTC = mean(sgd.dateUTC,1,'omitnan');
        dived.hourUTC = dived.dateUTC - fix(dived.dateUTC);
        dived.dateHST = mean(sgd.dateHST,1,'omitnan');
        dived.hourHST = dived.dateHST - fix(dived.dateHST);
        % UP CAST
        diveu.lon = mean(sgu.lon,1,'omitnan');
        diveu.lat = mean(sgu.lat,1,'omitnan');
        diveu.dive = dives;
        diveu.dateUTC = mean(sgu.dateUTC,1,'omitnan');
        diveu.hourUTC = diveu.dateUTC - fix(diveu.dateUTC);
        diveu.dateHST = mean(sgu.dateHST,1,'omitnan');
        diveu.hourHST = diveu.dateHST - fix(diveu.dateHST);

    end


end % function end