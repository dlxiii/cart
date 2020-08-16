%%
%%%------------------------------------------------------------------------
%%%                       DIRECTORY CONFIGURATION
%%%                       目录设置
%%%------------------------------------------------------------------------

% If you are developing the code
% "1" means create file and take time
% "2" means load developed mat files
% "3" means load model object mat files
develop_mode = 1;

% Which system am I using?
if ismac    % On Mac
    basedir = '/Users/yulong/GitHub/';
    addpath([basedir,'fvcomtoolbox/']);
    addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
    addpath([basedir,'fvcomtoolbox/utilities/']);
    addpath([basedir,'fvcomtoolbox/custom/']);
    addpath([basedir,'TMD/']);
    addpath([basedir,'TMD/FUNCTIONS/']);
    addpath([basedir,'TMD/DATA/']);
elseif isunix       % Unix?
    basedir = '/home/usr0/n70110d/github/';
    addpath([basedir,'fvcomtoolbox/']);
    addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
    addpath([basedir,'fvcomtoolbox/utilities/']);
    addpath([basedir,'TMD/']);
    addpath([basedir,'TMD/FUNCTIONS/']);
    addpath([basedir,'TMD/DATA/']);
elseif ispc     % Or Windows?
    basedir = 'C:/Users/Yulong WANG/Documents/GitHub/';      
    addpath([basedir,'fvcomtoolbox/']);
    addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
    addpath([basedir,'fvcomtoolbox/utilities/']);
    addpath([basedir,'fvcomtoolbox/custom/']);
    addpath([basedir,'TMD']);
    addpath([basedir,'TMD/FUNCTIONS/']);
    addpath([basedir,'TMD/DATA/']);
end

% Output directory: dat and nc files.
inputConf.outbase = pwd;

% Write out diary file.
% Clear the diary file if it does exist.
if exist(fullfile(inputConf.outbase, 'diary'),'file')==2
    try
        diary off;
        delete('diary');
    catch
        delete('diary');
    end
end
diary on;

% Write out all the required files.
% Make the output directory if it doesnt exist
if exist(inputConf.outbase, 'dir')~=7
    mkdir(inputConf.outbase)
end

% Working directory: grads folder and current folder
inputConf.base = [basedir,'cart/examples_canal/mesh/canal/'];

% Which version of FVCOM are we using (for the forcing file formats)?
inputConf.FVCOM_version = '4.0';

% Location of grads file
inputConf.grid = [inputConf.base,...
    'Mesh/Area Property Mesh/Area Property Mesh.2dm'];

% Case name for the model inputs and outputs
% Change this to whatever you want
inputConf.casename = 'canal';

% Case name for the model inputs and outputs
% Change this to whatever you want
inputConf.flag.sph = 'no';
inputConf.flag.pro = 'no';
inputConf.flag.tid = 'no';
inputConf.flag.riv = 'yes';
inputConf.flag.dye = 'yes';
inputConf.flag.sta = 'no';

%%
%%%------------------------------------------------------------------------
%%%                           Spatial stuff
%%%                           空间设置
%%%------------------------------------------------------------------------

% input coordinates (what's my input bathy in?)
inputConf.coordInput = 'spherical'; % 'spherical' or 'cartesian'
% output coordinates (FVCOM only likes cartesian at the moment)
if strcmpi(inputConf.flag.sph, 'no')
    inputConf.coordType = 'cartesian'; % 'spherical' or 'cartesian' 
else
    inputConf.coordType = 'spherical'; % 'spherical' or 'cartesian' 
end

% Input grid UTM Zone (if applicable)
% See: https://upload.wikimedia.org/wikipedia/commons/e/ed/Utm-zones.jpg
% As Utm-zones indicated, here should be 54-S, but here the utmZone should
% be tmzone (UTM longitudinal zone) and utmhemi (UTM hemisphere as array of
% 'N' or 'S' characters)
inputConf.utmZone = {'54 N'};

% Option to smooth the bathymetry data.
inputConf.smoothBathy = 'no'; % 'yes' or 'no'.
if strcmpi(inputConf.smoothBathy, 'yes')
    % Set the smoothing factor and number of iterations (see smoothmesh).
    inputConf.smoothFactors = [0.5, 4]; % [factor, iterations]
end

% vertical coordinates type: sigma or hybrid
inputConf.verticalCoordType = 'sigma';

%%
%%%------------------------------------------------------------------------
%%%                     Time
%%%                     时间设置
%%%------------------------------------------------------------------------

% Model time ([Y,M,D,h,m,s])
inputConf.modelYear = 2014;
inputConf.startDate = [inputConf.modelYear,01,01,00,00,00];
inputConf.endDate = [inputConf.modelYear+1,01,01,00,00,00];

% Convert times to Modified Julian Date
inputConf.startDateMJD = greg2mjulian(inputConf.startDate(1),...
    inputConf.startDate(2),inputConf.startDate(3),inputConf.startDate(4),...
    inputConf.startDate(5),inputConf.startDate(6));
inputConf.endDateMJD = greg2mjulian(inputConf.endDate(1),...
    inputConf.endDate(2),inputConf.endDate(3),inputConf.endDate(4),...
    inputConf.endDate(5),inputConf.endDate(6));

% The number of months in the period of data.
inputConf.mm = inputConf.startDate(2):inputConf.endDate(2);
if inputConf.mm == 1
    inputConf.dOffsets = [0, 4];
elseif inputConf.mm == 12
    inputConf.dOffsets = [2, 0];
else
    inputConf.dOffsets = [2, 4];
end

%%
%%%------------------------------------------------------------------------
%%%                     Model constants
%%%                     模型常数设置
%%%------------------------------------------------------------------------

% Sponge layer parameters
inputConf.spongeRadius = 1; % in metres, or -1 for variable
inputConf.spongeCoeff = 0.001;

% z0 value in metres
inputConf.bedRoughness = 0.015; % or 0.015, 0.025 or 0.03 - Davies and Furnes (1980) shelf model

% Estimated velocity (m/s) and tidal range (m) for time step estimate
inputConf.estVel = 1.0;
inputConf.estRange = 2.0;

%%
%%%------------------------------------------------------------------------
%%%                       Forcing and stuff
%%%                       驱动力基本设置
%%%------------------------------------------------------------------------

% Model time type ('non-julian' and 'julian')
inputConf.datetype = 'julian';
% Increment used for tide (days)
inputConf.datetide = 1/24;
% Model time is decided depend on the datetype.
if strcmpi(inputConf.datetype,'non-julian') 
    inputConf.modelTime = [...
        inputConf.startDateMJD, ...
        inputConf.endDateMJD];
elseif strcmpi(inputConf.datetype,'julian')
    inputConf.modelTime = [...
        inputConf.startDateMJD - inputConf.dOffsets(1), ...
        inputConf.endDateMJD + inputConf.dOffsets(2)];
end
% Open boundary forcing nodal forcing type (drived by OTPS).
inputConf.obcForcing = 'z'; 
inputConf.tidesMJD = inputConf.startDateMJD:inputConf.datetide:inputConf.endDateMJD;

% % Increment used for open boundary ST (days)
% inputConf.dateobs = 1/24;
% % Open boundary temperatures and salinities (string for source or number for constant).
% % The data is avaliable from 1992-10-02 00:00:00
% inputConf.obc_temp = 'HYCOM';
% inputConf.obc_salt = 'HYCOM';
% inputConf.obctsMJD = [inputConf.startDateMJD, inputConf.endDateMJD + 1];

% % Increment used for surface forcing (days)
% inputConf.dateforcing = 1/24;
% % The surface forcing from NCEP
% % The data is avaliable from 1948 - present
% inputConf.doForcing = 'NCEP-CALCULATED';
% inputConf.doForcing = 'GWO';
% if strcmpi(inputConf.doForcing, 'GWO')
%     inputConf.forceMJD = inputConf.startDateMJD:inputConf.dateforcing:inputConf.endDateMJD;
% elseif strcmpi(inputConf.doForcing, 'NCEP')
%     inputConf.forceMJD = [inputConf.startDateMJD, inputConf.endDateMJD];
% elseif strcmpi(inputConf.doForcing, 'NCEP-CALCULATED')
%     inputConf.forceMJD = [inputConf.startDateMJD, inputConf.endDateMJD];
% end

% The river forcing from flux, salinity and temperature
inputConf.riverForcing = 'FLUX';
% River information
inputConf.river.infos = {...
    'test01',...
    'test02',...
    'test03',...
    'test04',...
    };
% Location of river file
inputConf.river.flux = [basedir,'cart/examples_canal/inp/river_flux.csv'];
inputConf.river.temp = [basedir,'cart/examples_canal/inp/river_temp.csv'];
inputConf.river.salt = [basedir,'cart/examples_canal/inp/river_salt.csv'];
inputConf.river.location = [basedir,'cart/examples_canal/inp/river_location.csv'];

% % New information of river, combine Nakagawa and Arakawa data to one river.
% % River information
% inputConf.river.infos = {...
%     'Edogawa',...
%     'Nakagawa_Arakawa',...
%     'Sumidagawa',...
%     'Tamagawa',...
%     'Tsurumigawa',...
%     'Ebigawa'};
% % Location of river file
% inputConf.river.flux = [basedir,'river/data_q/river_flux_combine.csv'];
% inputConf.river.temp = [basedir,'river/data_q/river_temp_combine.csv'];
% inputConf.river.salt = [basedir,'river/data_q/river_salt_combine.csv'];
% inputConf.river.location = [basedir,'river/data_q/river_location_combine.csv'];

% Adjust river mouth location
% 139°55'57.84"	139°50'55.07"	139°46'29.73"	139°46'46.94"	139°40'53.63"	139°58'41.66"
%  35°41'56.29"	 35°38'36.31"	 35°38'49.41"	 35°31'45.11"	 35°28'25.39"	 35°40'52.25"

% Give some names to the boundaries. This must match the number of node
% strings defined in SMS. Ideally, the order of the names should match the
% order in which the boundaries were made in SMS.
inputConf.boundaryNames = {'open boundary'};

% Stations
% Read from kml or list file.
% inputConf.station = 'list';
% inputConf.station = 'kml';
% if strcmpi(inputConf.station,'kml')
%     inputConf.names = {};
%     inputConf.positions = [];
%     inputConf.fkml = fullfile(inputConf.outbase, ...
%         [inputConf.casename,'_stations.kml']);
%     kml = read_kml(inputConf.fkml);
%     for k = 1:size(kml,2)
%         if strcmpi(kml(k).Geometry,'Point')
%             inputConf.names{end+1} = kml(k).Name;
%             inputConf.positions(end+1,1) = kml(k).Lon;
%             inputConf.positions(end,2) = kml(k).Lat;
%         end
%     end
% elseif strcmpi(inputConf.station,'list')
%     inputConf.names = {...
%         'Tokyo',...
%         'Chiba',...
%         'Yokohamashinko',...
%         'Daini-Kaiho',...
%         'Yokosuka',...
%         'Kyurihamako',...
%         'Tidal-Station-1',...
%         'Tidal-Station-2',...
%         'Tidal-Station-3',...
%         'Tidal-Station-4',...
%         };
%     inputConf.positions = [...
%         139.7700000000000,35.648888888888889;...
%         140.0455555555556,35.568055555555556;...
%         139.6441666666666,35.454166666666667;...
%         139.7433333333333,35.308611111111105;...
%         139.6513888888889,35.288055555555556;...
%         139.7208333333333,35.227777777777778;...
%         139.6912777777778,35.062205555555556;...
%         139.7207861111112,35.076750000000004;...
%         139.7504166666667,35.093661111111111;...
%         139.7735055555556,35.127727777777778;...
%         ];
% else
%     inputConf.station = 'None';
% end
% 
% if ~strcmpi(inputConf.station,'None')
% UTMzone = regexpi(inputConf.utmZone,'\ ','split');
% for s = 1:length(inputConf.positions)
%     [inputConf.positions(s,3),inputConf.positions(s,4),~,~] = wgs2utm(...
%         inputConf.positions(s,2),inputConf.positions(s,1),...
%         str2double(char(UTMzone{1}(1))),char(UTMzone{1}(2)));
% end
% clear s UTMzone
% end

% The accepted distance error of cal and real, in the cartesian coordinate 
% system, the dist unit is meter.
inputConf.dist = 1000.00;

%%
%%%------------------------------------------------------------------------
%%%                      Model mesh generation
%%%                      网格生成
%%%------------------------------------------------------------------------

% Read the input mesh and bathymetry. Also creates the data necessary for
% the Coriolis correction in FVCOM.
mesh_str = strsplit(inputConf.grid,'.');
if mesh_str{1,2}=="2dm"
    Mobj = read_2dm_mesh(...
        '2dm',inputConf.grid,...
        'coordinate',inputConf.coordType,...
        'in_coord',inputConf.coordInput,...
        'project',true,...
        'zone',inputConf.utmZone,...
        'addCoriolis',true);
else
    Mobj = read_grid_mesh(...
        'grid',inputConf.grid,...
        'coordinate',inputConf.coordType,...
        'in_coord',inputConf.coordInput,...
        'project',true,...
        'zone',inputConf.utmZone,...
        'addCoriolis',true);
end
clear mesh_str

% Smooth the bathymetry if desired.
if strcmpi(inputConf.smoothBathy, 'yes')
    Mobj = setup_metrics(Mobj);
    % Backup Mobj.h as Mobj.h_backup
    Mobj.h_backup = Mobj.h;
    Mobj.h = smoothfield(Mobj.h, Mobj, ...
        inputConf.smoothFactors(1), inputConf.smoothFactors(2));
    % smoothfield2 is really inappropriate for bathymetry data.
    % Mobj.h = smoothfield2(Mobj.h,Mobj,inputconf.smoothFactors(2));
end

% Parse the open boundary nodes and add accordingly
% Add the sponge nodes
for i=1:size(Mobj.read_obc_nodes,2)
    nodeList = double(cell2mat(Mobj.read_obc_nodes(i)));
    Mobj = add_obc_nodes_list(Mobj,nodeList,inputConf.boundaryNames{i},1,1);
    if inputConf.spongeRadius < 0    % if we want a variable sponge radius
        if i==1
            % Create an array to store the radii
            Mobj.sponge_rad = zeros(size(Mobj.sponge_nodes));
        end
        % calculate the sponge radius
        spongeRadius = calc_sponge_radius(Mobj,nodeList);
        % Add the sponge nodes to the list
        Mobj = add_sponge_nodes_list(Mobj,nodeList,...
            [inputConf.boundaryNames{i},' sponge'],spongeRadius,...
            inputConf.spongeCoeff);
    else
        Mobj = add_sponge_nodes_list(Mobj,nodeList,...
            [inputConf.boundaryNames{i},' sponge'],inputConf.spongeRadius,...
            inputConf.spongeCoeff);
    end
    clear nodeList spongeRadius
end
clear i

% Model time is decided depend on the datetype.
if strcmpi(inputConf.verticalCoordType,'sigma')
    % Get the sigma depths in order to interpolate from the depths
    if exist(fullfile(inputConf.outbase, [inputConf.casename,'_sigma.dat']),'file')
        % If the sigma.dat file exists, read it
        Mobj = read_sigma(Mobj, fullfile(inputConf.outbase, [inputConf.casename,'_sigma.dat']));
    else
        % If we can't find the sigma.dat file, print an error message and finish
        error(['sigma.dat not found. Please put your sigma.dat file into ',...
            fullfile(inputConf.outbase),' and try again.'])
    end
elseif strcmpi(inputConf.verticalCoordType,'hybrid')
    inputConf.hybrid.sigma_file = 'coord_hybrid.dat';
    inputConf.hybrid.nlev = 11;                % number of vertical levels (layers + 1)
    inputConf.hybrid.H0 = 100;                 % transition depth of the hybrid coordinates
    inputConf.hybrid.KU = 2;                   % number of layers in the DU water column
    inputConf.hybrid.KL = 1;                   % number of layers in the DL water column
    inputConf.hybrid.DU = 20;                  % upper water boundary thickness (metres)
    inputConf.hybrid.DL = 10;                  % lower water boundary thickness (metres)
    Mobj = hybrid_coordinate(inputConf.hybrid, Mobj);
end

% Do the bed roughness
Mobj.z0 = ones(1,Mobj.nElems)*inputConf.bedRoughness;

% Generate center point coordination of elements
% Do the stations list
Mobj.xc = nodes2elems(Mobj.x, Mobj);
Mobj.yc = nodes2elems(Mobj.y, Mobj);
Mobj.lonc = nodes2elems(Mobj.lon, Mobj);
Mobj.latc = nodes2elems(Mobj.lat, Mobj);

% Add station
% if ~strcmpi(inputConf.station,'None')
%     Mobj = add_stations_list(Mobj,inputConf.positions,inputConf.names,inputConf.dist);
% end

% Estimate model time step. Supply estimated velocity (m/s) and tidal range
% (m) after the mesh object.
Mobj = estimate_ts(Mobj,inputConf.estVel,inputConf.estRange);
fprintf('Estimated time step:\t%.2f\n',min(Mobj.ts));

%%
%%%------------------------------------------------------------------------
%%%                     Output of basic configurations
%%%                     基本输出
%%%------------------------------------------------------------------------

% Grid
write_FVCOM_grid(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_grd.dat']));
% Bathymetry
write_FVCOM_bath(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_dep.dat']));
% Coriolis
write_FVCOM_cor(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_cor.dat']));
% Open boundaries
write_FVCOM_obc(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_obc.dat']))
% Sponge file
write_FVCOM_sponge(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_spg.dat']))
% Bed roughness (constant or variable (see above))
write_FVCOM_z0(Mobj.z0,fullfile(inputConf.outbase,[inputConf.casename,'_z0.nc']),'bottom roughness');
% % Time series wave stations
% if ~strcmpi(inputConf.station,'None')
%     write_FVCOM_stations(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_station.dat']));
% end
    
% Save Model object file
if exist('varb', 'dir')~=7
    mkdir('varb')
end
save('varb/Mobj_00.mat','Mobj','-v7.3','-nocompression');

%%
%
%%%------------------------------------------------------------------------
%%%                    Additional forcing: Tides
%%%                    潮汐驱动力 
%%%------------------------------------------------------------------------
tic
% Open boundary nodal forcing type
%   'z' for predicted surface elevation
%   'phase-amp' for amplitudes and phases
%   'model-output' for Tidal Model Driver output
fprintf('Calculating open boundary forcing data from OTPS...\n')
% Need to cd to TPXO directory or it doesn't work
inputConf.extractType = 'z'; 
% (Yes, this is inelegant but it's the easiest way for now)
here = pwd; % store the current working directory to return later
tpxo_dir = which('TMD');    % find the TPXO directory
tpxo_dir = tpxo_dir(1:end-5);   % remove TPXO.m
if strcmpi(inputConf.obcForcing, 'z')
    cd(tpxo_dir)    % go to TPXO directory
    % Location of the TMD model description file
    % OhS is not suitable for Kyushu area.
    % Use "China Seas & Indonesia 2016" instead.
    % inputConf.Model = [tpxo_dir,'DATA/Model_OhS'];
    inputConf.Model = [tpxo_dir,'DATA/Model_Ind_2016'];
elseif strcmpi(inputConf.obcForcing, 'phase-amp')
    cd(tpxo_dir)    % go to TPXO directory
    % Location of the TMD model description file
    % OhS is not suitable for Kyushu area.
    % Use "China Seas & Indonesia 2016" instead.
    % inputConf.Model = [tpxo_dir,'DATA/Model_OhS'];
    inputConf.Model = [tpxo_dir,'DATA/Model_Ind_2016'];
elseif strcmpi(inputConf.obcForcing, 'otps')
    fid=fopen([tpxo_dir,'LAT_LON/lat_lon_1st'],'w');
    fprintf('Rewriting %s...\n','lat_lon_1st');
    fprintf(fid,"   lat       lon             yy   mm   dd   hh   mi  sec  dt(min) TSLength\n");
    for i = 1:Mobj.nObcNodes
        fprintf(fid," %10.4f %10.4f %8d %4d %4d %4d %4d %4d %4d %10d\n",...
            Mobj.lat(i), Mobj.lon(i),...
            inputConf.startDate(1), inputConf.startDate(2),...
            inputConf.startDate(3), inputConf.startDate(4),...
            inputConf.startDate(5), inputConf.startDate(6),...
            inputConf.datetide*24*60, 1/inputConf.datetide*(inputConf.endDateMJD-inputConf.startDateMJD)+1);
    end
    fclose(fid);
    fprintf('Finishing rewriting %s...\n','lat_lon_1st');
end

% How many tidal constituents do we actually want to use at the model
% boundaries? Case sensitive (M2 != m2).
% Only relevant if using TPXO.
inputConf.tidalComponents = {'M2','S2','N2','K2','K1','O1','P1','Q1'};
% clear tpxo_dir

%%%------------------------------------------------------------------------
%%%                  Tides and its (non)julian output
%%%------------------------------------------------------------------------

% Generate a surface elevation time series for open boundary forcing.
if strcmpi(inputConf.obcForcing, 'z')
    if develop_mode == 3
        cd(here);
    	fprintf('Loading Model objet file...\n');
        load('varb/Mobj_01.mat');
        fprintf('Done!\n');
    else
        % Use tmd_tide_pred to predict surface elevations for a given time range.
        % Add the tidal components to the Mobj.
        % Mobj.Components = conf.obc_tides.components;
        Mobj.Components = inputConf.tidalComponents;
        % Get the indices to use the tidal constituents defined in
        % conf.obc_tides.components for TPXO (which requires a
        % numerical array of the constituents to be used). The order of the
        % TPXO constituents is M2, S2, N2, K2, K1, O1, P1, Q1, MF, MM, M4,
        % MS4, MN4.
        tpxoConsts = {'M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1'};
        tIndUse = nan(length(Mobj.Components), 1);
        tInd = 1:length(tpxoConsts);
        for i=1:length(Mobj.Components)
            tPos = tInd(strcmp(Mobj.Components{i}, tpxoConsts));
            if ~isempty(tPos)
                tIndUse(i) = tPos;
            else
                warning('Supplied constituent (%s) is not present in the TPXO data', Mobj.Components{i}); %#ok<WNTAG>
            end
        end
        % Tidy up a bit
        clear i c tpxoConsts tPos tInd
        tIndUse = tIndUse(~isnan(tIndUse));
        % We can't just use tmd_tide_pred to do all the surface elevations
        % at once. Instead, the useful approaches are:
        %   1. Time series at a single location
        %   2. Map of a given time step at all locations
        surfaceElevation = nan(Mobj.nObcNodes, size(inputConf.tidesMJD, 2), length(inputConf.boundaryNames));
        for i=1:length(inputConf.boundaryNames)
            for j=1:Mobj.nObcNodes
                % Get the current location (from the node ID)
                currLon = Mobj.lon(Mobj.obc_nodes(i,j));
                currLat = Mobj.lat(Mobj.obc_nodes(i,j));
                %if ftbverbose
                    fprintf('Position %i of %i (%.3f %.3f)... \n', j, Mobj.nObcNodes, currLon, currLat);
                %end
                [surfaceElevation(j,:,i), ~] = tmd_tide_pred(inputConf.Model, ...
                    inputConf.tidesMJD+678942.000000, currLat, currLon, 'z', tIndUse);
                if isnan(surfaceElevation(j,:))
                    % Try the global model instead.
                    [surfaceElevation(j,:,i), ~] = tmd_tide_pred(inputConf.Model, ...
                    inputConf.tidesMJD, currLat, currLon, 'z', tIndUse);
                end
            end
        end
        Mobj.surfaceElevation = surfaceElevation;
        % Tidy up some more
        clear i j tIndUse obc_lat obc_lon currLon currLat surfaceElevation
        cd(here);
        save('varb/Mobj_01.mat','Mobj','-v7.3','-nocompression');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TEST CODE
        % t.Model = '/home/usr0/n70110d/github/fvcomtoolbox/tmd/DATA/Model_OhS'
        % t.SDtime = [datenum([2019,05,14,15,00,00]):1/24/60*5:datenum([2019,05,14,15,00,00])+1]
        % t.lat = 35.5681
        % t.lon = 140.0456
        % t.ptype = 'z'
        % t.Cid = 
        % [t.TS,t.conList]=tmd_tide_pred(t.Model,t.SDtime,t.lat,t.lon,t.ptype)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
elseif strcmpi(inputConf.obcForcing, 'otps')
    if develop_mode == 3
    	fprintf('Loading Model objet file...\n');
        load('varb/Mobj_01.mat');
        fprintf('Done!\n');
    else
        surfaceElevation = nan(Mobj.nObcNodes, size(inputConf.tidesMJD, 2), length(inputConf.boundaryNames));
        for i=1:length(inputConf.boundaryNames)
            for j=1:Mobj.nObcNodes
                % Get variables from the generated mat format files
                filename = [tpxo_dir,'OUT/1st_',mat2str(j),'.mat'];
                fprintf('Getting %s from TMD\n',['1st_',mat2str(j),'.mat']);
                variable = {'SerialDay','TimeSeries'};
                load(filename,variable{:});
                [~,n] = find(SerialDay - 678942 == inputConf.tidesMJD(1));
                for k=1:length(inputConf.tidesMJD)
                    surfaceElevation(j,k,i) = TimeSeries(1,n+k-1);
                end
            end
        end
        Mobj.surfaceElevation = surfaceElevation;
        % Tidy up some more
        clear i j k m n filename variable TimeSeries SerialDay surfaceElevation
        save('varb/Mobj_01.mat','Mobj','-v7.3','-nocompression');
    end
elseif strcmpi(inputConf.obcForcing,'phase-amp')
    if develop_mode == 3
        cd(here);
    	fprintf('Loading Model objet file...\n');
        load('varb/Mobj_01.mat');
        fprintf('Done!\n');
    else
        % Boundary conditions from TPXO (for spectral tides or predicted
        % surface elevations)
        % Put the input list into the mesh object.
        Mobj.Components = inputConf.tidalComponents;
        % Set up the tidal struct. This contains the relevant information for up to
        % eight constituents, ordered as period, beta love number and equilibrium
        % amplitude.
        %                    period    beta    eq. amp.
        %                      (s)      (?)       (m)
        tideComponents.M2 = [44714.16, 0.693, 0.242334];
        tideComponents.S2 = [43200.00, 0.693, 0.112841];
        tideComponents.N2 = [45570.24, 0.693, 0.046398];
        tideComponents.K2 = [43082.28, 0.693, 0.030704];
        tideComponents.K1 = [86163.84, 0.736, 0.141565];
        tideComponents.O1 = [92949.84, 0.695, 0.100514];
        tideComponents.P1 = [86637.24, 0.706, 0.046843];
        tideComponents.Q1 = [96726.24, 0.695, 0.019256];
        %tideComponents.Mf = [1180260,  ?????, 0.041742];
        %tideComponents.Mm = [2380716,  ?????, 0.022026];
        %tideComponents.Ssa = [15778980, ????, 0.019446];
        % Extract the values for each tidal component into Mobj.period_obc,
        % Mobj.beta_love and Mobj.equilibrium_amp.
        for c=1:size(Mobj.Components,2)
            Mobj.period_obc(c) = tideComponents.(Mobj.Components{c})(1);
            Mobj.beta_love(c) = tideComponents.(Mobj.Components{c})(2);
            Mobj.equilibrium_amp(c) = tideComponents.(Mobj.Components{c})(3);
        end
        clear c
        % Provide amplitude and phase data for the boundary nodes. Use the TMD
        % function tmd_extract_HC.m to get harmonic constants at the boundary
        % nodes.
        amp=cell(1,Mobj.nObs);
        Gph=cell(1,Mobj.nObs);
        Depth=cell(1,Mobj.nObs);
        constList=cell(1,Mobj.nObs);
        for i=1:length(inputConf.boundaryNames)
            % It is possible to specify the indices of the constituents of interest
            % when calling tmd_extract_HC, but it requires knowing the order
            % they're stored in the file. Easier for me to extract the constituents
            % of interest separately. This makes it a bit slower (having to
            % interpolate all the constituents is slower than a select few), but
            % it's a bit easier to code up.
            if Mobj.have_lonlat
            [amp{i},Gph{i},Depth{i},constList{i}] = tmd_extract_HC(inputConf.Model,Mobj.lat(Mobj.read_obc_nodes{i}),Mobj.lon(Mobj.read_obc_nodes{i}),inputConf.extractType);
            else
                % Need to convert XY to latlon.
                try % to use the handy file exchange utm2deg function
                    % Make cell array of all the zones because utm2deg is a bit
                    % inflexible in that regard (size of utmZones must equal size
                    % of x and y).
                    % This is somewhat redundant now that the lat/long is added
                    % when generating the Coriolis values, but it's still
                    % worthwhile keeping it here just in case. No harm should
                    % come of it being here anyway.
                    utmZones=cellfun(@(x) repmat(x,length(Mobj.x(Mobj.read_obc_nodes{i})),1),inputConf.utmZone,'uni',false);
                    [tmpLat,tmpLon] = utm2deg(Mobj.x(Mobj.read_obc_nodes{i}),Mobj.y(Mobj.read_obc_nodes{i}),utmZones{1});
                    % Get the tidal data
                    [amp{i},Gph{i},Depth{i},constList{i}] = tmd_extract_HC(inputConf.Model,tmpLat,tmpLon,inputConf.extractType);
                catch %#ok<CTCH>
                    error('Can''t convert X/Y positions to lat/long, so can''t extract data from the TPXO data. Consider adding utm2deg to your PATH.')
                end
            end
            for j=1:numel(Mobj.Components)
                fprintf('Extracting %s... ',Mobj.Components{j});
                posIdx = strmatch(lower(Mobj.Components{j}),constList{i}); %#ok<MATCH2>
                Mobj.amp_obc{i}(j,:) = amp{i}(posIdx,:);
                Mobj.phase_obc{i}(j,:) = Gph{i}(posIdx,:); % Greenwich phase
                fprintf('done!\n');
            end
        end
        clear posIdx amp Gph Depth constList
        % Find NaNs in the boundaries
        for i=1:Mobj.nObs
            brokenBoundary=i;
            nanIdx = Mobj.read_obc_nodes{brokenBoundary}(isnan(Mobj.phase_obc{brokenBoundary}(1,:)));
            nanLon = Mobj.lon(nanIdx);
            nanLat = Mobj.lat(nanIdx);
            inputConf.doFig=0;
            if max(nanLon)-min(nanLon)==0
                minPos = min(nanLat);
                maxPos = max(nanLat);
                inputConf.doFig=1;
            elseif max(nanLat)-min(nanLat)==0
                minPos = min(nanLon);
                maxPos = max(nanLon);
                inputConf.doFig=1;
            elseif isempty(nanIdx)
                fprintf('No NaNs in %s boundary.\n',inputConf.boundaryNames{i});
                clear nanLon nanLat nanIdx
            else
                error('Boundaries are not linear. Won''t plot %s boundary NaNs',inputConf.boundaryNames{i});
            end
            if inputConf.doFig
                figure
                patch('Vertices',[Mobj.lon,Mobj.lat],'Faces',Mobj.tri,...
                    'Cdata',Mobj.h,'edgecolor','k','facecolor','interp');
                hold on;
                plot(Mobj.lon(nanIdx),Mobj.lat(nanIdx),'wo','LineWidth',3,'MarkerSize',12);
                plot(Mobj.lon(nanIdx),Mobj.lat(nanIdx),'ko','LineWidth',3,'MarkerSize',8);
                axis('equal','tight');
            end
        end
        save('varb/Mobj_01.mat','Mobj','-v7.3','-nocompression');
    end
end

if strcmpi(inputConf.obcForcing, 'z')
    % Write out the TPXO predicted surface elevation.
    write_FVCOM_elevtide(Mobj, ...
        inputConf.tidesMJD,...
        fullfile(inputConf.outbase, [inputConf.casename, '_julian_obc.nc']),...
        'Model surface elevation boundary input',...
        'floattime', true,...
        'julian', true);
elseif strcmpi(inputConf.obcForcing, 'otps')
    % Write out the TPXO predicted surface elevation.
    write_FVCOM_elevtide(Mobj, ...
        inputConf.tidesMJD,...
        fullfile(inputConf.outbase, [inputConf.casename, '_julian_obc.nc']),...
        'Model surface elevation boundary input',...
        'floattime', true,...
        'julian', true);
elseif strcmpi(inputConf.obcForcing,'phase-amp')  
    % Write out the TPXO predicted spectral tide.
    set_spectide(Mobj,...
        numel(Mobj.Components),...
        fullfile(inputConf.outbase,[inputConf.casename,'_non_julian_obc.nc']),...
        'TPXO spectral tidal boundary input');
end
fprintf('Tidal forcing working time: %.2f minutes\n', toc / 60);
%

%%
%%%------------------------------------------------------------------------
%%%                    HYCOM S&T forcing and staff
%%%------------------------------------------------------------------------
% tic
% % Now we need some boundary temperature and salinity conditions.
% if strcmpi('HYCOM', {inputConf.obc_temp, inputConf.obc_salt})
%     % Use HYCOM data for the boundary forcing.
%     % Offset the times to give us a bit of wiggle room.
%     if develop_mode == 3
%     	fprintf('Loading Model objet file...\n')
%         load('varb/Mobj_02.mat');
%         fprintf('Done!\n');
%     else
%     	if develop_mode == 1
%     		fprintf('Downloading daliy open boundary S&T forcing from HYCOM...\n');
%         	hycom = get_HYCOM_forcing(Mobj, inputConf.obctsMJD, {'temperature', 'salinity'}); 
%         	save(['hycom','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat'],'hycom','-v7.3','-nocompression');
%         	fprintf('Downloading daliy open boundary S&T forcing from HYCOM...Done!\n');
%         elseif develop_mode == 2
%     		fprintf('Loading daliy open boundary S&T forcing from local HYCOM database...\n')
%         	load(['hycom','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat']);
%         	fprintf('Downloading daliy open boundary S&T forcing from HYCOM...Done!\n');
%     	end
%     	% Interpolate the 4D HYCOM data on the FVCOM vertical grid at the open boundaries.
%     	Mobj = get_HYCOM_tsobc(Mobj, hycom, {'temperature'});
%     	Mobj = get_HYCOM_tsobc(Mobj, hycom, {'salinity'});
%         % Interpolate the 4D HYCOM data on the hourly time series
%         Mobj = get_HYCOM_series(Mobj, inputConf.dateobs);
%         save('varb/Mobj_02.mat','Mobj','-v7.3','-nocompression');
%     end
% elseif strcmpi('FRA-JCOPE', {inputConf.obc_temp, inputConf.obc_salt})
%     % [data,header]=read_grads(file_name,var_name,varargin)
%     % 
%     % file_name = ['C:\Users\Yulong WANG\Documents\GitHub\jcope-convert\fra_jcope\el.ctl'];
%     % var_name = ['all'];
%     % [data,header]=read_grads('C:\Users\Yulong WANG\Documents\GitHub\jcope-convert\fra_jcope\t.ctl','all'); 
% end
% fprintf('Open boundary ST forcing making time: %.2f minutes\n', toc / 60)
% 
% tic
% % Write the temperature and salinity.
% if strcmpi('HYCOM', {inputConf.obc_temp, inputConf.obc_salt})
%     fprintf('Writing daliy open boundary S&T forcing file.\n')
%     write_FVCOM_tsobc(fullfile(inputConf.outbase, inputConf.casename), ...
%         Mobj.ts_times, ...
%         size(Mobj.temperature, 2), ...
%         Mobj.temperature, ...
%         Mobj.salt, ...
%         Mobj, ...
%         'floattime', true,...
%         'julian', true)
% end
% fprintf('Open boundary ST forcing writing time: %.2f minutes\n', toc / 60)
% 
% % plot temp and salt
% % plot(datetime(Mobj.ts_times+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),squeeze(Mobj.temperature(1,1,:)));
% % plot(datetime(Mobj.ts_times+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),squeeze(Mobj.salt(1,1,:)));
% 
% %%%------------------------------------------------------------------------
% %%%                     Meteorology and output
% %%%------------------------------------------------------------------------
% % Get the surface heating data.
% if strcmpi(inputConf.doForcing, 'NCEP')
%     % Use the OPeNDAP NCEP script to get the following parameters:
%     %     - Potential evaporation rate (pevpr)[W/m^2]                   : for extracting land mask
%     %     - u wind component (uwnd)[m/s]                                : for wind
%     %     - v wind component (vwnd)[m/s]                                : for wind
%     %     - Precipitation rate (prate)[Kg/m^2/s]                        : for precipitation
%     %     - Latent Heat Net Flux at Surface (lhtfl)[W/m^2]              : for evaporation
%     %     - Sensible Heat Net Flux at Surface (shtfl)[W/m^2]
%     %     - Downward Longwave Radiation Flux at Surface (dlwrf) [W/m^2] : for heat flux
%     %     - Downward Solar Radiation Flux at Surface (dswrf) [W/m^2]    : for heat flux
%     %     - Upward Longwave Radiation Flux at Surface (ulwrf) [W/m^2]   : for heat flux
%     %     - Upward Solar Radiation Flux at Surface (uswrf) [W/m^2]      : for heat flux
%     %     - Sea level pressure (pres) [Pa]                              : for air pressure
%     %     - Air temperature at 2 m (air) [Kelvins]
%     %     - Relative Humidity on Pressure Levels (rhum) [%]
%     % The script also calculate the following parameters:
%     %     - Momentum flux (tau)
%     %     - Net solar radiation surface (nswrs = uswrf - dswrf)
%     %     - Net longwave radiation surface (nlwrs = ulwrf - dlwrf)
%     %     - Evaporation (Et)
%     %     - Precipitation-evaporation (P_E)
%     if develop_mode == 3
%     	fprintf('Loading Model objet file...\n')
%         load('varb/Mobj_03.mat');
%         load('varb/forcing_interp.mat');
%         fprintf('Done!\n');
%     else
%     	if develop_mode == 1
%         	% The script converts the NCEP data from the OPeNDAP server from longitudes 
%         	% in the 0 to 360 range to the latitudes in the -180 to 180 range. 
%         	% It also subsets for the right region (defined by Mobj.lon and Mobj.lat).
%         	% Uncomment the variables in get_NCEP_forcing as the varlist shows.
%         	fprintf('Downloading NCEP forcing from OPeNDAP server database...\n')
%         	forcing = get_NCEP_forcing(Mobj, inputConf.forceMJD, ...
%         	    'varlist', {'pevpr',...
%         	    'uwnd', 'vwnd',...
%         	    'uswrf', 'ulwrf', 'dswrf', 'dlwrf',...
%         	    'prate', 'lhtfl', 'shtfl', 'pres',...
%         	    'air', 'rhum'},...
%         	    'source', 'reanalysis2');
%         	forcing.domain_cols = length(forcing.lon);
%         	forcing.domain_rows = length(forcing.lat);
%         	if isfield(forcing, 'rhum')||isfield(forcing, 'pres')
%         	    forcing.domain_cols_alt = length(forcing.rhum.lon);
%         	    forcing.domain_rows_alt = length(forcing.rhum.lat);
%         	end
%         	% Convert the small subdomain into cartesian coordinates. We need
%         	% to do this twice because some of the NCEP data are on different
%         	% grids (e.g. sea level pressure, relative humidity etc.).
%         	tmpZone = regexpi(inputConf.utmZone,'\ ','split');
%         	[tmpLon, tmpLat] = meshgrid(forcing.lon, forcing.lat);
%         	[forcing.x, forcing.y] = wgs2utm(tmpLat(:), tmpLon(:), str2double(char(tmpZone{1}(1))), char(tmpZone{1}(2)));
%         	if isfield(forcing, 'rhum')||isfield(forcing, 'pres')
%         	    [tmpLon2, tmpLat2] = meshgrid(forcing.rhum.lon, forcing.rhum.lat);
%         	    [forcing.xalt, forcing.yalt] = wgs2utm(tmpLat2(:), tmpLon2(:), str2double(char(tmpZone{1}(1))), char(tmpZone{1}(2)));
%         	end
%         	clear tmpLon tmpLat tmpLon2 tmpLat2 tmpZone
%         	% Create arrays of the x and y positions.
%         	forcing.x = reshape(forcing.x, forcing.domain_rows, forcing.domain_cols);
%         	forcing.y = reshape(forcing.y, forcing.domain_rows, forcing.domain_cols);
%         	if isfield(forcing, 'rhum')||isfield(forcing, 'pres')
%         	    forcing.xalt = reshape(forcing.xalt, forcing.domain_rows_alt, forcing.domain_cols_alt);
%         	    forcing.yalt = reshape(forcing.yalt, forcing.domain_rows_alt, forcing.domain_cols_alt);
%         	end
%         	[forcing.lon, forcing.lat] = meshgrid(forcing.lon, forcing.lat);
%         	forcing = rmfield(forcing, {'domain_rows', 'domain_cols'});
%         	if isfield(forcing, 'rhum')||isfield(forcing, 'pres')
%         	    forcing = rmfield(forcing, {'domain_rows_alt', 'domain_cols_alt'});
%         	end
%         	fprintf('Saving NCEP forcing from OPeNDAP server database...')
%         	save(['../public/tokyobay/forcing','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat'],'forcing','-v7.3','-nocompression');
%         	fprintf('Done\n')
%         elseif develop_mode == 2
%     		fprintf('Loading NCEP forcing from the local database...');
%         	load(['../public/tokyobay/forcing','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat']);
%         	fprintf('Done!\n');
%     	end
%         % Interpolate the data onto the FVCOM unstructured grid.
%     	 
%     	save('varb/Mobj_03.mat','Mobj','-v7.3','-nocompression');
%         save('varb/forcing_interp.mat','forcing_interp','-v7.3','-nocompression');
%     end
%     % % adjust the values: temperature and wind speed
%     % % backup
%     % forcing_interp_backup = forcing_interp;
%     % % recover backup
%     % %forcing_interp = forcing_interp_backup;
%     % % adjust
%     % forcing_interp.uwnd.node = forcing_interp.uwnd.node * (1/3);
%     % forcing_interp.uwnd.data = forcing_interp.uwnd.data * (1/3);
%     % forcing_interp.vwnd.node = forcing_interp.vwnd.node * (1/3);
%     % forcing_interp.vwnd.data = forcing_interp.vwnd.data * (1/3);
%     % % adjust
%     % forcing_interp.air.node = forcing_interp.air.node + 10;
%     % forcing_interp.air.data = forcing_interp.air.data + 10;
%     % % adjust
%     % forcing_interp.nswrs.node = forcing_interp.nswrs.node * (1.40);
%     % forcing_interp.nswrs.data = forcing_interp.nswrs.data * (1.40);
%     % % adjust
%     % forcing_interp.nlwrs.node = forcing_interp.nlwrs.node * (1.00);
%     % forcing_interp.nlwrs.data = forcing_interp.nlwrs.data * (1.00);
%     % % adjust
%     % forcing_interp.lhtfl.node = forcing_interp.lhtfl.node * (0.50);
%     % forcing_interp.lhtfl.data = forcing_interp.lhtfl.data * (0.50);
%     % % adjust
%     % forcing_interp.shtfl.node = forcing_interp.shtfl.node * (0.00);
%     % forcing_interp.shtfl.data = forcing_interp.shtfl.data * (0.00);
% end
% 
% % Write out the surface forcing data to netCDF.
% fprintf('Writing meteorological forcing file...\n');
% if strcmpi(inputConf.doForcing, 'NCEP-CALCULATED')
%     forcing_interp_calculated.lon = forcing_interp.lon;
%     forcing_interp_calculated.lat = forcing_interp.lat;
%     forcing_interp_calculated.x = forcing_interp.x;
%     forcing_interp_calculated.y = forcing_interp.y;
%     forcing_interp_calculated.time = forcing_interp.time;
%     forcing_interp_calculated.uwnd = forcing_interp.uwnd;
%     forcing_interp_calculated.vwnd = forcing_interp.vwnd;
%     forcing_interp_calculated.prate = forcing_interp.prate;
%     forcing_interp_calculated.evap = forcing_interp.Et;
%     forcing_interp_calculated.nswrs = forcing_interp.nswrs;
%     forcing_interp_calculated.dlwrf = forcing_interp.dlwrf;
%     forcing_interp_calculated.air = forcing_interp.air;
%     forcing_interp_calculated.pres = forcing_interp.pres;
%     forcing_interp_calculated.rhum = forcing_interp.rhum;
%     plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.uwnd.node(1,:));
%     plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.vwnd.node(1,:));
%     % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.prate.node(1,:));
%     % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.evap.node(1,:));
%     plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.nswrs.node(1,:));
%     plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.dlwrf.node(1,:));
%     plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.air.node(1,:));
%     plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.pres.node(1,:));
%     plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.rhum.node(1,:));
%     write_FVCOM_forcing_calculated(Mobj, ...
%         fullfile(inputConf.outbase,inputConf.casename),...
%         forcing_interp_calculated,...
%         [inputConf.doForcing, 'atmospheric forcing data'],...
%         inputConf.FVCOM_version, ...
%         'floattime', true,...
%         'julian', true);
%     % fileprefix=fullfile(inputConf.outbase,inputConf.casename);
%     % data=forcing_interp_calculated;
%     % infos=[inputConf.doForcing, 'atmospheric forcing data'];
%     % fver=inputConf.FVCOM_version;
% end
% 
% % Write out the surface forcing data to netCDF.
% fprintf('Writing meteorological forcing file...\n');
% if strcmpi(inputConf.doForcing, 'NCEP')
%     write_FVCOM_forcing(Mobj, ...
%         fullfile(inputConf.outbase,inputConf.casename),...
%         forcing_interp,...
%         [inputConf.doForcing, 'atmospheric forcing data'],...
%         inputConf.FVCOM_version, ...
%         'floattime', true,...
%         'julian', true);
% %     fileprefix=fullfile(inputConf.outbase,inputConf.casename);
% %     data=forcing_interp;
% %     infos=[inputConf.doForcing, 'atmospheric forcing data'];
% %     fver=inputConf.FVCOM_version;
% end
% 
% if strcmpi(inputConf.doForcing, 'GWO')
%     Mobj.gwo.time = inputConf.forceMJD';
%     load(fullfile(basedir,'gwo\data\output\wnd.mat'));
%     Mobj.gwo.uwnd.node = [wnd(1,2);wnd(1:length(inputConf.forceMJD)-1,2)]';
%     Mobj.gwo.vwnd.node = [wnd(1,3);wnd(1:length(inputConf.forceMJD)-1,3)]';
%     Mobj.gwo.uwnd.data = [wnd(1,2);wnd(1:length(inputConf.forceMJD)-1,2)]';
%     Mobj.gwo.vwnd.data = [wnd(1,3);wnd(1:length(inputConf.forceMJD)-1,3)]';
%     load(fullfile(basedir,'gwo\data\output\nhf.mat'));
%     Mobj.gwo.hfx.node = [nhf(1,2);nhf(1:length(inputConf.forceMJD)-1,2)]';
%     Mobj.gwo.hfx.data = [nhf(1,2);nhf(1:length(inputConf.forceMJD)-1,2)]';
%     load(fullfile(basedir,'gwo\data\output\qr.mat'));
%     Mobj.gwo.nswrs.node = [qr(1,2);qr(1:length(inputConf.forceMJD)-1,2)]';
%     Mobj.gwo.nswrs.data = [qr(1,2);qr(1:length(inputConf.forceMJD)-1,2)]';
%     load(fullfile(basedir,'gwo\data\output\ql.mat'));
%     Mobj.gwo.nlwrs.node = -1*[ql(1,2);ql(1:length(inputConf.forceMJD)-1,2)]';
%     Mobj.gwo.nlwrs.data = -1*[ql(1,2);ql(1:length(inputConf.forceMJD)-1,2)]';
%     load(fullfile(basedir,'gwo\data\output\qh.mat'));
%     Mobj.gwo.shtfl.node = [qh(1,2);qh(1:length(inputConf.forceMJD)-1,2)]';
%     Mobj.gwo.shtfl.data = [qh(1,2);qh(1:length(inputConf.forceMJD)-1,2)]';
%     load(fullfile(basedir,'gwo\data\output\qe.mat'));
%     Mobj.gwo.lhtfl.node = [qe(1,2);qe(1:length(inputConf.forceMJD)-1,2)]';
%     Mobj.gwo.lhtfl.data = [qe(1,2);qe(1:length(inputConf.forceMJD)-1,2)]';
%     load(fullfile(basedir,'gwo\data\output\evp.mat'));
%     load(fullfile(basedir,'gwo\data\output\prate.mat'));
%     % convert kg/m2/s to m/s: 1 kg/m2/s = 0.001 m/s
%     Mobj.gwo.evap.node = [evp(1,2)*-0.001;evp(1:length(inputConf.forceMJD)-1,2)*-0.001]';
%     Mobj.gwo.evap.data = [evp(1,2)*-0.001;evp(1:length(inputConf.forceMJD)-1,2)*-0.001]';
%     Mobj.gwo.prate.node = [prate(1,2)*0.001;prate(1:length(inputConf.forceMJD)-1,2)*0.001]';
%     Mobj.gwo.prate.data = [prate(1,2)*0.001;prate(1:length(inputConf.forceMJD)-1,2)*0.001]';
%     load(fullfile(basedir,'gwo\data\output\pres.mat'));
%     Mobj.gwo.pres.node = 100*[pres(1,2);pres(1:length(inputConf.forceMJD)-1,2)]';
%     Mobj.gwo.pres.data = 100*[pres(1,2);pres(1:length(inputConf.forceMJD)-1,2)]';
%     % air 
%     load(fullfile(basedir,'gwo\data\output\air.mat'));
%     Mobj.gwo.air.node = [air(1,2);air(1:length(inputConf.forceMJD)-1,2)]';
%     Mobj.gwo.air.data = [air(1,2);air(1:length(inputConf.forceMJD)-1,2)]';
%     load(fullfile(basedir,'gwo\data\output\rhum.mat'));
%     Mobj.gwo.rhum.node = 100*[rhum(1,2);rhum(1:length(inputConf.forceMJD)-1,2)]';
%     Mobj.gwo.rhum.data = 100*[rhum(1,2);rhum(1:length(inputConf.forceMJD)-1,2)]';
%     %%%%
%     % cloud pe rhum sst vpres ghi
%     load(fullfile(basedir,'gwo\data\output\cloud.mat'));
%     Mobj.gwo.cloud = [cloud(1,2);cloud(1:length(inputConf.forceMJD)-1,2)]';
%     load(fullfile(basedir,'gwo\data\output\pe.mat'));
%     Mobj.gwo.pe = [pe(1,2);pe(1:length(inputConf.forceMJD)-1,2)]';
%     load(fullfile(basedir,'gwo\data\output\sst.mat'));
%     Mobj.gwo.sst = [sst(1,2);sst(1:length(inputConf.forceMJD)-1,2)]';
%     load(fullfile(basedir,'gwo\data\output\vpres.mat'));
%     Mobj.gwo.vpres = [vpres(1,2);vpres(1:length(inputConf.forceMJD)-1,2)]';
%     load(fullfile(basedir,'gwo\data\output\ghi.mat'));
%     Mobj.gwo.ghi = [ghi(1,2);ghi(1:length(inputConf.forceMJD)-1,2)]';
%     clear q* wnd nhf prate evp pres air cloud pe rhum sst vpres ghi
%     save('varb/Mobj_03.mat','Mobj','-v7.3','-nocompression');
% end
% 
% if strcmpi(inputConf.doForcing, 'GWO')
%     write_FVCOM_gwo_forcing(Mobj, ...
%         fullfile(inputConf.outbase,[inputConf.casename,'_gwo']),...
%         [inputConf.doForcing, ' atmospheric forcing data'],...
%         inputConf.FVCOM_version, ...
%         'floattime', true,...
%         'julian', true,...
%         'time dependent', true);
% end
% fprintf('Writing meteorological forcing file...done!\n');
% 
% 
% % uwnd and vwnd
% subplot(2,1,1);
% plot(forcing_interp.time,forcing_interp.uwnd.node(1,:));
% title('NCEP: uwind (m/s)');
% subplot(2,1,2);
% plot(Mobj.gwo.time,Mobj.gwo.uwnd.node);
% title('GWO: uwind (m/s)');
% % uwnd and vwnd
% subplot(2,1,1);
% plot(forcing_interp.time,forcing_interp.vwnd.node(1,:));
% title('NCEP: vwind (m/s)');
% subplot(2,1,2);
% plot(Mobj.gwo.time,Mobj.gwo.vwnd.node);
% title('GWO: vwind (m/s)');
% % nswrs
% subplot(2,1,1);
% plot(forcing_interp.time,forcing_interp.nswrs.node(1,:));
% title('NCEP: Net shortwave radiation surface (nswrs, W/m2)');
% subplot(2,1,2);
% plot(Mobj.gwo.time,Mobj.gwo.nswrs.node);
% title('GWO: ?????????? Qr (nswrs, W/m2)');
% % nlwrs
% subplot(2,1,1);
% plot(forcing_interp.time,forcing_interp.nlwrs.node(1,:));
% title('NCEP: Net longwave radiation surface (nlwrs, W/m2)');
% subplot(2,1,2);
% plot(Mobj.gwo.time,Mobj.gwo.nlwrs.node);
% title('GWO: ?????????? Ql (nlwrs, W/m2)');
% % shtfl
% subplot(2,1,1);
% plot(forcing_interp.time,forcing_interp.shtfl.node(1,:));
% title('NCEP: Sensible heat flux (shtfl, W/m2)');
% subplot(2,1,2);
% plot(Mobj.gwo.time,Mobj.gwo.shtfl.node);
% title('GWO: ????????????? Qh (shtfl, W/m2)');
% % lhtfl
% subplot(2,1,1);
% plot(forcing_interp.time,forcing_interp.lhtfl.node(1,:));
% title('NCEP: Latent heat flux (lhtfl, W/m2)');
% subplot(2,1,2);
% plot(Mobj.gwo.time,Mobj.gwo.lhtfl.node);
% title('GWO: ????????????? Qe (lhtfl, W/m2)');
% % nshf = nlwrs + nswrs - lhtfl - shtfl
% % QT = Qr + Ql - Qh - Qe
% subplot(2,1,1);
% plot(forcing_interp.time,forcing_interp.nswrs.node(1,:)+forcing_interp.nlwrs.node(1,:)-forcing_interp.shtfl.node(1,:)-forcing_interp.lhtfl.node(1,:));
% title('NCEP: Surface net heat flux (nshf, W/m2)');
% subplot(2,1,2);
% plot(Mobj.gwo.time,Mobj.gwo.nswrs.node+Mobj.gwo.nlwrs.node-Mobj.gwo.shtfl.node-Mobj.gwo.lhtfl.node);
% title('GWO: ?????? QT (hfx, W/m2)');
% % air
% subplot(2,1,1);
% plot(forcing_interp.time,forcing_interp.air.node(1,:));
% title('NCEP: Air temperature (air, degC)');
% subplot(2,1,2);
% plot(Mobj.gwo.time,Mobj.gwo.air);
% title('GWO: Air temperature (air, degC)');
% % pres
% subplot(2,1,1);
% plot(forcing_interp.time,forcing_interp.pres.node(1,:));
% title('NCEP: Air pressure (pres, unit)');
% subplot(2,1,2);
% plot(Mobj.gwo.time,Mobj.gwo.pres.node);
% title('GWO: Air pressure (pres, unit)');
% % rhum
% subplot(2,1,1);
% plot(forcing_interp.time,forcing_interp.rhum.node(1,:));
% title('NCEP: Relative humudity (rhum, %)');
% subplot(2,1,2);
% plot(Mobj.gwo.time,Mobj.gwo.rhum);
% title('GWO: Relative humudity (rhum, unit)');
% % prate
% subplot(2,1,1);
% plot(forcing_interp.time,(forcing_interp.prate.node(1,:)));
% title('NCEP: prate (prate, m/s)');
% subplot(2,1,2);
% plot(Mobj.gwo.time,Mobj.gwo.prate.node);
% title('GWO: prate (prate, kg/m2/s)');
% % % evap
% % subplot(2,1,1);
% % plot(forcing_interp.time,(forcing_interp.Et.node(1,:)));
% % title('NCEP: Evap (evap, m/s)');
% % subplot(2,1,2);
% % plot(Mobj.gwo.time,Mobj.gwo.evap.node);
% % title('GWO: Evap (evap, kg/m2/s)');
% % % P-E
% % subplot(2,1,1);
% % plot(forcing_interp.time,(forcing_interp.prate.node(1,:)+forcing_interp.Et.node(1,:)));
% % title('NCEP: P-E (P-E, m/s)');
% % subplot(2,1,2);
% % plot(Mobj.gwo.time,Mobj.gwo.prate.node+Mobj.gwo.evap.node);
% % title('GWO: P-E (P-E, kg/m2/s)');
%%
%%%------------------------------------------------------------------------
%%%                     River discharge and output
%%%------------------------------------------------------------------------
% FVCOM river.
if strcmpi(inputConf.riverForcing, 'FLUX')
    fprintf('Reading river forcing file...');
    Mobj = get_COSTUM_river_location(inputConf, Mobj);
    Mobj = get_COSTUM_river_variable(Mobj, ...
        {inputConf.river.flux,inputConf.river.temp,inputConf.river.salt}, ...
        {'flux','temp','salt'},...
        inputConf.river.infos,...
        'time', true);
    Mobj.nRivers = Mobj.river.number;
    inc = 1;
    Mobj.rivermouth = cell(1); % don't preallocate as we don't know how many we'll have
    for s = 1:Mobj.nRivers
        [node, ~] = find_nearest_pt(Mobj.river.location(s, 3), Mobj.river.location(s, 4), Mobj);
        [~, elem] = min(abs(sqrt((Mobj.xc - Mobj.river.location(s, 3)).^2 + Mobj.yc - Mobj.river.location(s, 4)).^2));
        Mobj.rivermouth{inc} = {inc, Mobj.river.location(s, 3), Mobj.river.location(s, 4), node, Mobj.h(node), Mobj.river.name(s), elem};
        riverList(s,1) = Mobj.rivermouth{inc}(1,4);
        Mobj.riverList = cell2mat(riverList);
        inc = inc + 1;
    end
    clear node elem inc s
    Mobj = add_river_nodes_list(Mobj,Mobj.riverList,Mobj.river.name,1);
    fprintf('Done!\n');
end

% river file
fprintf('Writing river forcing file...');
write_FVCOM_river(fullfile(inputConf.outbase,...
    [inputConf.casename,'_river.nc']),...
     inputConf.river.infos,...
     Mobj.river.timeMJD,...
     Mobj.river.flux*0.1,...
     Mobj.river.temp,...
     Mobj.river.salt,...
    'Tokyo Bay rivers',...
    'Model river boundary input');
write_FVCOM_river_nml(Mobj, ...
    fullfile(inputConf.outbase,'rivers_namelist.nml'), ...
    [inputConf.casename,'_river.nc'],...
    '''uniform''');
fprintf('Done!\n');
save('varb/Mobj_04.mat','Mobj','-v7.3','-nocompression');

%%
% clear ans tpxo_dir
% fprintf('All done!\n')
% diary off;
% 
% %%
% %%%------------------------------------------------------------------------
% %%%                     NML file and output
% %%%------------------------------------------------------------------------
% % FVCOM running setup in nml file
% % these cannot be left to the default values
% inputConf.estimate = min(Mobj.ts);
% inputConf.timestep = floor(inputConf.estimate) - 2;
% inputConf.EXTSTEP_SECONDS = inputConf.timestep;
% inputConf.ramp = 1.0; %day
% inputConf.isplit = 10;
% inputConf.IRAMP = floor(inputConf.ramp*24*60*60/(inputConf.EXTSTEP_SECONDS*inputConf.isplit)); % ramp over one day
% inputConf.START_DATE=datestr(inputConf.startDate,'yyyy-mm-dd HH:MM:SS');%           '2046-02-01 00:00:00';
% inputConf.END_DATE=datestr(inputConf.endDate,'yyyy-mm-dd HH:MM:SS');% '2046-03-01 00:00:00';
% inputConf.RST_FIRST_OUT=inputConf.START_DATE;
% % Change sigma file
% inputConf.report = 6/24;
% inputConf.IREPORT = floor(inputConf.report*24*60*60/(inputConf.timestep*inputConf.isplit));
% inputConf.NC_FIRST_OUT=inputConf.START_DATE;
% inputConf.NCAV_FIRST_OUT=inputConf.START_DATE;
% % NC file
% inputConf.NC_FIRST_OUT=inputConf.START_DATE;
% inputConf.NCAV_FIRST_OUT=inputConf.START_DATE;
% inputConf.PROJECTION_REFERENCE='+proj=utm +zone=54 +ellps=bessel +units=m +no_defs';
% inputConf.TS_nudge = inputConf.dateobs;
% if isfield(inputConf,'TS_nudge')
%     inputConf.OBC_TEMP_NUDGING_TIMESCALE = 1/(inputConf.TS_nudge*3600/inputConf.timestep);
%     inputConf.OBC_SALT_NUDGING_TIMESCALE = 1/(inputConf.TS_nudge*3600/inputConf.timestep);
% end
% conf = inputConf;
% [fmt, nml] = make_default_nml(inputConf);
% res = write_model_nml(inputConf, nml, fmt);
% %% Mesh
% % Plot the mesh, nodes and water depth.
% x = Mobj.lon;
% y = Mobj.lat;
% nodeList = double(cell2mat(Mobj.read_obc_nodes(1)));
% plot_MESH(01, [Mobj.lon, Mobj.lat], Mobj.tri, Mobj.h, ...
%     x(nodeList), y(nodeList), x(Mobj.riverList), y(Mobj.riverList));
% saveas(gcf,'plot/00_mesh.png')
% clear nodeList x y
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % %% TS surface and profiles
% % % Plot TS at open boundary
% % nn = 20;   % open boundary index
% % tt = 1;    % time index
% % fvz = 1;   % fvcom depth index (currently 1-20)
% % hyz = 1;   % coarse depth index (1-33)
% % % Open boundary index.
% % oNodes = [Mobj.read_obc_nodes{:}];
% % % Number of sigma layers.
% % fz = size(Mobj.siglayz, 2);
% % fields = fieldnames(hycom);
% % % Find the first 4D array and use it to get the number of vertical levels
% % % and time steps.
% % for ff = 1:length(fields)
% %     if isfield(hycom.(fields{ff}), 'data') && ndims(hycom.(fields{ff}).data) > 3
% %         [nx, ny, nz, nt] = size(hycom.(fields{ff}).data);
% %         break
% %     end
% % end
% % hdepth = permute(repmat(hycom.Depth.data, [1, nx, ny]), [2, 3, 1]);
% % % Find the coarse seabed indices
% % % [~, hyz] = nanmax(hdepth, [], 3);
% % % Use the existing rectangular arrays for the nearest point lookup.
% % [lon, lat] = deal(hycom.lon, hycom.lat);
% % fvlon = Mobj.lon(oNodes);
% % fvlat = Mobj.lat(oNodes);
% % % Get the corresponding indices for the coarse data
% % [~, idx] = min(sqrt((lon(:) - fvlon(nn)).^2 + (lat(:) - fvlat(nn)).^2));
% % [xidx, yidx] = ind2sub(size(lon), idx);
% % zidx = isfinite(hdepth(xidx, yidx, :));
% % hz = 1:nz;
% % dx = mean(diff(hycom.lon(:)));
% % dy = mean(diff(hycom.lat(:)));
% % vartoplot = {'temperature','salinity'};
% % Mobj.salinity = Mobj.salt;
% % % Plot the surface temperature and salinity
% % for i =1:2
% %     var = hycom.(vartoplot{i}).data(:, :, :, tt);
% %     plot_TS_surface(20+i, ...
% %         hycom.lon - (dx / 2), ...
% %         hycom.lat - (dy / 2), ...
% %         squeeze(var(:, :, hyz)), ...
% %         Mobj.lon(oNodes), Mobj.lat(oNodes), 40, Mobj.(vartoplot{i})(:, fvz, tt), ...
% %         lon(xidx, yidx), lat(xidx, yidx), ...
% %         Mobj.lon(oNodes(nn)), Mobj.lat(oNodes(nn)),...
% %         vartoplot{i});
% %     name = ['plot_02_TS_surface_',vartoplot{i},'.png'];
% %     saveas(gcf, name)
% % end
% % % To check the intepolation quality
% % for i =1:2
% %     plot_TS_profile(22+i, ...
% %         Mobj.(vartoplot{i})(nn, :, tt), Mobj.siglayz(oNodes(nn), :), ...
% %         squeeze(hycom.(vartoplot{i}).data(xidx, yidx, zidx, tt)), squeeze(-hdepth(xidx, yidx, zidx)), ...
% %         Mobj.(vartoplot{i})(nn, :, tt), 1:fz, ...
% %         squeeze(hycom.(vartoplot{i}).data(xidx, yidx, zidx, tt)), hz(zidx),...
% %         vartoplot{i});
% %     name = ['plot_02_TS_profile_',vartoplot{i},'.png'];
% %     saveas(gcf, name)
% % end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % %% Animation of TS surface
% % vartoplot = {'temperature','salinity'};
% % Num = 10;
% % plot_nc_layer_TS(Num, ncfile, 4, [1], vartoplot);
% % 
% % %% Animation of TS open boundary
% % vartoplot = {'temperature','salinity'};
% % Num = 20;
% % plot_obj_obc_TS(Num, Mobj, vartoplot);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% inputConf.sst_pattern = 'MUR-JPL-L4-GLOB-v4.1.nc';
% inputConf.sst_dir = 'C:\Users\Yulong WANG\Documents\GitHub\MUR-JPL-L4-GLOB-v4.1';
% i = 2003;
% inputConf.year=i;
% Mobj = interp_sst_assimilation(Mobj, inputConf, ['sst_',num2str(i),'.nc']);
% for i=2004:2018
%     inputConf.year=i;
%     Mobj = interp_sst_assimilation_update(Mobj, inputConf, ['sst_',num2str(i),'.nc']);
% end
% % conf=inputConf;
% % output_file='sst.nc'