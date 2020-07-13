function nc = loadNetCDF(ncfile,nc)

%     nc.kspe_dye        = size(nc.k_specify,1);   % Number of sigma layer for specify dye release
%     nc.mspe_dye        = size(nc.m_specify,1);   % Number of node for specify dye release
% 
%     nc.dti       = ncread(ncfile,'dti');
% 
%     nc.i_obc_n   = ncread(ncfile,'obc_nodes');          % [OK][nobc]:   Open boundary node list for fvcom
%     nc.ntrg      = ncread(ncfile,'ntrg');               % [OK][ncv]:    Element associated with this control volume edge
%     nc.dt1       = ncread(ncfile,'dt1');                % [OK][nt]:     Depth at previous time step
%     nc.dltxe     = ncread(ncfile,'dltxe');              % [OK][ncv]:    X length of nodal control volume edges
%     nc.dltye     = ncread(ncfile,'dltye');              % [OK][ncv]:    Y length of nodal control volume edges
%     nc.uout      = ncread(ncfile,'uout');               % [NC][nt,kb]:  X velocity
%     nc.vout      = ncread(ncfile,'vout');               % [NC][nt,kb]:  Y velocity
%     nc.dz1       = ncread(ncfile,'dz1');                % [OK][nt,kb]:  Delta-sigma value
%     nc.ntsn      = ncread(ncfile,'ntsn');               % [NC][m]:      Number of nodes surrounding each node
%     nc.nbsn      = ncread(ncfile,'nbsn');               % [NC][m,8]:    Indices of nodes surrounding each node
%     nc.dltxtrie  = ncread(ncfile,'dltxtrie');           % [OK][m,kb+1]: Delta x triangle edge
%     nc.dltytrie  = ncread(ncfile,'dltytrie');           % [OK][m,kb+1]: Delta y triangle edge
    nc.dye       = ncread(ncfile,'DYE');                % [NC][m,kb+1]: Dye concentration at node
    nc.dye_age   = ncread(ncfile,'DYE_AGE');            % [NC][m,kb+1]: Dye concentration at node
%     nc.art2      = ncread(ncfile,'art2');               % [NC][m]:      Area of elements around node
%     nc.nn_hvc    = ncread(ncfile,'nn_hvc');             % [OK][m]:      Variable horizontal viscosity coefficents
%     nc.dltxncve  = ncread(ncfile,'dltxncve');           % [OK][ncv,2]:  Delta x node to control volume edge
%     nc.dltyncve  = ncread(ncfile,'dltyncve');           % [OK][ncv,2]:  Delta y node to control volume edge
%     nc.viscofh   = ncread(ncfile,'viscofh');            % [NC][m,kb]:   Horizontal Turbulent Eddy Viscosity For Scalars
%     nc.niec      = ncread(ncfile,'niec');               % [OK][ncv,2]:  Node indices of each volume edge
%     nc.iswetn    = ncread(ncfile,'wet_nodes');          % [NC][m]:      (wet_nodes) Node porosity at nodes for time n
%     nc.iswetnt   = ncread(ncfile,'wet_nodes_prev_int'); % [NC][m]:      (wet_nodes_prev_int) Node porosity at nodes for time n-1 internal
%     nc.wts       = ncread(ncfile,'omega');              % [NC][m,kb+1]: (omega) Vertical velocity in sigma system
%     nc.dz        = ncread(ncfile,'dz');                 % [OK][m,kb]:   Delta-sigma value
%     nc.art1      = ncread(ncfile,'art1');               % [NC][m]:      Area of node-base control volume
%     nc.isonb     = ncread(ncfile,'isonb');              % [OK][m]:      Node maker = 0,1,2
%     nc.dtfa      = ncread(ncfile,'dtfa');               % [OK][m]:      Adjusted depth for mass conservation
%     nc.dt        = ncread(ncfile,'dt');                 % [NC][m]:      Depth at previous time step
% 
%     nc.kh        = ncread(ncfile,'kh');                 % [NC][m,kb]:   Turbulent diffusivity for salinity/temp
%     nc.d         = ncread(ncfile,'h');                  % [NC][m]:      Current depth
%     nc.dzz       = ncread(ncfile,'dzz');                % [OK][m,kbm2]: Delta of intra level sigma
%     nc.sita_gd   = ncread(ncfile,'sita_gd');            % [OK][m]:      Bottom depth gradients？
%     nc.ah_bottom = ncread(ncfile,'ah_bottom');          % [OK][m]
%     nc.phpn      = ncread(ncfile,'phpn');               % [OK][m]
%     nc.next_obc  = ncread(ncfile,'next_obc');           % [OK][nobc]: Interior neighbor of open boundary node
%     nc.uard_obcn = ncread(ncfile,'uard_obcn');          % [OK][nobc]: Nonlinear Velocity Open Boundary Condition Arrays
    nc.time      = ncread(ncfile,'time');
    nc.timestart = nc.time(1);
    nc.timestop  = nc.time(end);

%     nc.kb        = size(nc.kh,2);
%     nc.kbm1      = nc.kb-1;
%     nc.kbm2      = nc.kb-2;
%     nc.m         = size(nc.kh,1);
%     nc.ncv       = size(nc.ntrg,1);                     % NUMBER OF INTERNAL CONTROL VOLUMES (EXTENDED LOCAL ONLY)
%     nc.ncv_i     = nc.ncv;
%     nc.iobcn     = size(nc.i_obc_n,1);                  % Local number of open boundary nodes for fvcom
%     nc.nlid      = [1:size(nc.kh,1)];                   % [OK][m]
    
    nc.lon       = ncread(ncfile,'lon');
    nc.lon       = double(nc.lon);
    nc.lat       = ncread(ncfile,'lat');
    nc.lat       = double(nc.lat);
    nc.x         = ncread(ncfile,'x');
    nc.x         = double(nc.x);
    nc.y         = ncread(ncfile,'y');
    nc.y         = double(nc.y);
    nc.nv        = ncread(ncfile,'nv');
    nc.nv        = double(nc.nv);
    
    % make matrix
    nc.xy        = [nc.x, nc.y];
    nc.lonlat    = [nc.lon, nc.lat];

    nc.xy        = double(nc.xy);
    nc.lonlat    = double(nc.lonlat);
    
    % make trianglua mesh
    nc.tri_xy    = triangulation(nc.nv,nc.xy);
    nc.tri_lonlat= triangulation(nc.nv,nc.lonlat);
    
    nc.tri       = nc.tri_xy.ConnectivityList;
    
    nc.tri(:,[2 3]) =nc.tri(:,[3 2]);
    
%     nc.time = ncread(ncfile,'time');
    nc.Times = ncread(ncfile,'Times');
    nc.Times = nc.Times';

end