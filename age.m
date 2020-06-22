function dyef_age = age(nc,dye_age,inttime,c)

%%
    % INPUT
    % nc        : Parameters from FVCOM calculation.
    % dye       : Dye concentration at nodes. [m,kb]
    % inttime   : Current time.
    
    % OUTPUT
    % dyef      : Dye concentration from previous time. [m,kb]
    
    % 在平面二维方向上使用有限体积法
    % 垂向计算使用的是有限差分法
    % Use the FVM in the 2D plane
    % Uses the FDM in vertical plane
    
    dye_age_source_term = 0;
    % dye_source_term = nc.dye_source_term;
    umol            = nc.umol;         % Vertical mixing coefficient (1E-5)    
    fact            = nc.fact;
    fm1             = nc.fm1;

    k_specify       = nc.k_specify;    % NO. of sigma layer for specify dye release
    m_specify       = nc.m_specify;    % NO. of node for specify dye release
    kspe_dye        = nc.kspe_dye;     % Number of sigma layer for specify dye release
    mspe_dye        = nc.mspe_dye;     % Number of node for specify dye release

    dti       = nc.dti;

    i_obc_n   = nc.i_obc_n;            % [OK][nobc]:   Open boundary node list for fvcom
    ntrg      = nc.ntrg;               % [OK][ncv]:    Element associated with this control volume edge
    dt1       = nc.dt1;                % [OK][nt]:     Depth at previous time step
    dltxe     = nc.dltxe;              % [OK][ncv]:    X length of nodal control volume edges
    dltye     = nc.dltye;              % [OK][ncv]:    Y length of nodal control volume edges
    u         = nc.u;                  % [NC][nt,kb]:  X velocity
    v         = nc.v;                  % [NC][nt,kb]:  Y velocity
    dz1       = nc.dz1;                % [OK][nt,kb]:  Delta-sigma value
    ntsn      = nc.ntsn;               % [NC][m]:      Number of nodes surrounding each node
    nbsn      = nc.nbsn;               % [NC][m,8]:    Indices of nodes surrounding each node
    dltxtrie  = nc.dltxtrie;           % [OK][m,kb+1]: Delta x triangle edge
    dltytrie  = nc.dltytrie;           % [OK][m,kb+1]: Delta y triangle edge
    art2      = nc.art2;               % [NC][m]:      Area of elements around node
    nn_hvc    = nc.nn_hvc;             % [OK][m]:      Variable horizontal viscosity coefficents
    dltxncve  = nc.dltxncve;           % [OK][ncv,2]:  Delta x node to control volume edge
    dltyncve  = nc.dltyncve;           % [OK][ncv,2]:  Delta y node to control volume edge
    viscofh   = nc.viscofh;            % [NC][m,kb]:   Horizontal Turbulent Eddy Viscosity For Scalars
    niec      = nc.niec;               % [OK][ncv,2]:  Node indices of each volume edge
    iswetn    = nc.iswetn;             % [NC][m]:      (wet_nodes) Node porosity at nodes for time n
    iswetnt   = nc.iswetnt;            % [NC][m]:      (wet_nodes_prev_int) Node porosity at nodes for time n-1 internal
    wts       = nc.wts;                % [NC][m,kb+1]: (omega) Vertical velocity in sigma system
    dz        = nc.dz;                 % [OK][m,kb]:   Delta-sigma value
    art1      = nc.art1;               % [NC][m]:      Area of node-base control volume
    isonb     = nc.isonb;              % [OK][m]:      Node maker = 0,1,2
    dtfa      = nc.dtfa;               % [OK][m]:      Adjusted depth for mass conservation
    dt        = nc.dt;                 % [NC][m]:      Depth at previous time step

    kh        = nc.kh;                 % [NC][m,kb]:   Turbulent diffusivity for salinity/temp
    d         = nc.d;                  % [NC][m]:      Current depth
    dzz       = nc.dzz;                % [OK][m,kbm2]: Delta of intra level sigma
    sita_gd   = nc.sita_gd;            % [OK][m]:      Bottom depth gradients？
    ah_bottom = nc.ah_bottom;          % [OK][m]
    phpn      = nc.phpn;               % [OK][m]
    next_obc  = nc.next_obc;           % [OK][nobc]: Interior neighbor of open boundary node
    uard_obcn = nc.uard_obcn;          % [OK][nobc]: Nonlinear Velocity Open Boundary Condition Arrays
    dyestart  = nc.dyestart;
    dyestop   = nc.dyestop;

    kb        = size(kh,2);
    kbm1      = kb-1;
    kbm2      = kb-2;
    m         = size(kh,1);
    ncv       = size(ntrg,1);          % NUMBER OF INTERNAL CONTROL VOLUMES (EXTENDED LOCAL ONLY)
    ncv_i     = ncv;
    iobcn     = size(i_obc_n,1);       % Local number of open boundary nodes for fvcom
    nlid      = [1:size(kh,1)];        % [OK][m]

    %% SUBROUTINE ADV_DYE 
    %-------------------------------------------------------
    % Initialize Fluxes and dyef
    %-------------------------------------------------------
    xflux     = zeros(m,kbm1);
    xflux_adv = zeros(m,kbm1);
    dyef_age      = zeros(m,kbm1);

    %-------------------------------------------------------
    % Loop Over Control Volume Sub-Edges And 
    % Calculate Normal Velocity
    %-------------------------------------------------------
    % FVCOM计算时候也不是按照控制体数量进行循环，而是按照控制边的数量循环。
    % 控制体边上法相速度计算。
    % Loop is not  the on the number of control volume.
    % Loop is on the number of control volume edges.
    for i=1:ncv
      i1=ntrg(i);
      for k=1:kbm1
          dtij(i,k)=dt1(i1)*dz1(i1,k);
          uvn(i,k) = v(i1,k)*dltxe(i) - u(i1,k)*dltye(i);
      end
    end

    %-------------------------------------------------------
    % Calculate the Advection and Horizontal Diffusion Terms
    %-------------------------------------------------------

    % 按照分层水体循环。
    % Sigma loop.
    for k=1:kbm1
       pdyepx  = zeros(m,1);
       pdyepy  = zeros(m,1);
       pdyepxd = zeros(m,1);
       pdyepyd = zeros(m,1);
       % 按照节点循环。
       % node loop.
       for i=1:m
           % 按照节点为中心的控制体上所有边循环。
           % i1 和 i2 是边的两个端点。
           % Control volume edge loop.
           % i1 and i2 are two points of an edge.
           for j=1:ntsn(i)-1
               i1=nbsn(i,j);
               i2=nbsn(i,j+1);
               % 根据干湿条件确定 ffd 和 ff1。
               % iswetn == 0 表示为干。
               % ffd 和 ff1 是控制体边上中心点的浓度，或平均浓度。
               % Calculate ffd and ff1 based on wet/dry condition.
               % iswetn == 0 means dry.
               % ffd and ff1 is the central-point concentration.
               if iswetn(i1) == 0 && iswetn(i2) == 1
                   ffd=0.5*(dye_age(i,k)+dye_age(i2,k));
                   ff1=0.5*(dye_age(i,k)+dye_age(i2,k));
               elseif iswetn(i1) == 1 && iswetn(i2) == 0
                   ffd=0.5*(dye_age(i1,k)+dye_age(i,k));
                   ff1=0.5*(dye_age(i1,k)+dye_age(i,k));
               elseif iswetn(i1) == 0 && iswetn(i2) == 0
                   ffd=0.5*(dye_age(i,k)+dye_age(i,k));
                   ff1=0.5*(dye_age(i,k)+dye_age(i,k));
               else
                   ffd=0.5*(dye_age(i1,k)+dye_age(i2,k));
                   ff1=0.5*(dye_age(i1,k)+dye_age(i2,k));
               end
               % 一阶偏导计算也采用格林公式将面积分降维化为线积分计算。
               % 平面积分所采用的控制体和上面不同，采用和节点相邻的所有三角形单元计算。
               % 对于边界点来说，只取相邻的两个单元。
               % The 1st-order PDE uses the Green's formula to reduce the dimensionality,
               % to reduce the area to linear integral calculations.
               % This use all neibour trianle elements around the node.
               pdyepx(i) = pdyepx(i) + ff1 * dltytrie(i,j);
               pdyepy(i) = pdyepy(i) + ff1 * dltxtrie(i,j);
               pdyepxd(i)= pdyepxd(i)+ ffd * dltytrie(i,j);
               pdyepyd(i)= pdyepyd(i)+ ffd * dltxtrie(i,j);
           end
           % 将相邻控制体上的线积分累加为一阶偏导到节点上。
           % gather all neighboring control volumes connecting at node.
           % ∫cidy = (∂ci/∂x)*Ai, ∫-cidx = (∂ci/∂y)*Ai
           % ∂ci/∂x = (∫cidy)/Ai, ∂ci/∂y = (∫-cidx)/Ai
           % ∂ci/∂x 和 ∂ci/∂y 为节点周围所有三角形以内的平均值。
           pdyepx(i)  =pdyepx(i)/art2(i); % 
           pdyepy(i)  =pdyepy(i)/art2(i);
           pdyepxd(i) =pdyepxd(i)/art2(i);
           pdyepyd(i) =pdyepyd(i)/art2(i);
       end

       if k == kbm1
           for i=1:m
               pfpxb(i) = pdyepx(i);
               pfpyb(i) = pdyepy(i);
           end
       end

       for i=1:m
           viscoff(i)=viscofh(i,k);
       end
      
       % Sea bottom
       if k == kbm1
           ah_bottom(1:m) = (fact*viscoff(1:m) + fm1) * nn_hvc(1:m);
       end
       
       % control volume loop.
       for i=1:ncv_i
           ia=niec(i,1);
           ib=niec(i,2);
           
           % 采用迎风格式。其中，控制边上 fij1 和 fij2 计算采用泰勒公式2阶精度。
           % 对于水平对流项中控制边上Ci的二阶精度计算，只需按照泰勒公式计算即可。
           % Horizontal convection term,
           % The 2nd-order PDE of Ci on the control volume edge, 
           % it is okay to calculate according to the Taylor formula.
           % ci1 = ci,a + ∆x,nc * ∂ci,a/∂x,a + ∆y,nc * ∂ci,a/∂y,a
           % ci2 = ci,b + ∆x,nc * ∂ci,b/∂x,b + ∆y,nc * ∂ci,b/∂y,b
           fij1=dye_age(ia,k) + dltxncve(i,1)*pdyepx(ia) + dltyncve(i,1)*pdyepy(ia);
           fij2=dye_age(ib,k) + dltxncve(i,1)*pdyepx(ib) + dltyncve(i,1)*pdyepy(ib);

           s1min=min(dye_age(nbsn(ia,1:ntsn(ia)-1),k));
           s1min=min(s1min, dye_age(ia,k));
           s1max=max(dye_age(nbsn(ia,1:ntsn(ia)-1),k));
           s1max=max(s1max, dye_age(ia,k));
           s2min=min(dye_age(nbsn(ib,1:ntsn(ib)-1),k));
           s2min=min(s2min, dye_age(ib,k));
           s2max=max(dye_age(nbsn(ib,1:ntsn(ib)-1),k));
           s2max=max(s2max, dye_age(ib,k));

           % if u1 < 0, c1,a = s1min, Upwind。
           % if u1 > 0, c1,a = s1max, Upwind。
           % if u2 < 0, c2,a = s2min, Upwind。
           % if u2 > 0, c2,a = s2max, Upwind。
           if fij1 < s1min, fij1=s1min; end
           if fij1 > s1max, fij1=s1max; end
           if fij2 < s2min, fij2=s2min; end
           if fij2 > s2max, fij2=s2max; end

           un = uvn(i,k);

           % 其中 viscof 为水平扩散系数，一阶偏导采用控制边相邻两个节点平均值。
           % viscof is horizontal diffusion coefficient.
           % 1st order average.
           viscof=(...
               fact*0.5*(viscoff(ia)*nn_hvc(ia)+viscoff(ib)*nn_hvc(ib)) + ...
               fm1*0.5*(nn_hvc(ia)+nn_hvc(ib))...
               );
           
           % fxx 和 fyy 中还包含了水深 dtij，后面会在计算流量时除去。
           % fxx and fyy are horizontal diffusion term for x and y.
           % but the following term contain water depth dtij.
           % this will be divided later after updating dye.
           txx=0.5*(pdyepxd(ia)+pdyepxd(ib))*viscof;
           tyy=0.5*(pdyepyd(ia)+pdyepyd(ib))*viscof;

           fxx=-dtij(i,k)*txx*dltye(i);
           fyy= dtij(i,k)*tyy*dltxe(i);

           % 水平对流项离散：ci*(udy-vdx)
           % Horizontal advection term: ci*(udy-vdx)
           % ci = ci,0 + (∂ci/∂x)*∆x + (∂ci/∂y)*∆y
           % ci 在控制边上计算二阶精度。
           % ci are 2nd order value of control volume edge.
           % 其中，ci,0 的选取由距离控制边最近的节点获取。
           % and ci,0 is calculated from nearest one nodes.
           % 水平扩散项离散：
           % descrete ci*(udy-vdx)
           % un           : normal velocity
           % dtij         : water depth
           % fij1 and fij2: ci at one nearest nodes
           % -un*dtij(i,k)*((1.0+sign(1.0,un))*fij2+(1.0-sign(1.0,un))*fij1)*0.5
           
           % exflux=-un*dtij(i,k)*((1.0+sign(1.0,un))*fij2+...
           %     (1.0-sign(1.0,un))*fij1)*0.5+...
           %     fxx+fyy;
           % did not find sign in matlab, use if.
           % 没有找到替换 sign 函数，于是改写如下。
           if un >= 0
               exflux=-un*dtij(i,k)*...
                   ((1.0+1.0)*fij2+(1.0-1.0)*fij1)*0.5+...
                   fxx+fyy;
           else
               exflux=-un*dtij(i,k)*...
                   ((1.0-1.0)*fij2+(1.0+1.0)*fij1)*0.5+...
                   fxx+fyy;
           end
           
           xflux(ia,k)=xflux(ia,k) + exflux;
           xflux(ib,k)=xflux(ib,k) - exflux;
           xflux_adv(ia,k)=xflux_adv(ia,k) + (exflux-fxx-fyy);
           xflux_adv(ib,k)=xflux_adv(ib,k) - (exflux-fxx-fyy);
       end % control volume loop end.
    end % sigma loop end.

    %-------------------------------------------------------
    % Accumulate Fluxes at Boundary Nodes
    %-------------------------------------------------------

    xflux_obc = zeros(iobcn,kbm1);
    for k=1:kbm1
     if iobcn > 0
         for i=1:iobcn
             i1=i_obc_n(i);
             xflux_obc(i,k)=xflux_adv(i1,k);
         end
     end
    end

    %-------------------------------------------------------
    % The central difference scheme in vertical advection？
    %-------------------------------------------------------
    
    for i=1:m
       if iswetn(i)*iswetnt(i) == 1
           for k=1:kbm1
               if k == 1
                   temp=-wts(i,k+1)*(dye_age(i,k)*dz(i,k+1)...
                       +dye_age(i,k+1)*dz(i,k))/(dz(i,k)+dz(i,k+1));
               elseif k == kbm1
                       temp= wts(i,k)*(dye_age(i,k)*dz(i,k-1)...
                           +dye_age(i,k-1)*dz(i,k))/(dz(i,k)+dz(i,k-1));
               else
                   temp= wts(i,k)*(dye_age(i,k)*dz(i,k-1)+dye_age(i,k-1)*dz(i,k))...
                       /(dz(i,k)+dz(i,k-1))...
                       -wts(i,k+1)*(dye_age(i,k)*dz(i,k+1)+dye_age(i,k+1)*dz(i,k))...
                       /(dz(i,k)+dz(i,k+1));
               end

               if isonb(i) == 2
                   xflux(i,k)=temp*art1(i);    %/dz(k)
               else
                   xflux(i,k)=xflux(i,k)+temp*art1(i);    %/dz(k)
               end
           end
       end
    end

    %-------------------------------------------------------
    % Update Dye
    %------------------------------------------------------- 

    for i=1:m
     if iswetn(i)*iswetnt(i) == 1
         for k=1:kbm1
             % 最终更新标量值Ci的程序。
             % Final code for updating Ci.
             % dti : 为时间步长
             % dti : time step.
             % dt  : 上个时间水深
             % dt  : water depth at previous time step.
             % dtfa: 新的计算水深
             % dtfa: new water depth.
             dyef_age(i,k)=(dye_age(i,k)-xflux(i,k)/art1(i)*(dti/(dt(i)*dz(i,k))))...
                 *(dt(i)/dtfa(i))+c(i,k);
         end
     else
         for k=1:kbm1
             dyef_age(i,k)=dye_age(i,k)+c(i,k);
         end
     end
    end
    % END SUBROUTINE ADV_DYE
    
    %% SUBROUTINE VDIF_DYE(F)

    umolpr = umol*1.e0;

    %-------------------------------------------------------
    % The following section solves the equation
    % dti*(kh*f')'-f=-fb
    %-------------------------------------------------------

    for k = 2:kbm1
        for i = 1:m
            if iswetn(i) == 1
                fkh = kh(i,k);
                if k == kbm1
                    khbottom(i)=fkh;
                end
                af(i,k-1)=-dti*(fkh+umolpr)/(dz(i,k-1)*dzz(i,k-1)*d(i)*d(i));
                cf(i,k)=-dti*(fkh+umolpr)/(dz(i,k)*dzz(i,k-1)*d(i)*d(i));
            end
        end
    end

    %-------------------------------------------------------
    % The net heat flux input.
    % The method shown below can be used when we turn off the
    % body force in subroutine advt. Be sure this method could
    % cause the surface overheated if the vertical resolution
    % is not high enough.
    %-------------------------------------------------------

    swradf = zeros(m);
    wfsurf = zeros(m);
    rad    = zeros(m,kb);

    %-------------------------------------------------------
    % surface bcs; wfsurf
    %-------------------------------------------------------

    for i = 1:m
        if iswetn(i) == 1
            vhf(i,1) = af(i,1) / (af(i,1)-1.);
            vhpf(i,1) = -dti *(wfsurf(i)-swradf(i)+...
                rad(i,1)-rad(i,2)) / (-dz(i,1)*d(i)) - dyef_age(i,1);
            vhpf(i,1) = vhpf(i,1) / (af(i,1)-1.);
        end
    end

    for k = 2:kbm2
        for i = 1:m
            if iswetn(i) == 1
                vhpf(i,k)=1./ (af(i,k)+cf(i,k)*(1.-vhf(i,k-1))-1.);
                vhf(i,k) = af(i,k) * vhpf(i,k);
                vhpf(i,k) = (cf(i,k)*vhpf(i,k-1)-dyef_age(i,k)+...
                    dti*(rad(i,k)-rad(i,k+1))/(d(i)*dz(i,k)))*vhpf(i,k);
            end
        end
    end

    for k = 1:kbm1
        for i = 1:m
            if iswetn(i) == 1
                ff(i,k) = dyef_age(i,k);
            end
        end
    end

    for i = 1:m
        if iswetn(i) == 1 && isonb(i) ~= 2
            tmp1=pfpxb(i)*cos(sita_gd(i))+pfpyb(i)*sin(sita_gd(i));
            tmp2=ah_bottom(i)*phpn(i);
            tmp3=khbottom(i)+umolpr+ah_bottom(i)*phpn(i)*phpn(i);
            tmp=tmp1*tmp2/tmp3*(khbottom(i)+umolpr);

            tmp=0.0;
            gw=0.0;

            ff(i,kbm1) = ((cf(i,kbm1)*vhpf(i,kbm2)-ff(i,kbm1)-gw+...
                dti*(rad(i,kbm1)-rad(i,kb)-tmp)/(d(i)*dz(i,kbm1)))/...
                (cf(i,kbm1)*(1.-vhf(i,kbm2))-1.));
        end
    end

    for k = 2:kbm1
        ki = kb - k;
        for i = 1:m
            if iswetn(i) == 1 && isonb(i) ~= 2
                ff(i,ki) = (vhf(i,ki)*ff(i,ki+1)+vhpf(i,ki));
            end
        end
    end

    for i = 1:m
        if iswetn(i)*iswetnt(i) == 1
            for k = 1:kbm1
                dyef_age(i,k) = ff(i,k);
            end
        end
    end

    %------------------------------------------------------- 
    % Specify the source term
    %-------------------------------------------------------

    if inttime >= dyestart && inttime <= dyestop
        for kk=1:kspe_dye
            k=k_specify(kk);
            for j=1:mspe_dye
                i = nlid(m_specify(j));
                dyef_age(i,k)= dye_age_source_term;
            end
        end
    end

    % SUBROUTINE VDIF_DYE

    %% subroutine bcond_dye

    if iobcn > 0

    %------------------------------------------------------- 
    % set dye conditions on outer boundary
    %------------------------------------------------------- 

        for i=1:iobcn
            j=i_obc_n(i);
            j1=next_obc(i);
            s2d=0.0;
            s2d_next=0.0;
            xflux2d=0.0;
            for k=1:kbm1
                s2d=s2d+dye_age(j,k)*dz(j,k);
                s2d_next=s2d_next+dyef_age(j1,k)*dz(j1,k);
                xflux2d=xflux2d+xflux_obc(i,k);%*dz(k)
            end

            if uard_obcn(i) > 0.0 % if the flow is out of domain
                tmp=xflux2d+s2d*uard_obcn(i);
                s2d_obc=(s2d*dt(j)-tmp*dti/art1(j))/d(j);
                for k=1:kbm1
                    dyef_age(j,k)=dyef_age(j1,k);
                end
                for k=1:kbm1
                    smax = max(dye_age(nbsn(j,1:ntsn(j)),k));
                    smin = min(dye_age(nbsn(j,1:ntsn(j)),k));
                    if k == 1
                        smax = max(smax,(dye_age(j,k)*dz(j,k+1)+dye_age(j,k+1)*dz(j,k))/...
                            (dz(j,k)+dz(j,k+1)));
                        smin = min(smin,(dye_age(j,k)*dz(j,k+1)+dye_age(j,k+1)*dz(j,k))/...
                            (dz(j,k)+dz(j,k+1)));
                    elseif k == kbm1
                        smax = max(smax,(dye_age(j,k)*dz(j,k-1)+dye_age(j,k-1)*dz(j,k))/...
                            (dz(j,k)+dz(j,k-1)));
                        smin = min(smin,(dye_age(j,k)*dz(j,k-1)+dye_age(j,k-1)*dz(j,k))/...
                            (dz(j,k)+dz(j,k-1)));
                    else
                        smax = max(smax,(dye_age(j,k)*dz(j,k-1)+dye_age(j,k-1)*dz(j,k))/...
                            (dz(j,k)+dz(j,k-1)),...
                            (dye_age(j,k)*dz(j,k+1)+dye_age(j,k+1)*dz(j,k))/...
                            (dz(j,k)+dz(j,k+1)));
                        smin = min(smin,(dye_age(j,k)*dz(j,k-1)+dye_age(j,k-1)*dz(j,k))/...
                            (dz(j,k)+dz(j,k-1)),...
                            (dye_age(j,k)*dz(j,k+1)+dye_age(j,k+1)*dz(j,k))/...
                            (dz(j,k)+dz(j,k+1)));
                    end
                    if smin-dyef_age(j,k) > 0.0, dyef_age(j,k) = smin; end
                    if dyef_age(j,k)-smax > 0.0, dyef_age(j,k) = smax; end
                end
            else % if the flow is into the domain
                for k=1:kbm1
                    dyef_age(j,k)=dye_age(j,k);
                end
            end
        end
    end
end