clear all;
clc;

% ncfile = "test_0001.nc";
% nc.dye_source_term = 1;
% nc.umol            = 1.0E-05;                % Vertical mixing coefficient (1E-5)    
% nc.fact            = 1;
% nc.fm1             = 0;
% 
% nc.k_specify       = [1;2;3;4;5];            % NO. of sigma layer for specify dye release
% nc.m_specify       = [6];                    % NO. of node for specify dye release
% nc.dyestart        = mjuliandate(2015,01,01,00,00,00);
% nc.dyestop         = mjuliandate(2015,01,07,00,00,00);

% ncfile = "long_box_0001.nc";
% nc.dye_source_term = 1;
% nc.umol            = 1.0E-05;                % Vertical mixing coefficient (1E-5)    
% nc.fact            = 1;
% nc.fm1             = 0;
% 
% nc.k_specify       = [1;2];            % NO. of sigma layer for specify dye release
% nc.m_specify       = [988;987;986;1032;1034;1033;939;940];                   % NO. of node for specify dye release
% nc.dyestart        = mjuliandate(2015,01,01,00,00,00);
% nc.dyestop         = mjuliandate(2015,01,07,00,00,00);
% 

ncfile = "test_0001.nc";
nc.none = 0;
% nc.dye_source_term = 1;
% nc.umol            = 1.0E-05;                % Vertical mixing coefficient (1E-5)    
% nc.fact            = 1;
% nc.fm1             = 0;
% 
% nc.k_specify       = [1;2;3;4;5];            % NO. of sigma layer for specify dye release
% nc.m_specify       = [6];                    % NO. of node for specify dye release
% nc.dyestart        = mjuliandate(2015,01,01,00,00,00);
% nc.dyestop         = mjuliandate(2015,01,07,00,00,00);

nc                 = loadNetCDF(ncfile,nc);



% tic;
% bar = waitbar(0,'Calculation begin ...');
% contr = zeros(size(nc.dye));
% alpha = zeros(size(nc.dye));
% % inttime = [nc.timestart:nc.dti/3600/24:nc.timestop];
% inttime = [nc.timestart:nc.dti/(3600*24):nc.timestop];
% for i = 1:length(inttime)-1
%     nc.u = nc.uout(:,:,i);
%     nc.v = nc.vout(:,:,i);
%     dyef = dye(nc,contr(:,:,i),inttime(i));
%     contr(:,:,i+1)=dyef;
%     alpf = age(nc,alpha(:,:,i),inttime(i),contr(:,:,i));
%     alpha(:,:,i+1)=alpf;
%     
%     str=['Progress: ',num2str(100*i/(length(inttime)-1)),'%'];
%     waitbar(i/(length(inttime)-1),bar,str)
% end
% clear dyef alpf
% 
% wa = alpha ./ contr;
% wa(find(isnan(wa)==1)) = 0;
% 
% %close(bar)
% toc;

% plot(squeeze(nc.vout(6,1,2:end)))
% plot(squeeze(dye(6,1,2:end)))
% plot(squeeze(wa(6,1,2:end)))

% save('out_dye.mat','contr','-v7.3','-nocompression');
% save('out_age.mat','alpha','-v7.3','-nocompression');
% save('out_wa.mat','wa','-v7.3','-nocompression');


if ismac    % On Mac
    basedir = '/Users/yulong/GitHub/';
    addpath([basedir,'fvcomtoolbox/']);
    addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
    addpath([basedir,'fvcomtoolbox/utilities/']);
    addpath([basedir,'fvcomtoolbox/custom/']);
elseif isunix       % Unix?
    basedir = '/home/usr0/n70110d/github/';
    addpath([basedir,'fvcomtoolbox/']);
    addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
    addpath([basedir,'fvcomtoolbox/utilities/']);
    addpath([basedir,'fvcomtoolbox/custom/']);
elseif ispc     % Or Windows?
    basedir = 'C:/Users/Yulong WANG/Documents/GitHub/';      
    % basedir = 'C:/Users/Yulong/Documents/GitHub/';      
    addpath([basedir,'fvcomtoolbox/']);
    addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
    addpath([basedir,'fvcomtoolbox/utilities/']);
    addpath([basedir,'fvcomtoolbox/custom/']);
end
%%
% fig = figure(01);
% for tt = 1:length(nc.time)
%     plotMesh(01, [nc.lon, nc.lat], nc.tri, nc.dye_age(:,1,tt)./nc.dye(:,1,tt),...
%         'time',nc.Times(tt,1:19),...
%         'dye age',{'Time (s)',[0,10000]});
%     drawnow
%     frame = getframe(fig);
%     im{tt} = frame2im(frame);
% end
% 
%     plotMesh(01, [nc.lon, nc.lat], nc.tri, nc.dye(:,1,tt),...
%         'time',nc.Times(tt,1:19),...
%         'dye age',{'Time (s)',[0,0.1]});
%     
%     plotMesh(01, [nc.lon, nc.lat], nc.tri, nc.dye_age(:,1,tt),...
%         'time',nc.Times(tt,1:19),...
%         'dye age',{'Time (s)',[0,1000]});
% 
% fig = figure(01);
% for tt = 1:length(nc.time)
%     plotMesh(01, [nc.lon, nc.lat], nc.tri, nc.dye_age(:,1,tt),...
%         'time',nc.Times(tt,1:19),...
%         'dye age',{'Time (s)',[0,500]});
%     drawnow
%     frame = getframe(fig);
%     im{tt} = frame2im(frame);
% end
% filename = '01_dye.gif'; 
% for tt = 1:length(nc.time)
%     [A,map] = rgb2ind(im{tt},256);
%     if tt == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
%     end
% end

fig = figure(01);
for tt = 1:length(nc.time)
    plotMesh(01, [nc.lon, nc.lat], nc.tri, nc.dye(:,1,tt),...
        'time',nc.Times(tt,1:19),...
        'dye',{'Concentration (-)',[0,1]});
    drawnow
    frame = getframe(fig);
    im{tt} = frame2im(frame);
end
filename = '101_dye.gif'; 
for tt = 1:length(nc.time)
    [A,map] = rgb2ind(im{tt},256);
    if tt == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end

fig = figure(02);
for tt = 1:73%length(nc.time)
    plotMesh(02, [nc.lon, nc.lat], nc.tri, nc.dye_age(:,1,tt),...
        'time',nc.Times(tt,1:19),...
        'dye age',{'Time (s)',[0,5*length(nc.time)]});
    drawnow
    frame = getframe(fig);
    im{tt} = frame2im(frame);
end
filename = '02_dye_age.gif'; 
for tt = 1:73%length(nc.time)
    [A,map] = rgb2ind(im{tt},256);
    if tt == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end

water_mask = squeeze(nc.dye(:,1,:)<=1E-6);
water_age = squeeze(nc.dye_age(:,1,:)./nc.dye(:,1,:));
water_age(water_mask)=nan;

% nc.dye_age(:,1,tt)./nc.dye(:,1,tt)

fig = figure(03);
for tt = 1:length(nc.time)
    plotMesh(03, [nc.lon, nc.lat], nc.tri, water_age(:,tt)/24/60,...
        'time',nc.Times(tt,1:19),...
        'dye age',{'Time (day)',[0,3]});
    drawnow
    frame = getframe(fig);
    im{tt} = frame2im(frame);
end
filename = '101_03_water_age.gif'; 
for tt = 1:length(nc.time)
    [A,map] = rgb2ind(im{tt},256);
    if tt == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end

plotWaterAge('test_0001.nc',1);
%%
list_node = [609,605,597,591,585,579,573,567,561,555,549,543,537,531,525,519,...
    513,507,501,495,489,483,477,471,465,459,453,447,441,435,429,423,417,...
    411,405,399,393,387,381,375,369,363,357,351,345,339,333,327,321,315,...
    309,303,297,291,285,279,273,267,261,255,249,243,237,231,225,219,213,...
    207,201,195,189,183,177,171,165,159,153,147,141,135,129,123,117,111,...
    105,99,93,87,81,75,69,63,57,51,45,39,33,27,21,14,7,2];

nc.ua       = ncread(ncfile,'ua');
nc.va       = ncread(ncfile,'va');

for i = 1:length(list_node)
    j = list_node(i);
    list_lat(i) = nc.lat(j);
    list_dye(i,:,:) = nc.dye(j,:,:);
    list_dye_age(i,:,:) = nc.dye_age(j,:,:);
    list_ua(i,:) = nc.ua(j,:);
    list_va(i,:) = nc.va(j,:);
end

for i = 1:length(list_node)
    list_x(i) = (max(list_lat)-list_lat(i))/(max(list_lat)-min(list_lat));
end

list_water_mask = list_dye<=1E-9;
list_dye(list_water_mask)=nan;
list_water_age = list_dye_age./list_dye;
list_water_age(list_water_mask)=nan;
list_water_age(list_water_age<0)=nan;
list_water_age_day = list_water_age/24;
list_water_age_day_surface = squeeze(list_water_age_day(:,1,:));

% list_water_age_day_surface(end,end)
fig = figure(01);
for tt = 1:length(nc.time)
    yyaxis right;
    scatter(list_x,list_water_age_day_surface(:,tt));
%     plot(list_x,list_water_age_day_surface(:,tt),'-o');
%     axis([0 1 0 200]);% 001 case
%     axis([0 1 0 400]);% 001 case
%     axis([0 1 0 2000]);% 001 case
    axis([0 1 0 200]);% 001 case
    ylabel('Water age (day)','FontSize',14);
    yyaxis left;
    scatter(list_x,list_dye(:,1,tt));
%     plot(list_x,list_dye(:,1,tt),'-o');
    axis([0 1 0 1.2]);
    ylabel('Concentration (-)','FontSize',14);
    xlabel('x (-)','FontSize',14);
    title(['C and WA (',nc.Times(tt,1:19),')'],'FontSize',14);
    speed = -1 * mean(list_va(:,tt));
    speed = num2str(speed, '%.6f');
    txt = 'mean velocity';
    text(0.05,1.1,[txt,' = ',speed,' m/s'],'FontSize',14);
    drawnow
    frame = getframe(fig);
    im{tt} = frame2im(frame);
end
filename = '04_concentration_water_age.gif'; 
for tt = 1:73%length(nc.time)
    [A,map] = rgb2ind(im{tt},256);
    if tt == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end

fig = figure(02);
for tt = 1:length(nc.time)
    list_spd(tt) = -1 * mean(list_va(:,tt));
end
scatter([1:length(nc.time)],list_spd);
axis([0 length(nc.time) 0.00 0.1]);
ylabel('Mean Velocity (m/s)','FontSize',14);
xlabel('time (hour)','FontSize',14);
title(['Mean velocity (m/s)'],'FontSize',14);

% check how long time to flow cannal
ans = list_dye(:,1,:);
ans = squeeze(ans);
plot(ans(end,:));
axis([0 75 0 1]);
% axis([0 150 0 1]);
% axis([0 750 0 1]);

% 1000m3/s * 4 / 910.3730m / 50m = 0.0879m/s
% 70*3600*0.0879=22185m(*0.75=16639)
% at about 66 hour flow go throgh canal.
% 0500m3/s * 4 / 910.3730m / 50m = 0.0440m/s
% 140*3600*0.0437=22024m(*0.75=16639)
% at about 131 hour flow go throgh canal.
% 0100m3/s * 4 / 910.3730m / 50m = 0.0088m/s
% 701*3600*0.0088=22024m(*0.75=16639)
% at about 701 hour flow go throgh canal.

% 22193.75499782758m