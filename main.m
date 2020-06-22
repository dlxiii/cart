clear all;
clc;

% ncfile = "test_0001.nc";
% nc.dye_source_term = 1;
% nc.umol            = 1.0E-05;                % Vertical mixing coefficient (1E-5)    
% nc.fact            = 1;
% nc.fm1             = 0;
% 
% nc.k_specify       = [1;2;3;4;5];            % NO. of sigma layer for specify dye release
% nc.m_specify       = [6];                   % NO. of node for specify dye release
% nc.dyestart        = mjuliandate(2015,01,01,00,00,00);
% nc.dyestop         = mjuliandate(2015,01,07,00,00,00);

ncfile = "long_box_0001.nc";
nc.dye_source_term = 1;
nc.umol            = 1.0E-05;                % Vertical mixing coefficient (1E-5)    
nc.fact            = 1;
nc.fm1             = 0;

nc.k_specify       = [1;2];            % NO. of sigma layer for specify dye release
nc.m_specify       = [988;987;986;1032;1034;1033;939;940];                   % NO. of node for specify dye release
nc.dyestart        = mjuliandate(2015,01,01,00,00,00);
nc.dyestop         = mjuliandate(2015,01,07,00,00,00);

nc                 = loadNetCDF(ncfile,nc);

%%

tic;
bar = waitbar(0,'Calculation begin ...');
contr = zeros(size(nc.dye));
alpha = zeros(size(nc.dye));
% inttime = [nc.timestart:nc.dti/3600/24:nc.timestop+1];
inttime = [nc.timestart:3600/3600/24:nc.timestop+1];
for i = 1:length(inttime)-1
    nc.u = nc.uout(:,:,i);
    nc.v = nc.vout(:,:,i);
    dyef = dye(nc,contr(:,:,i),inttime(i));
    contr(:,:,i+1)=dyef;
    alpf = age(nc,alpha(:,:,i),inttime(i),contr(:,:,i));
    alpha(:,:,i+1)=alpf;
    
    str=['Progress: ',num2str(100*i/(length(inttime)-1)),'%'];
    waitbar(i/(length(inttime)-1),bar,str)
end
clear dyef alpf

wa = alpha ./ contr;
wa(find(isnan(wa)==1)) = 0;

%close(bar)
toc;

% plot(squeeze(nc.vout(6,1,2:end)))
% plot(squeeze(dye(6,1,2:end)))
% plot(squeeze(wa(6,1,2:end)))

% save('out_dye.mat','contr','-v7.3','-nocompression');
% save('out_age.mat','alpha','-v7.3','-nocompression');
% save('out_wa.mat','wa','-v7.3','-nocompression');
%%

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
% 
fig = figure(01);
for tt = 1:length(nc.time)
    plotMesh(01, [nc.lon, nc.lat], nc.tri, contr(:,1,tt),...
        'time',nc.Times(tt,1:19),...
        'dye',{'Concentration (-)',[0,1]});
    drawnow
    frame = getframe(fig);
    im{tt} = frame2im(frame);
end
% filename = '01_dye.gif'; 
% for tt = 1:length(nc.time)
%     [A,map] = rgb2ind(im{tt},256);
%     if tt == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
%     end
% end
