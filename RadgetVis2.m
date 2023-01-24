%Script for visualization of radish
close all
clear
clc

%% IO specification
sampleFolder = 'Sample\plants\';
xrayFolder = 'Sample\plants\';
expName = 'radish'; % plant species
expDate = '2022-01-13'; % experiment dates
fiberName = 'FS02021LUNA063307';
videoS = false;
videoTrimI = 3*3600/trimI; % interval video output

% Time domain
trimI = 10*60; % trimmed interval in second

% Load reformat strain data
load([sampleFolder,expName,expDate,fiberName,'_FOSreformat'])
% Load X-ray CT data
load([sampleFolder,expName,expDate,fiberName,'_CT'])

%% Plot heatmap of strain data
[et,L] = size(sMvalid);
iTime = 1:et; % time index
iTime = iTime*trimI/3600/24; % unit in day
figure('Name','2D heatmap of processed strain data','Position',[0 0 800 600])
imagesc(iTime,1:L,sMvalid');
colorbar;
axis tight; axis ij;
title('Heatmap of strain (\mu\epsilon)')
xlabel('Time (day)'); ylabel('FOS index')
set(gca,'fontsize',18)
clear iTime

%% 3D visualization
timeFrame = et;
sThresh = 50;
nColors = 101;
minM = 0; % minimum strain
maxM = 200; % maximum strain
nStrains = linspace(minM,maxM,nColors);
colorSet = parula(nColors);
XS = 40;
[EmeshX,EmeshY] = meshgrid(-XS:XS,-XS:XS);
colors = zeros([2*XS+1,2*XS+1,3]);
colors(:,:,1) = ones(2*XS+1).*linspace(0.5,1,2*XS+1);
colors(:,:,2) = 0;
colors(:,:,3) = 0;
if videoS
    tps = 1:videoTrimI:et;
    writerObj = VideoWriter([expName,fiberName,expDate,'_kinetics.avi']);
    writerObj.FrameRate = 4;
    open(writerObj);
else
    tps = timeFrame;
end
clear videoTrimI
f3d = figure('Name',['3D plot of t =',num2str(timeFrame)],'Position',[0 0 600 600]);
for i = 1:length(tps)
    plotSlice = sMvalid(tps(i),:)';
    rootFOS = zeros(2*XS+1,2*XS+1);
    if (min(plotSlice)<minM)
        minM = min(plotSlice);
        nStrains = linspace(minM,maxM,nColors); % update colormap interval
    end
    if (max(plotSlice)>maxM)
        maxM = max(plotSlice);
        nStrains = linspace(minM,maxM,nColors); % update colormap interval
    end
    for k = 1:nColors
        iplots = find(plotSlice>nStrains(k));
        plot3(rs(iplots).*cos(thetas(iplots)),rs(iplots).*sin(thetas(iplots)),...
                    zs(iplots),'.','markers',30,'Color',colorSet(k,:));
        hold on;
    end
    plotSlice(plotSlice<=sThresh) = nan;
    est_r = rs.*sqrt(plotSlice.^2+2*plotSlice)*1e-3; % radish size
    % FOS radish reconstruction
    for j = 1:L
        if ~isnan(est_r(j))
            xrange = ceil(est_r(j)*cos(thetas(j)))+XS+1:XS+1;
            if isempty(xrange)
                xrange = XS:ceil(est_r(j)*cos(thetas(j)))+XS;
            end
            yrange = ceil(est_r(j)*sin(thetas(j)))+XS+1:XS;
            if isempty(yrange)
                yrange = XS:ceil(est_r(j)*sin(thetas(j)))+XS;
            end
            rootFOS(xrange,yrange) = zs(j);
        end
    end
    if videoS
        text(0,20,10,sprintf('t = day %.1f',tps(i)/6/24),'FontSize',18)
    end
    s = surf(EmeshX,EmeshY,rootFOS,colors,'FaceAlpha',0.6);
    s.EdgeColor = 'red';
    axis equal;grid minor;
    axis([min(rs.*cos(thetas))-5 max(rs.*cos(thetas))+5 ...
        min(rs.*cos(thetas))-5 max(rs.*cos(thetas))+5 min(zs)-5 -0.1]);
    colorbar; caxis([minM,maxM]);
    view(60,8)
    set(gca,'fontsize',18)
    if videoS
        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
        hold off;
    end
    clear plotSlice iplots est_r xrange yrange frame j k
end
if videoS
    close(writerObj);
end
clear videoS tps nColors nStrains colorSet minM maxM s colors i

% %% X-ray CT reconstruction
colors = zeros([2*XS+1,2*XS+1,3])+0.2;
rootCT = zeros(2*XS+1,2*XS+1);
for i = 1:size(cM,1)
    x = cM(i,1);
    y = cM(i,2);
    z = cM(i,3);
    w = cM(i,4);
    h = cM(i,5);
    xrange = XS+round(x)+1:XS+round(x+w)+1;
    yrange = XS+round(y)+1:XS+round(y+h)+1;
    rootCT(xrange,yrange) = z;
end
hold on;
s = surf(EmeshX,EmeshY,rootCT,colors,'FaceAlpha',0.6);
s.EdgeColor = 'black';
clear c x y z w h xrange yrange colors XS s i
