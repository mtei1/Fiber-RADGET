%%Script for signal processing, data formatting, visualization for radish
close all
clear
clc

%% IO specification
sampleFolder = 'Sample\Plants\';
configFolder = 'Sample\Configs\';
expName = 'radish'; % plant species
fiberName = 'FS02021LUNA063307';
expDate = '2022-01-13'; % experiment dates
trimI = 10*60; % trimmed interval in second
gageI = 0.65; % FOS interval in mm
Tc = 1.7e-4; % thermal expansion coefficient for PTFE
Kt = 0.6800000071525574;
Ke = 6.716787338256836;
videoS = true;

%% Fiber-RADGET dimensions
load([configFolder,fiberName,'_650um_gages'])
load([configFolder,'Fiber-RADGET_dim_',expName])
gages = gages(2:end);
S = min(gages);
E = max(gages);
sensorRanges = gages-S+1;

%% Import full experiment data
%Load strain data
d = dir([sampleFolder,expName,expDate,fiberName,'*.mat']);
ets = zeros(size(d));
for i = 1:length(d)
    load([sampleFolder,d(i).name]);
    if i == 1
        sM = mMall(:,S:E);
        [et,L] = size(sM);
    else
        tmp = sM;
        et = size(mMall,1);
        sM = nan(sum(ets)+et,L);
        sM(1:sum(ets),:) = tmp;
        sM(sum(ets)+1:end,:) = mMall(:,S:E);
    end
    ets(i) = et;
end
%Load CT data
load([sampleFolder,'CT_',expName,expDate,fiberName])

%% Process strain data
T = sum(ets);
%Band-pass filter dimensions
deltaL = ceil(15/gageI); % positional threshold for low-pass filter in index
deltaT = ceil(24*3600/trimI); % notch filter cutoff time period in index
dw = 5; % window size for time cutoff frequency
%Band-pass filter
kLP = L/deltaL; % positional cutoff frequency for low-pass filter
lNF = T/deltaT-dw; % time cutoff frequency for notch filter
uNF = T/deltaT+dw; % time cutoff frequency for notch filter
C = zeros(size(sM)); % filter
C(:,[1:ceil(kLP),end-ceil(kLP)+1:end]) = 1;
for j = 1:8
    C([ceil(j*lNF):ceil(j*uNF),end-ceil(j*lNF)+1:end-ceil(j*uNF)+1],:) = 0;
end
clear lNF uNF dw kLP
sMhat = ifft2(fft2(sM).*C);

%% Background subtraction
sMhat = real(sMhat);
sMhat = bsxfun(@minus,sMhat,mean(sMhat,2)); % Subtract time effect
sMhat = bsxfun(@minus,sMhat',mean(sMhat(1:deltaT,:))')'; % Subtract jig effect
sMhat = sMhat - mean(mean(sMhat));

%% Reformat strain data
Ls = abs(sensorRanges(1:5)-sensorRanges(6:10))+1;
sMvalid = nan(sum(ets),sum(Ls));
thetas = nan(sum(Ls),1);
zs = thetas;
for i = 1:length(Ls)
    tr = sensorRanges(i):sensorRanges(i+5);
    currTh = MAXFTH;
    swirl = -1;
    if isempty(tr) % FOS is wrapped in theta direction 0 -> pi
        tr = sensorRanges(i+5):sensorRanges(i);
        currTh = MINFTH;
        swirl = 1;
    end
    if i > 1
        cumL = sum(Ls(1:i-1));
    else
        cumL = 0;
    end
    for j = 1:Ls(i)
        sMvalid(:,cumL+j) = sMhat(:,tr(j));
        thetas(cumL+j) = currTh;
        zs(cumL+j) = sensorDepths(i);
        if swirl < 0
            currTh = currTh - (a*currTh+MINR-sqrt((a*currTh+MINR)^2-2*a*gageI))/a;
        else
            currTh = currTh + (sqrt((a*currTh+MINR)^2+2*a*gageI)-(a*currTh+MINR))/a;
        end
    end
    clear tr currTh swirl cumL
end
rs = a*(thetas-MINTH)+MINR;

%% Plot heatmap of strain data
iTime = 1:sum(ets); % time index
iTime = iTime*trimI/3600/24; % unit in day
%Original plot
figure('Name','2D heatmap of original strain data','Position',[0 0 800 600])
imagesc(iTime,1:size(mMall,2),mMall');
colorbar;
axis tight; axis ij;
title('Heatmap of strain (\mu\epsilon)')
xlabel('Time (day)'); ylabel('FOS index')
set(gca,'fontsize',18)
%Processed plot
figure('Name','2D heatmap of processed strain data','Position',[0 0 800 600])
imagesc(iTime,1:sum(Ls),sMvalid');
colorbar;
axis tight; axis ij;
title('Heatmap of strain (\mu\epsilon)')
xlabel('Time (day)'); ylabel('FOS index')
set(gca,'fontsize',18)
clear iTime

%% 3D visualization
videoTrimI = 3*3600/trimI; % video output interval
finalFrame = T;
sThresh = 50; % microstrain threshold for radish visualization
nColors = 101; % color discreteness
minStrain = 0; % minimum strain
maxStrain = 200; % maximum strain
strain2color = linspace(minStrain,maxStrain,nColors);
colorSet = parula(nColors);
XS = 40; % virtual radish mesh size
[EmeshX,EmeshY] = meshgrid(-XS:XS,-XS:XS); % backbone mesh
colors = zeros([2*XS+1,2*XS+1,3]); % fake color for Fiber-RADGET reconstruction
colors(:,:,1) = ones(2*XS+1).*linspace(0.5,1,2*XS+1);
colors(:,:,2) = 0;
colors(:,:,3) = 0;
if videoS
    tps = 1:videoTrimI:T;
    writerObj = VideoWriter([expName,fiberName,expDate,'_kinetics.avi']);
    writerObj.FrameRate = 4;
    open(writerObj);
else
    tps = finalFrame;
end
f3d = figure('Name',['3D plot of t =',num2str(finalFrame)],'Position',[0 0 600 600]);
for i = 1:length(tps)
    currStrain = sMvalid(tps(i),:)';
    rootFOS = zeros(2*XS+1,2*XS+1); % reconstruction from FOS
    if (min(currStrain)<minStrain)
        minStrain = min(currStrain);
        strain2color = linspace(minStrain,maxStrain,nColors); % update colormap interval
    end
    if (max(currStrain)>maxStrain)
        maxStrain = max(currStrain);
        strain2color = linspace(minStrain,maxStrain,nColors); % update colormap interval
    end
    for k = 1:nColors
        iplots = find(currStrain>strain2color(k));
        plot3(rs(iplots).*cos(thetas(iplots)),rs(iplots).*sin(thetas(iplots)),...
                    zs(iplots),'.','markers',30,'Color',colorSet(k,:));
        hold on;
    end
    currStrain(currStrain<=sThresh) = nan;
    est_r = rs.*sqrt(currStrain.^2+2*currStrain)*1e-3; % radish size
    for j = 1:sum(Ls)
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
    text(0,42,10,sprintf('t = day %.1f',tps(i)/6/24),'FontSize',18)
    s = surf(EmeshX,EmeshY,rootFOS,colors,'FaceAlpha',0.6);
    s.EdgeColor = 'red';
    axis equal;grid minor;
    axis([-MAXR MAXR -MAXR MAXR MIND -1]);
    colorbar; clim([minStrain,maxStrain]);
    view(60,8)
    set(gca,'fontsize',18)
    if videoS
        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
        hold off;
    end
    clear currStrain iplots est_r xrange yrange frame
end
if videoS
    close(writerObj);
end
clear tps nColors strain2color colorSet minStrain maxStrain s colors

%% Plot X-ray CT result
colors = zeros([2*XS+1,2*XS+1,3])+0.2; % fake color for CT reconstruction
rootCT = zeros(2*XS+1,2*XS+1); % reconstruction from X-ray CT
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
clear x y z w h xrange yrange colors
