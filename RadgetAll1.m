%Script for signal processing, data formatting, visualization for rice
close all
clear
clc

%% IO specification
sampleFolder = 'Sample\Plants\';
configFolder = 'Sample\Configs\';
expName = 'kp'; % plant species
fiberName = 'FS02020LUNA063328';
expDate = '2021-07-*'; % experiment dates
trimI = 10*60; % trimmed interval in second
videoS = false;

%% Fiber-RADGET dimensions
load([configFolder,fiberName,'_650um_gages'])
load([configFolder,'Fiber-RADGET_dim_',expName])
Tc = 1.7e-4; % thermal expansion coefficient for PTFE
Kt = 0.6800000071525574;
Ke = 6.716787338256836;
gageI = 0.65; % FOS interval in mm
S = min(gages);
E = max(gages);

%% Simulation parameters
P = 5e5;
nu = 0.5; % Poisson ratio 0-0.5
G = .1; % GPa

%% Import full experiment data
d = dir([sampleFolder,expName,expDate,fiberName,'*.mat']);
ets = zeros(size(d));
cMall = cell(size(d)); % CT data storage
for i = 1:length(d)
    load([sampleFolder,d(i).name]);
    if i == 1
        mMall = mM(:,S:E);
        [et,L] = size(mMall);
    else
        et = size(mM,1);
        tmp = mMall;
        mMall = nan(sum(ets)+et,L);
        mMall(1:sum(ets),:) = tmp;
        mMall(sum(ets)+1:end,:) = mM(:,S:E);
        clear tmp
    end
    ets(i) = et;
    % Load x-ray CT data
    try
        load([configFolder,'CT_dim_',d(i).name]);
    catch
    end
    try
        load([sampleFolder,'CT_',d(i).name]);
        cM = [bsxfun(@minus,rinfo(:,2),x0)*pix2mm,... % x-coordinate
              bsxfun(@plus,-rinfo(:,3),y0)*pix2mm,... % y-coordinate
              bsxfun(@minus,rinfo(:,1),z0)*z2mm]; % z-coordinate
        cMall{i} = cM;
    catch
    end
end
clear et mM rinfo cM

%% Signal processing
%Band-pass filter sizes
deltaL = ceil(15/gageI); % positional threshold for low-pass filter in index
kLP = L/deltaL; % positional cutoff frequency for low-pass filter
deltaT = ceil(24*3600/trimI); % notch filter cutoff time period in index
dw = 2; % window size for time cutoff frequency
fMall = mMall;
for i = 1:length(ets)
    T = ets(i);
    if i ==1
        currRange = 1:ets(1);
    else
        currRange = sum(ets(1:i-1))+1:sum(ets(1:i));
    end
    sM = mMall(currRange,:);
    %Band-pass filter
    C = zeros(size(sM)); % filter
    C(:,[1:ceil(kLP),end-ceil(kLP)+1:end]) = 1; % low-pass
    lNF = T/deltaT-dw; % time cutoff frequency for notch filter
    uNF = T/deltaT+dw; % time cutoff frequency for notch filter
    for j = 1:24 % notch-filter
        C([ceil(j*lNF):ceil(j*uNF),end-ceil(j*lNF)+1:end-ceil(j*uNF)+1],:) = 0;
    end
    sMhat = real(ifft2(fft2(sM).*C));
    % Background subtraction
    sMhat = bsxfun(@minus,sMhat',mean(sMhat(1:deltaT,:))')'; % Subtract jig effect
    sMhat = bsxfun(@minus,sMhat,mean(sMhat,2)); % Subtract time effect
    sMhat = sMhat - mean(mean(sMhat));
    fMall(currRange,:) = sMhat;
    clear lNF uNF C T currRange sMhat
end

%% Reformat strain data
Ls = abs(gages(2:8)-gages(9:15))+1;
sMvalid = nan(sum(ets),sum(Ls));
thetas = nan(sum(Ls),1);
zs = thetas;
for i = 1:length(Ls)
    tr = gages(i+1):gages(i+8);
    currTh = MINFTH;
    swirl = 1;
    if isempty(tr) % FOS is wrapped in reverse pi -> 0
        tr = gages(i+8):gages(i+1);
        currTh = MAXFTH;
        swirl = -1;
    end
    if i ==1
        cumL = 0;
    else
        cumL = sum(Ls(1:i-1));
    end
    for j = 1:length(tr)
        sMvalid(:,cumL+j) = fMall(:,tr(j));
        thetas(cumL+j) = currTh;
        zs(cumL+j) = sensorDepths(i);
        if swirl<0
            currTh = currTh - (a*currTh+MINR-sqrt((a*currTh+MINR)^2-2*a*gageI))/a; % approximate angle
        else
            currTh = currTh + (sqrt((a*currTh+MINR)^2+2*a*gageI)-(a*currTh+MINR))/a; % approximate angle
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
imagesc(iTime,1:L,mMall');
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
iExp = 2;
finalFrame = sum(ets(1:iExp));
nColors = 101; % color discreteness
minStrain = 0; % minimum strain
maxStrain = 200; % maximum strain
strain2color = linspace(minStrain,maxStrain,nColors);
colorSet = parula(nColors);
if videoS
    tps = 1:videoTrimI:finalFrame;
    writerObj = VideoWriter([expName,fiberName,expDate,'_kinetics.avi']);
    writerObj.FrameRate = 4;
    open(writerObj);
else
    tps = finalFrame;
end
f3d = figure('Name','3D plot of FOS kinetics','Position',[0 0 600 600]);
for i = 1:length(tps)
    currStrain = sMvalid(tps(i),:);
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
            zs(iplots),'.','markers',30,'Color',colorSet(k,:))
        hold on;
    end
    text(0,42,10,sprintf('t = day %.1f',tps(i)*trimI/3600/24),'FontSize',18)
    axis equal; grid minor;
    axis([-MAXR*1.25 MAXR*1.25 -MAXR*1.25 MAXR*1.25 MIND 0]);
    colorbar; clim([minStrain,maxStrain])
    view(-10,10);
    set(gca,'fontsize',18);
    clear currStrain iplots frame k
end
clear tps nColors strain2color colorSet minStrain maxStrain

% Initial solution
sThresh = 20; % strain threshold for root detection
numRoots = 30; % maximum number of roots per layer
iPeaks = nan(length(Ls),numRoots);
for i = 1:length(Ls)
    if i ==1
        cumL = 0;
    else
        cumL = sum(Ls(1:i-1));
    end
    tr = sMvalid(finalFrame,cumL+1:cumL+Ls(i));
    grad = tr(2:end) - tr(1:end-1);
    ip = find((sign(grad(1:end-1))-sign(grad(2:end)))==2); % de/dx=0, d^2e/(dx)^2 <0 
    ip = ip(tr(ip)>sThresh); % peak strain value > threshold
    for j = 1:length(ip)
        iPeaks(i,j) = cumL+ip(j);
    end
    clear tr grad ip cumL
end
iPeaks = iPeaks(~isnan(iPeaks));
numRoots = numel(iPeaks);
x0 = nan(3*numRoots,1);
P0 = MINR^2*sThresh; % strain at the pressure origin
for i = 1:numRoots
    r_est = rs(iPeaks(i))-sqrt(P0./sMvalid(finalFrame,iPeaks(i))); % estimated r of root
    x0(3*i-2:3*i) = [r_est.*cos(thetas(iPeaks(i)));r_est.*sin(thetas(iPeaks(i)));...
        zs(iPeaks(i))];
end
clear r_est
% Kelvin simulations
obs = sMvalid(finalFrame,:);
obs(obs<sThresh) = 0;
% Define error function
pxs = repmat(rs.*cos(thetas),1,numRoots);
pys = repmat(rs.*sin(thetas),1,numRoots);
pzs = repmat(zs,1,numRoots);
len = length(obs);
Kelvinn = @(x) sum(P/16/pi./...
    ( (pxs-repmat(x(1:3:end)',len,1)).^2 + ...
    (pys-repmat(x(2:3:end)',len,1)).^2 + ...
    (pzs-repmat(x(3:3:end)',len,1)).^2 ).^(3/2).*...
    abs(repmat(x(3:3:end)',len,1) - pzs ) ,2 );
Boussy = @(x) sum(P/4/pi/G*...
    ((2*nu-1)./ (sqrt( (repmat(x(3:3:end),len,1)-pzf).^2 ).*sqrt...
    (repmat(x(1:3:end),len,1).^2+prf.^2-...
    2*repmat(x(1:3:end),len,1).*prf.*cos(repmat(x(2:3:end),len,1)-ptf)...
    +(repmat(x(3:3:end),len,1)-pzf).^2) +...
    repmat(x(1:3:end),len,1).^2+prf.^2 -...
    2*repmat(x(1:3:end),len,1).*prf.*cos(repmat(x(2:3:end),len,1)-ptf)...
    +(repmat(x(3:3:end),len,1)-pzf).^2) + sqrt( (repmat(x(3:3:end),len,1)-pzf).^2 )./(...
    (repmat(x(1:3:end),len,1).^2+prf.^2-...
    2*repmat(x(1:3:end),len,1).*prf.*cos(repmat(x(2:3:end),len,1)-ptf)...
    +(repmat(x(3:3:end),len,1)-pzf).^2).^(3/2)) ),2);
% Optimize solution
errorFun = @(x) sum((obs'-Kelvinn(x)).^2);
lb = -MINR*ones(3*numRoots,1);
ub = MINR*ones(3*numRoots,1);
lb(3:3:end) = sensorDepths(end-1);
ub(3:3:end) = 0;
A = eye(3*numRoots);
b = [MINR;MINR;0];
options = optimoptions(@fmincon,'Algorithm','sqp','MaxFunctionEvaluations',1e5);
estE = Kelvinn(x0);
errorR = ( sum((obs'-estE).^2) / sum(obs.^2) ).^(1/2);
sprintf('After optimization normalized sum of squared error: %f', errorR)
figure('Name','Fit comparison','Position',[0.1 0.1 800 600])
plot(obs,'b-','linewidth',1.5);
hold on
plot(estE,'k-','linewidth',1.5);
x = fmincon(errorFun,x0,A,repmat(b,numRoots,1),[],[],lb,ub,[],options);
estE = Kelvinn(x);
errorR = ( sum((obs'-estE).^2) / sum(obs.^2) ).^(1/2);
sprintf('After optimization normalized sum of squared error: %f', errorR)
x_est = x(1:3:end);
y_est = x(2:3:end);
z_est = x(3:3:end);
plot(estE,'r-','linewidth',1.5)
legend('Observed','Initial fit','Optimized fit')
grid minor
ylabel('Strain (\mu\epsilon)')
xlabel('Sensor position index')
set(gca,'FontSize',18)
set(gcf,'PaperPositionMode','auto')
clear len A b lb ub options pxs pys pzs x estE
% Plot solution
r_est = sqrt(x_est.^2+y_est.^2);
theta_est = acos(x_est./r_est);
theta_est(y_est<0) = -theta_est(y_est<0);
figure(f3d)
hold on;
plot3(x_est,y_est,z_est,'ro','linewidth',1.5)
for i = 1:numRoots
    sqrr = linspace(0,sqrt(r_est(i)),10);
    zz = linspace(0,z_est(i),10);
    plot3(sqrr.^2.*cos(theta_est(i)),sqrr.^2.*sin(theta_est(i)),zz,'r-','linewidth',1.5)
end

%% Overlay CT data
cM = cMall{iExp};
idx0 = find(cM(:,end)==0); % surface index
figure(f3d);
hold on;
plot3(cM(1:idx0(1),1),cM(1:idx0(1),2),cM(1:idx0(1),3),...
  '-','color',[0.5 0.5 0.5],'LineWidth',2)
for i = 2:length(idx0)
  plot3(cM(idx0(i-1)+1:idx0(i),1),cM(idx0(i-1)+1:idx0(i),2),...
      cM(idx0(i-1)+1:idx0(i),3),'-','color',[0.5 0.5 0.5],'LineWidth',2)
end

%% Comparison
nbins = 12;
NE = nan(nbins,1);
NC = NE;
for i = 1:nbins
    NE(i) = sum(z_est<MIND/nbins*(i-1));
    NC(i) = sum(cM([1;idx0(1:end-1)+1],3)<MIND/nbins*(i-1));
end
figure('Name','Reconstruction methods comparison','Position',[0 0 600 400])
plot(-MIND/nbins*linspace(0,nbins-1,12),NC,'-s','linewidth',2)
hold on
plot(-MIND/nbins*linspace(0,nbins-1,12),NE,'-o','linewidth',2)
axis([0 120 0 20])
legend('From CT','From FOS')
grid minor
ylabel('Number of roots')
xlabel('Depth (mm)')
set(gca,'FontSize',18)