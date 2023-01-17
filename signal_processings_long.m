close all
clear
clc

%% Import data
% Load strain data
sampleDir = 'Sample\Plants\';
expName = 'kp';
expDate = '2021-08-13';
fiberName = 'FS02020LUNA063028';
filename = [expName,'_',expDate,'_',fiberName];
load([sampleDir,filename]);
% Load FOS info
fiberDir = 'Sample\fiberInfo\';
load([fiberDir,fiberName,'_650um.mat']); % FOS index<->position
load([fiberDir,fiberName,'_650um_gages'])
S = min(gages);
E = max(gages);
pos = pos(S:E);
L = E-S+1;
n = 2; % number of measurements per trigger
triggerI = 360; % trigger interval in second

%% Process strain data
cutOff = 1e3; % > 1e-3 deformations are likely errors
nfTime = 24*3600; % notch filter cutoff time period in second
deltaT = nfTime/triggerI; % notch filter cutoff time period in index
dw = 1; % window size for time cutoff frequency
deltaL = 24; % positional threshold for low-pass filter, > 15 mm period
et = floor(length(allM(:,1))/n);
M = allM(1:et*n,:);
M(M.^2>cutOff^2) = nan;
% Time average
if (n>1)
    mM = reshape(M',[L,n,et]);
    mM = reshape(mean(mM,2),L,et)';
else
    mM = M;
end
mM = fillmissing(mM,'previous');
% Low-pass filter
kLP = L/deltaL; % positional cutoff frequency for low-pass filter
wNF = et/deltaT; % time cutoff frequency for notch filter
C = zeros(size(mM)); % filter
C(:,[1:ceil(kLP),end-ceil(kLP)+1:end]) = 1;
C([ceil(wNF)-dw:ceil(wNF)+dw,end-ceil(wNF)-dw+1:end-ceil(wNF)+dw+1],:) = 0;
mMhat = ifft2(fft2(mM).*C);
rmMhat = real(mMhat);
% Background subtraction
pM = bsxfun(@minus,rmMhat,mean(rmMhat,2));
pM = bsxfun(@minus,pM',mean(pM(1:deltaT,:))')';
pM = pM - mean(mean(pM));
clear kLP wNF C mM mMhat rmMhat
%% Plot heatmaps
% All
iTime = 1:et*n;
figure('Name','Heatmap of all','Position',[0 0 800 600]);
imagesc(iTime*triggerI/3600,pos,allM');
colorbar
axis tight; axis ij;
title(['Heatmap all for ',expDate])
xlabel('Time (h)'); ylabel('Position (m)')
set(gca,'FontSize',18);
clear iTime
% Processed
iTime = 1:et;
figure('Name','Heatmap after processing','Position',[0 0 800 600]);
imagesc(iTime*triggerI/3600,pos,pM');
colorbar
axis tight; axis ij;
title(['Heatmap processed for ',expDate])
xlabel('Position'); ylabel('Time points')
set(gca,'FontSize',18);
clear iTime