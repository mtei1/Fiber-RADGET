close all
clear
clc

%% Import data
% Load strain data
sampleDir = 'Sample\Abiotic\';
expName = 'FGB30_12e-1';
expDate = '2022-07-27_1';
fiberName = 'FS02021LUNA010924';
filename = [expName,'_',expDate,'.csv'];
M = readmatrix([sampleDir,filename]);
% Load FOS information
fiberDir = 'Sample\fiberInfo\';
load([fiberDir,fiberName,'_650um.mat']); % FOS index<->position
hertz = 20; % Hz, measurement frequency
triggerI = 1/hertz; % sec, trigger interval
load([fiberDir,fiberName,'_650um_gages'])
S = min(gages);
E = max(gages);
L = E-S+1;
sensorRanges = gages-min(gages)+1;

%% Process strain data
et = floor(length(M(:,1))/hertz);
M = M(1:et*hertz,S:E)/10; % csv file is misisng decimal
mM = reshape(M',[L,hertz,et]);
mM = reshape(mean(mM,2),L,et)';
cutOff = 1000; % > 1e-3 deformations are likely errors
mM(mM.^2>cutOff^2) = nan;
mM = fillmissing(mM,'previous');
% Low-pass filter
deltaL = 24; % position threshold for low-pass filter, > 15 mm period
fLU = L/deltaL; % upper limit of positional frequency
C = zeros(size(mM)); % filter
C(:,:) = 1;
mMhat = ifft2(fft2(mM).*C);
rmMhat = real(mMhat);
% Jig bias subtraction
pM = bsxfun(@minus,rmMhat,rmMhat(1,:));
clear C mM mMhat rmMhat

%% Plot heatmap for overall
% Original
iTime = 1:et*hertz;
figure('Name','Heatmap of original','Position',[0 0 800 600]);
imagesc(iTime*triggerI,pos(S:E),M');
colorbar
axis tight; axis ij;
title('Heatmap original')
xlabel('Time (s)'); ylabel('Position (m)')
set(gca,'FontSize',18);
clear iTime
% Processed
iTime = 1:et;
figure('Name','Heatmap after signal processing','Position',[0 0 800 600]);
imagesc(iTime,pos(S:E),pM');
colorbar
axis tight; axis ij;
title('Heatmap processed')
xlabel('Time (s)'); ylabel('Position (m)')
set(gca,'FontSize',18);
clear iTime