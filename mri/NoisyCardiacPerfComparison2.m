%%
% Multicoil cardiac perfusion MRI reconstruction algorithms in the presence
% of noise
%

rng(2);

%--------------------------------------------------------------------------
% Knobs
%--------------------------------------------------------------------------
NSigma = 13;
sigmalist = linspace(1,25,13);

NParams = 5;
Plist = logspace(log10(0.1),log10(10),NParams);

Pbest_ist =  [3 3 3 4 4 4 4 5 5 5 5 5 5];
Pbest_slrn = [1 1 2 2 2 2 3 3 3 3 3 3 3];

lambdaList = 0.01;
lambdaSist = 0.01;

rslrn = 1;
lambdaSslrn = 0.1;

tol = 0.0025;
MaxIters = 40;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Load data
%--------------------------------------------------------------------------
% Data
load('cardiac_perf_R8.mat');
[nx ny nt nc] = size(kdata);
nd = nx * ny;
E = Emat_xyt(kdata(:,:,:,1) ~= 0,b1); % Aquisition operator
T = TempFFT(3); % Sparsifying transformation : 3D temporal DFT
%d = kdata; % Sub-sampled k-space data
%idx = (d ~= 0);

% True images
load('./figures/LS_cardiac_ist.mat');
Mtrue1 = abs(Ltrue + Strue);
load('./figures/LS_cardiac_SLRN.mat');
Mtrue2 = abs(Ltrue + Strue);
Mtrue = (Mtrue1 + Mtrue2) / 2;

% True sub-sampled k-space data
d = E * Mtrue;
idx = (d ~= 0);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Data storage
%--------------------------------------------------------------------------
ISTn.sigma = nan(NSigma,1);
ISTn.P = nan(NSigma,1);
ISTn.L = cell(NSigma,1);
ISTn.S = cell(NSigma,1);
ISTn.cost = cell(NSigma,1);
ISTn.deltaM = cell(NSigma,1);
ISTn.time = cell(NSigma,1);
ISTn.its = nan(NSigma,1);
ISTn.NRMSE = cell(NSigma,1);
ISTn.NRMSEend = nan(NSigma,1);
ISTn.SNR = cell(NSigma,1);
ISTn.SNRend = nan(NSigma,1);

SLRNn.sigma = nan(NSigma,1);
SLRNn.P = nan(NSigma,1);
SLRNn.L = cell(NSigma,1);
SLRNn.S = cell(NSigma,1);
SLRNn.cost = cell(NSigma,1);
SLRNn.deltaM = cell(NSigma,1);
SLRNn.time = cell(NSigma,1);
SLRNn.its = nan(NSigma,1);
SLRNn.NRMSE = cell(NSigma,1);
SLRNn.NRMSEend = nan(NSigma,1);
SLRNn.SNR = cell(NSigma,1);
SLRNn.SNRend = nan(NSigma,1);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Perform simulations
%--------------------------------------------------------------------------
for i = 1:NSigma
    % Add noise
    sigma = sigmalist(i);
    noise = (sigma / sqrt(2 * nd)) * (randn(size(d)) + 1i * randn(size(d)));
    dn = zeros(size(d));
    dn(idx) = d(idx) + noise(idx);
    C = max(svd(reshape(E' * dn,[nd nt]))); % scaling constant for Otazo's SVT
    
    % Display progress
    fprintf('\nTrial %i/%i\n',i,NSigma);
    
    % Otazo's L + S algo
    P1 = Plist(Pbest_ist(i));
    [L S cost deltaM time its NRMSE SNR] = lps_ist(E,T,dn,P1 * C * lambdaList,P1 * lambdaSist,tol,MaxIters,Mtrue);
    ISTn.sigma(i) = sigma;
    ISTn.P(i) = P1;
    ISTn.L{i} = L;
    ISTn.S{i} = S;
    ISTn.cost{i} = cost;
    ISTn.deltaM{i} = deltaM;
    ISTn.time{i} = time;
    ISTn.its(i) = its;
    ISTn.NRMSE{i} = NRMSE;
    ISTn.NRMSEend(i) = NRMSE(end);
    ISTn.SNR{i} = SNR;
    ISTn.SNRend(i) = SNR(end);
    
    % Our SLRN algo
    P2 = Plist(Pbest_slrn(i));
    [L S cost deltaM time its NRMSE SNR] = SLRN_MRI(E,T,dn,rslrn,P1 * lambdaSslrn,tol,MaxIters,Mtrue);
    SLRNn.sigma(i) = sigma;
    SLRNn.P(i) = P1;
    SLRNn.L{i} = L;
    SLRNn.S{i} = S;
    SLRNn.cost{i} = cost;
    SLRNn.deltaM{i} = deltaM;
    SLRNn.time{i} = time;
    SLRNn.its(i) = its;
    SLRNn.NRMSE{i} = NRMSE;
    SLRNn.NRMSEend(i) = NRMSE(end);
    SLRNn.SNR{i} = SNR;
    SLRNn.SNRend(i) = SNR(end);
end
%--------------------------------------------------------------------------

%% Plot NRMSEs over a specific ROI

% Subset of [1 nx 1 ny] = [1 128 1 128]
ROI = [1 128 1 128];

% Ground truth
MtrueROI = abs(Mtrue(ROI(1):ROI(2),ROI(3):ROI(4),:));
NRMSE_fcn = @(M,Mtrue) 100 * norm(M(:) - Mtrue(:),2) / norm(Mtrue(:),2);

% Compute NRMSE
nrmse_ist = nan(NSigma,nt);
nrmse_slrn = nan(NSigma,nt);
for i = 1:NSigma
    Mist = abs(ISTn.L{i}(ROI(1):ROI(2),ROI(3):ROI(4),:) + ISTn.S{i}(ROI(1):ROI(2),ROI(3):ROI(4),:));
    Mslrn = abs(SLRNn.L{i}(ROI(1):ROI(2),ROI(3):ROI(4),:) + SLRNn.S{i}(ROI(1):ROI(2),ROI(3):ROI(4),:));
    for j = 1:nt
        nrmse_ist(i,j) = NRMSE_fcn(Mist(:,:,j),MtrueROI(:,:,j));
        nrmse_slrn(i,j) = NRMSE_fcn(Mslrn(:,:,j),MtrueROI(:,:,j));
    end
end

% Plot NRMSEs
cm = hsv(NSigma);
figure;
subplot(1,2,1);
for i = 1:NSigma
    plot(1:nt,nrmse_ist(i,:),'-o','Color',cm(i,:));
    hold on;
end
title('IST');
xlabel('Frame');
ylabel('NRMSE');
subplot(1,2,2);
hold on;
for i = 1:NSigma
    plot(1:nt,nrmse_slrn(i,:),'-o','Color',cm(i,:));
    hold on;
end
title('SLRN');
xlabel('Frame');
ylabel('NRMSE');

% Plot NRMSE percent decrease 
figure;
for i = 1:NSigma
    plot(1:nt,100 * (nrmse_ist(i,:) - nrmse_slrn(i,:)) ./ nrmse_ist(i,:),'-o','Color',cm(i,:));
    hold on;
end
title('Percent Reduction in NRMSE');
xlabel('Frame');
ylabel('100 * (NRMSE(IST) - NRMSE(SLRN)) / NRMSE(IST)');
