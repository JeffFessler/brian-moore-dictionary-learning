%% Noiseless
%
% L+S reconstruction of undersampled multicoil cardiac perfusion MRI
%
% Ricardo Otazo (2013)
%

% Load undersampled data 
load('cardiac_perf_R8.mat');
[nx ny nt nc] = size(kdata);
nd = nx * ny;

%--------------------------------------------------------------------------
% L + S reconstruction
%--------------------------------------------------------------------------
E = Emat_xyt(kdata(:,:,:,1) ~= 0,b1);
T = TempFFT(3);
d = kdata; % sub-sampled k-space data

% Regularization parameters
C = max(svd(reshape(E' * d,[nd nt]))); % scaling constant for Otazo's SVT
lambdaL_ist = 0.01;
lambdaS_ist = 0.01;

r = 1;
lambdaS = 0.1;

%r = 5;
%lambdaS = 0.001;

%lambda = 1 / sqrt((nnz(d) / numel(d)) * nd);
%epsilon = 18.25;
%NInnerIters = 7;

tol = 0.0025;
MaxIters = 50;

[L S cost deltaM time its] = lps_ist(E,T,d,C * lambdaL_ist,lambdaS_ist,tol,MaxIters);
%[L S cost deltaM time its] = SLRN_MRI(E,T,d,r,lambdaS,tol,MaxIters);
%[L S cost deltaM time its] = RPCA_DN_MRI(E,T,d,lambda,epsilon,tol,MaxIters,NInnerIters);

M = L + S;
%--------------------------------------------------------------------------

% Extract 4 frames
Md = M(33:96,33:96,2);
Md = cat(2,Md,M(33:96,33:96,8));
Md = cat(2,Md,M(33:96,33:96,14));
Md = cat(2,Md,M(33:96,33:96,24));
Ld = L(33:96,33:96,2);
Ld = cat(2,Ld,L(33:96,33:96,8));
Ld = cat(2,Ld,L(33:96,33:96,14));
Ld = cat(2,Ld,L(33:96,33:96,24));
Sd = S(33:96,33:96,2);
Sd = cat(2,Sd,S(33:96,33:96,8));
Sd = cat(2,Sd,S(33:96,33:96,14));
Sd = cat(2,Sd,S(33:96,33:96,24));

%--------------------------------------------------------------------------
% Plot images
%--------------------------------------------------------------------------
figure;

subplot(3,1,1);
imshow(abs(Md),[0 1]);
title('L + S');

subplot(3,1,2);
imshow(abs(Ld),[0 1]);
title('L');

subplot(3,1,3);
imshow(abs(Sd),[0 1]);
title('S');
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot algorithm convergence
%--------------------------------------------------------------------------
figure;

subplot(1,3,1);
plot(1:its,cost,'b-o');
title('Cost Function');
xlabel('Iteration');
set(gca,'XLim',[1 its]);

subplot(1,3,2);
semilogy(1:its,deltaM,'b-o');
title('Convergence Criterion');
xlabel('Iteration');
set(gca,'XLim',[1 its]);

subplot(1,3,3);
semilogy(1:its,cumsum(time),'b-o');
title('Time (Seconds)');
xlabel('Iteration');
set(gca,'XLim',[1 its]);
%--------------------------------------------------------------------------

%% Noisy
%
% L+S reconstruction of undersampled multicoil cardiac perfusion MRI
%
% Ricardo Otazo (2013)
%

% Load undersampled data 
load('cardiac_perf_R8.mat');
[nx ny nt nc] = size(kdata);
nd = nx * ny;

% k-spase noise variance
sigma = 0 / sqrt(nd);

% Constants
tol = 0.0025;
MaxIters = 50;

% Format data
E = Emat_xyt(kdata(:,:,:,1) ~= 0,b1);
T = TempFFT(3);
d = kdata; % sub-sampled k-space data

% Add noise
idx = (d ~= 0);
noise = (sigma / sqrt(2)) * (randn(size(d)) + 1i * randn(size(d)));
dn = zeros(size(d));
dn(idx) = d(idx) + noise(idx);

%--------------------------------------------------------------------------
% IST
%--------------------------------------------------------------------------
% Regularization parameters
C = max(svd(reshape(E' * dn,[nd nt]))); % scaling constant for Otazo's SVT
lambdaL_ist = 0.01;
lambdaS_ist = 0.01;

% True images
load('LS_ist.mat');
Mtrue_ist = Ltrue + Strue;

% Perform reconstruction
[L_ist S_ist cost_ist deltaM_ist time_ist its_ist NRMSE_ist SNR_ist] = lps_ist(E,T,dn,C * lambdaL_ist,lambdaS_ist,tol,MaxIters,Mtrue_ist);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% SLRN
%--------------------------------------------------------------------------
% Regularization parameters
r_slrn = 1;
lambdaS_slrn = 0.1;

% True images
load('LS_SLRN.mat');
Mtrue_slrn = Ltrue + Strue;

% Perform reconstruction
[L_slrn S_slrn cost_slrn deltaM_slrn time_slrn its_slrn NRMSE_slrn SNR_slrn] = SLRN_MRI(E,T,dn,r_slrn,lambdaS_slrn,tol,MaxIters,Mtrue_slrn);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot NRMSEs and SNRs
%--------------------------------------------------------------------------
figure;

subplot(1,2,1);
plot(1:its_ist,NRMSE_ist,'b-o');
hold on;
plot(1:its_slrn,NRMSE_slrn,'r-o');
title('NRMSE (%)');
xlabel('Iteration');
set(gca,'XLim',[1 max(its_ist,its_slrn)]);

subplot(1,2,2);
plot(1:its_ist,SNR_ist,'b-o');
hold on;
plot(1:its_slrn,SNR_slrn,'r-o');
title('SNR (dB)');
xlabel('Iteration');
set(gca,'XLim',[1 max(its_ist,its_slrn)]);
%--------------------------------------------------------------------------

%%

% Extract 4 frames
Md = M(33:96,33:96,2);
Md = cat(2,Md,M(33:96,33:96,8));
Md = cat(2,Md,M(33:96,33:96,14));
%Md = cat(2,Md,M(33:96,33:96,24));

% Plot images
figure;
imshow(abs(Md),[0 1]);
%imshow(abs(M(:,:,14)),[0 1]);
