%
% Multicoil cardiac perfusion MRI reconstruction algorithms in the presence
% of noise
%

%--------------------------------------------------------------------------
% Knobs
%--------------------------------------------------------------------------
MinSigma = 1;
MaxSigma = 25;
NSigma = 13;

MinP = 0.1;
MaxP = 10;
NParams = 5;

lambdaList = 0.01;
lambdaSist = 0.01;

rslrn = 1;
lambdaSslrn = 0.1;

tol = 0.0025;
MaxIters = 40;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Simulation lists
%--------------------------------------------------------------------------
sigmalist = linspace(MinSigma,MaxSigma,NSigma);
Plist = logspace(log10(MinP),log10(MaxP),NParams);
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
load('LS_ist.mat');
Mtrue = Ltrue + Strue;
%load('LS_SLRN.mat');
%Mtrue = Ltrue + Strue;

% True sub-sampled k-space data
d = E * Mtrue;
idx = (d ~= 0);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Data storage
%--------------------------------------------------------------------------
ISTn.sigma = nan(NSigma,NParams);
ISTn.P = nan(NSigma,NParams);
ISTn.L = cell(NSigma,NParams);
ISTn.S = cell(NSigma,NParams);
ISTn.cost = cell(NSigma,NParams);
ISTn.deltaM = cell(NSigma,NParams);
ISTn.time = cell(NSigma,NParams);
ISTn.its = nan(NSigma,NParams);
ISTn.NRMSE = cell(NSigma,NParams);
ISTn.NRMSEend = nan(NSigma,NParams);
ISTn.SNR = cell(NSigma,NParams);
ISTn.SNRend = nan(NSigma,NParams);

SLRNn.sigma = nan(NSigma,NParams);
SLRNn.P = nan(NSigma,NParams);
SLRNn.L = cell(NSigma,NParams);
SLRNn.S = cell(NSigma,NParams);
SLRNn.cost = cell(NSigma,NParams);
SLRNn.deltaM = cell(NSigma,NParams);
SLRNn.time = cell(NSigma,NParams);
SLRNn.its = nan(NSigma,NParams);
SLRNn.NRMSE = cell(NSigma,NParams);
SLRNn.NRMSEend = nan(NSigma,NParams);
SLRNn.SNR = cell(NSigma,NParams);
SLRNn.SNRend = nan(NSigma,NParams);
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
    
    for j = 1:NParams
        % Display progress
        fprintf('\nTrial %i/%i\n',(i - 1) * NParams + j,NSigma * NParams);
        
        % Parameter multiplier
        P = Plist(j);
        
        % Otazo's L + S algo
        [L S cost deltaM time its NRMSE SNR] = lps_ist(E,T,dn,P * C * lambdaList,P * lambdaSist,tol,MaxIters,Mtrue);
        ISTn.sigma(i,j) = sigma;
        ISTn.P(i,j) = P;
        ISTn.L{i,j} = L;
        ISTn.S{i,j} = S;
        ISTn.cost{i,j} = cost;
        ISTn.deltaM{i,j} = deltaM;
        ISTn.time{i,j} = time;
        ISTn.its(i,j) = its;
        ISTn.NRMSE{i,j} = NRMSE;
        ISTn.NRMSEend(i,j) = NRMSE(end);
        ISTn.SNR{i,j} = SNR;
        ISTn.SNRend(i,j) = SNR(end);
        
        % Our SLRN algo
        [L S cost deltaM time its NRMSE SNR] = SLRN_MRI(E,T,dn,rslrn,P * lambdaSslrn,tol,MaxIters,Mtrue);
        SLRNn.sigma(i,j) = sigma;
        SLRNn.P(i,j) = P;
        SLRNn.L{i,j} = L;
        SLRNn.S{i,j} = S;
        SLRNn.cost{i,j} = cost;
        SLRNn.deltaM{i,j} = deltaM;
        SLRNn.time{i,j} = time;
        SLRNn.its(i,j) = its;
        SLRNn.NRMSE{i,j} = NRMSE;
        SLRNn.NRMSEend(i,j) = NRMSE(end);
        SLRNn.SNR{i,j} = SNR;
        SLRNn.SNRend(i,j) = SNR(end);
    end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Save data
%--------------------------------------------------------------------------
%save('C:\Users\brimoor\Desktop\NoisyCardiacPerfComparison.mat','ISTn','SLRNn','-v7.3');
save('NoisyCardiacPerfComparison.mat','ISTn','SLRNn','-v7.3');
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
[~,Pbest_ist] = max(ISTn.SNRend,[],2);
[~,Pbest_slrn] = max(SLRNn.SNRend,[],2);
MaxIts = max(max(ISTn.its(sub2ind([NSigma NParams],(1:NSigma)',Pbest_ist))),max(SLRNn.its(sub2ind([NSigma NParams],(1:NSigma)',Pbest_slrn))));
cm = hsv(NSigma);

figure;
%NMarkers = 10;
gap1 = 0.1;
gap2 = 0.2;
ha = tight_subplot(1,2,[0 gap2],[gap1 gap1],[gap1 gap1]);

% NRMSE
axes(ha(1));
phndl = zeros(1,2 * NSigma);
pstr = cell(1,2 * NSigma);
for i = NSigma:-1:1
    phndl(2 * i - 1) = plot(1:ISTn.its(i,Pbest_ist(i)),ISTn.NRMSE{i,Pbest_ist(i)},':*','Color',cm(i,:));
    %phndl(2 * i - 1) = line_fewer_markers(1:ISTn.its(i,Pbest_ist(i)),ISTn.NRMSE{i,Pbest_ist(i)},NMarkers,':*','Color',cm(i,:));
    pstr{2 * i - 1} = sprintf('IST - \\sigma = %.1f',sigmalist(i));
    hold on;
    phndl(2 * i) = plot(1:SLRNn.its(i,Pbest_slrn(i)),SLRNn.NRMSE{i,Pbest_slrn(i)},'-o','Color',cm(i,:));
    %phndl(2 * i) = line_fewer_markers(1:SLRNn.its(i,Pbest_slrn(i)),SLRNn.NRMSE{i,Pbest_slrn(i)},NMarkers,'-o','Color',cm(i,:));
    pstr{2 * i} = sprintf('SLRN - \\sigma = %.1f',sigmalist(i));
    hold on;
end
title('NRMSE (%)');
xlabel('Iteration');
legend(phndl,pstr{:});
set(gca,'XLim',[1 MaxIts]);

% SNR
axes(ha(2));
for i = NSigma:-1:1
    plot(1:ISTn.its(i,Pbest_ist(i)),ISTn.SNR{i,Pbest_ist(i)},':*','Color',cm(i,:));
    %line_fewer_markers(1:ISTn.its(i,Pbest_ist(i)),ISTn.SNR{i,Pbest_ist(i)},NMarkers,':*','Color',cm(i,:));
    hold on;
    plot(1:SLRNn.its(i,Pbest_slrn(i)),SLRNn.SNR{i,Pbest_slrn(i)},'-o','Color',cm(i,:));
    %line_fewer_markers(1:SLRNn.its(i,Pbest_slrn(i)),SLRNn.SNR{i,Pbest_slrn(i)},NMarkers,'-o','Color',cm(i,:));
    hold on;
end
title('SNR (dB)');
xlabel('Iteration');
set(gca,'XLim',[1 MaxIts]);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Extract NRMSEs and SNRs from ISTn and SLRNn
%--------------------------------------------------------------------------
nrmse_ist = inf(NSigma,1);
snr_ist = zeros(NSigma,1);
nrmse_slrn = inf(NSigma,1);
snr_slrn = zeros(NSigma,1);

for i = 1:NSigma
    nrmse_ist(i) = min(ISTn.NRMSE{i,Pbest_ist(i)});
    snr_ist(i)   = max(ISTn.SNR{i,Pbest_ist(i)});
    
    nrmse_slrn(i) = min(SLRNn.NRMSE{i,Pbest_slrn(i)});
    snr_slrn(i)   = max(SLRNn.SNR{i,Pbest_slrn(i)});
    
    %{
    for j = 1:NParams
        nrmse_ist(i) = min(min(ISTn.NRMSE{i,j}),nrmse_ist(i));
        snr_ist(i) = max(max(ISTn.SNR{i,j}),snr_ist(i));
        
        nrmse_slrn(i) = min(min(SLRNn.NRMSE{i,j}),nrmse_slrn(i));
        snr_slrn(i) = max(max(SLRNn.SNR{i,j}),snr_slrn(i));
    end
    %}
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot more results
%--------------------------------------------------------------------------
figure;

subplot(1,2,1);
plot(sigmalist,nrmse_ist,'b-o');
hold on;
plot(sigmalist,nrmse_slrn,'r-o');
title('NRMSE (%)');
xlabel('\sigma');

subplot(1,2,2);
plot(sigmalist,snr_ist,'b-o');
hold on;
plot(sigmalist,snr_slrn,'r-o');
title('SNR (dB)');
xlabel('\sigma');
%--------------------------------------------------------------------------
