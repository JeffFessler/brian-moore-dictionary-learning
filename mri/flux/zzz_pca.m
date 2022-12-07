%% Compare rank-robustness of OptShrink and PCA

%--------------------------------------------------------------------------
% Knobs
%--------------------------------------------------------------------------
% Data generation
nc = 8;
nt = 21;
SNR = 55;
seed = 1;
%inpath = 'dynobj.mat';
%inpath = 'dynobj2.mat';
inpath = 'dynobj3.mat';

%{
% L + S
proxS = 'soft';
algoStr = 'L+S PGM';
NRMSEpath = './nrmse_3/NRMSE_LpS_soft_3.mat';
T = 1; % {1,TempFFT(2)}
datapath = './data/pca_comp_LpS_soft_3.mat';
%}

%{
% L + S
proxS = 'soft';
algoStr = 'L+S PGM';
NRMSEpath = './nrmse_3/NRMSE_LpS_soft_T_3.mat';
T = TempFFT(2); % {1,TempFFT(2)}
datapath = './data/pca_comp_LpS_soft_T_3.mat';
%}

%{
% L + S
proxS = 'mixed21';
algoStr = 'L+S PGM';
NRMSEpath = './nrmse_3/NRMSE_LpS_mixed21_3.mat';
T = 1; % {1,TempFFT(2)}
datapath = './data/pca_comp_LpS_mixed21_3.mat';
%}

%{
% L & S
proxS = 'soft';
algoStr = 'L&S ADMM';
NRMSEpath = './nrmse_3/NRMSE_LaS_soft_3.mat';
T = 1; % {1,TempFFT(2)}
datapath = './data/pca_comp_LaS_soft_3.mat';
%}

%{
% L & S
proxS = 'soft';
algoStr = 'L&S ADMM';
NRMSEpath = './nrmse_3/NRMSE_LaS_soft_T_3.mat';
T = TempFFT(2); % {1,TempFFT(2)}
datapath = './data/pca_comp_LaS_soft_T_3.mat';
%}

%{
% L & S
proxS = 'mixed21';
algoStr = 'L&S ADMM';
NRMSEpath = './nrmse_3/NRMSE_LaS_mixed21_3.mat';
T = 1; % {1,TempFFT(2)}
datapath = './data/pca_comp_LaS_mixed21_3.mat';
%}

% Test parameters
r = 1:10; % Test ranks
metric = [1 5]; % {1 = X, 2 = lesion #1, ..., 5 = ROI}
%--------------------------------------------------------------------------

% Load data
load(NRMSEpath);
gtData = load(inpath);
Xtrue_full = gtData.dyn_obj;
dce = gtData.dce;
clear gtData;

% Parse parameters
isLpS = strcmp(algoStr,'L+S PGM');
Nr = length(r);
Nm = length(metric);
N3 = 1 + isLpS;

% Get FFT NRMSE data
nrmsefft = zeros(1,Nm,1);
[~,nc_fft_idx] = min(abs(labels_fft.nc - nc));
[~,nt_fft_idx] = min(abs(labels_fft.nt - nt));
[~,SNR_fft_idx] = min(abs(labels_fft.SNR - SNR));
for j = 1:Nm
    nrmsefft(1,j,1) = NRMSE_fft{metric(j)}(nc_fft_idx,nt_fft_idx,SNR_fft_idx);
end

% Get optimal SVT parameters
[~,nc_svt_idx] = min(abs(labels_svt.nc - nc));
[~,nt_svt_idx] = min(abs(labels_svt.nt - nt));
[~,SNR_svt_idx] = min(abs(labels_svt.SNR - SNR));
inds1 = idxopt_svt{metric(j)}{nc_svt_idx,nt_svt_idx,SNR_svt_idx};
lambdaL = labels_svt.lambdaL(inds1(1));
lambdaS1 = labels_svt.lambdaS(inds1(2));

% Get optimal OptShrink lambdaS
[~,nc_opt_idx] = min(abs(labels_opt.nc - nc));
[~,nt_opt_idx] = min(abs(labels_opt.nt - nt));
[~,SNR_opt_idx] = min(abs(labels_opt.SNR - SNR));
inds2 = idxopt_opt{metric(j)}{nc_opt_idx,nt_opt_idx,SNR_opt_idx};
lambdaS2 = labels_opt.lambdaS(inds2(2));

% Perform simulations
nrmsesvt = rand(1,Nm,N3);
nrmseopt = rand(Nr,Nm,N3);
nrmsepca = rand(Nr,Nm,N3);
for j = 1:Nm
    % SVT simulation
    [~,recon_svt,~,Xfft,mask,ROIs,~,splineInterp] = run_svt_algo(nc,nt,SNR,seed,lambdaL,proxS,lambdaS1,algoStr,T,Xtrue_full,dce);
    Xhat_svt = splineInterp(abs(embed(recon_svt.X,mask)));
    nrmsesvt(1,j,1) = computeNRMSE(Xhat_svt,Xtrue_full,metric(j),ROIs);
    if isLpS
        Lhat_svt = splineInterp(abs(embed(recon_svt.L,mask)));
        nrmsesvt(1,j,2) = computeNRMSE(Lhat_svt,Xtrue_full,metric(j),ROIs);
    end
    
    % Loop over ranks
    for i = 1:Nr
        % OptShrink simulation
        [~,recon_opt,~,Xfft,mask,ROIs,~,splineInterp] = run_optshrink_algo(nc,nt,SNR,seed,r(i),proxS,lambdaS2,algoStr,T,Xtrue_full,dce);
        Xhat_opt = splineInterp(abs(embed(recon_opt.X,mask)));
        nrmseopt(i,j,1) = computeNRMSE(Xhat_opt,Xtrue_full,metric(j),ROIs);
        if isLpS
            Lhat_opt = splineInterp(abs(embed(recon_opt.L,mask)));
            nrmseopt(i,j,2) = computeNRMSE(Lhat_opt,Xtrue_full,metric(j),ROIs);
        end
        
        % PCA simulation
        [~,recon_pca,~,Xfft,mask,ROIs,~,splineInterp] = run_pca_algo(nc,nt,SNR,seed,r(i),proxS,lambdaS2,algoStr,T,Xtrue_full,dce);
        Xhat_pca = splineInterp(abs(embed(recon_pca.X,mask)));
        nrmsepca(i,j,1) = computeNRMSE(Xhat_pca,Xtrue_full,metric(j),ROIs);
        if isLpS
            Lhat_pca = splineInterp(abs(embed(recon_pca.L,mask)));
            nrmsepca(i,j,2) = computeNRMSE(Lhat_pca,Xtrue_full,metric(j),ROIs);
        end
    end
end

% Save data
if exist('datapath','var') && ~isempty(datapath)
    save(datapath,'nrmsefft','nrmsesvt','nrmseopt','nrmsepca','lambdaL','r','lambdaS1','lambdaS2','metric','nc','nt','SNR','seed','proxS','algoStr');
end

fprintf('DONE\n');

%% Plot results

% Knobs
datapath = './data/pca_comp_LpS_soft_3.mat';
%datapath = './data/pca_comp_LpS_soft_T_3.mat';
%datapath = './data/pca_comp_LpS_mixed21_3.mat';
%datapath = './data/pca_comp_LaS_soft_3.mat';
%datapath = './data/pca_comp_LaS_soft_T_3.mat';
%datapath = './data/pca_comp_LaS_mixed21_3.mat';
gapy = 0.1;
pos = [597 370 729 462];
%pos = [597 473 728 257];

% Load data
load(datapath);
isLpS = strcmp(algoStr,'L+S PGM');
Nr = length(r);
Nm = length(metric);
N3 = 1 + isLpS;

% Plot results
figure;
cm = linspecer(4);
phndl = zeros(1,4);
for j = 1:Nm
    % Get y-label string
    switch metric(j)
            case 1
                region = 'X';
            case 2
                region = 'ROI1';
            case 3
                region = 'ROI2';
            case 4
                region = 'ROI3';
            case 5
                region = 'ROI';
    end
    for i = 1:N3
        % Plot results (i,j)
        subplot(N3,Nm,j + (i - 1) * Nm);
        hold on;
        phndl(1) = plot(r,nrmseopt(:,j,i),'-o','Color',cm(1,:));
        phndl(2) = plot(r,nrmsepca(:,j,i),'-o','Color',cm(2,:));
        phndl(3) = plot(r([1 Nr]),nrmsesvt(1,j,i) * [1 1],'-','Color',cm(3,:));
        phndl(4) = plot(r([1 Nr]),nrmsefft(1,j,1) * [1 1],'-','Color',cm(4,:));
        xlabel('r');
        ylabel(sprintf('NRMSE(%s)',region));
        if (i == 1)
            title('Xhat');
        else
            title('Lhat');
        end
        axis tight;
        padAxis(gca,0,'x');
        padAixs(gca,gapy,'y');
    end
end
legend(phndl,'OptShrink(r)','PCA(r)','SVT','FFT','Location','Best');
set(gcf,'Position',pos);

%export_fig -pdf -transparent pca_comp_LpS_soft_3
%export_fig -pdf -transparent pca_comp_LpS_soft_T_3
%export_fig -pdf -transparent pca_comp_LpS_mixed21_3
%export_fig -pdf -transparent pca_comp_LaS_soft_3
%export_fig -pdf -transparent pca_comp_LaS_soft_T_3
%export_fig -pdf -transparent pca_comp_LaS_mixed21_3
