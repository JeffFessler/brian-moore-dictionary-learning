%% Compute fully-sampled ground truth

% Knobs
inpath = '../cardiac_perf_full.mat';
outpath = './otazo_true.mat';

% Load data
load(inpath);
[ny nx nt nc] = size(kdata);
nd = ny * nx;
samp = b1;
mask = kdata(:,:,:,1) ~= 0;
Ytrue = kdata;
A = Emat_xyt(mask,samp);

% Compute norm(A) for fully-sampled data
% NOTE: This is an upperbound on norm(A) when undersampling exists
%{
m = ny * nx * nt;
AAtfcn = @(x) reshape(A' * (A * reshape(x,[ny nx nt])),m,1);
opts = struct('issym',true,'isreal',false,'maxit',5,'disp',1);
normA = sqrt(eigs(AAtfcn,m,1,'LM',opts)) %#ok
%}
normA = 1;

% Compute ground truth
% NOTE: This works because A' * A == I
Xtrue = A' * Ytrue;

% Save results
if exist('outpath','var') && ~isempty(outpath)
    save(outpath,'Ytrue','Xtrue','samp','normA','ny','nx','nd','nt','nc');
end

%% Generate column sampling distribution

% Knobs
inpath = './otazo.mat';
%outpath = './otazo_true.mat';
w = 1; % Smoothing filter width
kappa = 0.75; % Constant percentage
beta = 0.05; % Minimum sampling probability
Ntrials = 40;

% Load raw 8x undersampled data
load(inpath);
Ysamp = abs(Y(:,:,:,1)) ~= 0;

% Generate sampling mask image
[nh nw] = bestSubplotShape(nt);
C = cell(nh,nw);
for i = 1:nt
    C{i} = Ysamp(:,:,i);
end
M = cell2mov(C,30,1);

% Parameteric sampling probabilities
nhalf = 0.5 * nx;
x1 = round(kappa * nhalf);
x2 = nhalf + 0.5;
%2 = @(x) (x >= x1) .* ((1 - beta) * ((x - x1) / (x2 - x1)).^2) + beta;
%ps = p([1:nhalf,nhalf:-1:1]);
p = @(x) (x >= x1) .* (beta * exp(-log(beta) * (x - x1) / (x2 - x1))) + beta .* (x < x1);
ps = p([1:nhalf,nhalf:-1:1]);

% Otazo sampling probabilities
ps1 = squeeze(sum(double(Ysamp(1,:,:)),3)) / nt;
if w > 1
    % Smooth distribution
    ps1 = conv(ps1,ones(1,w) / w,'same');
    %ps1 = medfilt1(ps1,w);
end

% Empirical samples from parametric curve
ps2 = zeros(Ntrials,nx);
for i = 1:Ntrials
    ps2(i,sampleCols(round(0.125 * nx),ps)) = 1;
end
ps2 = mean(ps2,1);

% Save empirical sampling probabilities
if exist('outpath','var') && ~isempty(outpath)
    save(outpath,'ps','-append');
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
% Sampling masks
figure;
imshow(M,[]);
%{
export_fig -png -transparent samp_im
%}

% Column sampling probabilities
figure;
PS = [ps; ps1; ps2];
pstr = {'parametric','otazo','samples'};
[N, nx] = size(PS);
cm = linspecer(N);
subplot(1,2,1);
phndl = zeros(1,N);
for i = 1:N
    phndl(i) = plot(PS(i,:),'-o','Color',cm(i,:));
    hold on;
end
legend(phndl,pstr);
xlabel('Column');
ylabel('PDF');
axis tight; padAxis();
subplot(1,2,2);
for i = 1:N
    CMFplot(CMFgen({1:nx,cumsum(PS(i,:))}),'Color',cm(i,:));
    hold on;
end
title('');
xlabel('Column');
ylabel('CDF');
axis tight; padAxis();

mtitle(gcf,'Pr(sampled)');

%{
export_fig -pdf -transparent samp_prob
%}
%--------------------------------------------------------------------------

%% Format Otazo's 8x undersampled data

% Knobs
inpath = '../cardiac_perf_R8.mat';
%outpath = 'otazo.mat';

% Load raw data
load(inpath);
[ny nx nt nc] = size(kdata);
nd = ny * nx;
samp = b1;
mask = kdata(:,:,:,1) ~= 0;
Y = kdata;

% Compute norm(A)
%{
A = Emat_xyt(mask,samp);
m = ny * nx * nt;
AAtfcn = @(x) reshape(A' * (A * reshape(x,[ny nx nt])),m,1);
opts = struct('issym',true,'isreal',false,'maxit',5,'disp',1);
normA = sqrt(eigs(AAtfcn,m,1,'LM',opts)) %#ok
%}
normA = 1;

% Initial iterate: naive
A = Emat_xyt(mask,samp);
X0 = A' * Y;

%{
% Initial iterate: LS
A = Emat_xyt(mask,samp);
X0 = A' * Y;
tau = 1.5; % tau in (0,2 / sigma1(A)^2)
Niters = 20;
for i = 1:Niters
    X0last = X0;
    X0 = X0 - tau * (A' * (A * X0 - Y));
    delta = norm(X0(:) - X0last(:)) / norm(X0last(:));
    fprintf('Iteration %i/%i (delta = %.3g)\n',i,Niters,delta);
end
%}

%{
% Initial iterate: fft
kfull = data_share_fill(reshape(Y(abs(Y) ~= 0),[],nc),mask); % Undersampled
%kfull = kdata; % Fully sampled
xc = permute(ir_fftshift2(ifft2(ir_fftshift2(kfull))),[1 2 4 3]);
Xfft = squeeze(sum(xc .* conj(repmat(samp,[1 1 1 nt])),3)) ./ ...
       repmat(sum(abs(samp).^2,3),[1 1 nt]);
X0 = flipdim(flipdim(Xfft,1),2);
%}

% Save data
if exist('outpath','var') && ~isempty(outpath)
    save(outpath,'Y','normA','mask','samp','X0','ny','nx','nd','nt','nc');
end

%% Generate masks

% Knobs
path = 'otazo.mat';
tau = 0.24;

% Load data
load(path);
Xhat = Xhat / max(Xhat(:));

% Region #1
%M1 = Xhat(:,:,7) > tau;
%P = Paint(M1);
%M1 = P.Close(); save('mask1.mat','M1');

% Region #2
%M2 = Xhat(:,:,13) > tau;
%P = Paint(M2);
%M2 = P.Close(); save('mask2.mat','M2');

% Region #3
%M3 = Xhat(:,:,19) > tau;
%P = Paint(M3);
%M3 = P.Close(); save('mask3.mat','M3');

% Region #4
%M4 = Xhat(:,:,21) > tau;
%P = Paint(M4);
%M4 = P.Close(); save('mask4.mat','M4');

% Region #5
%M5 = Xhat(:,:,14) > tau;
%P = Paint(M5);
%M5 = P.Close(); save('mask5.mat','M5');

% Get static values
%{
w = 5;
H = ones(w,w) / w^2;
X1 = imfilter(Xhat(:,:,1),H,'replicate','same');
imshow(X1);
[x y] = ginput();
inds = sub2ind(size(X1),round(y),round(x));
static_vals = X1(inds);
%}
%static_vals = [0.17029; 0.11321; 0.037934];

% Save data
%save('otazo_masks.mat','Xhat','M1','M2','M3','M4','M5','static_vals');

%% Plot contrast curves

% Knobs
path = 'otazo_masks.mat';

% Load data
load(path);
Ysamp = M1 | M2 | M3 | M4 | M5;
cm = linspecer(5);
[ny nx nt] = size(Xhat);

% Scale to [0 1]
Xhat = Xhat / max(Xhat(:));

% Mean image
Xbar = mean(Xhat,3);
XX = repmat(im2uint8(Xbar),[1 1 3]);

% Colored mask image
CM = uint8(floor(255 * cm));
%MM = 255 * ones(ny,nx,3,'uint8');% White background
%MM = zeros(ny,nx,3,'uint8'); % Black background
MM = XX; % Mean background
for i = 1:5
    Mi = eval(sprintf('M%i;',i));
    for j = 1:3
        MMj = MM(:,:,j);
        MMj(Mi) = CM(i,j);
        MM(:,:,j) = MMj;
    end
end

% Constrast curves
Tbar = zeros(5,nt);
Tmin = zeros(5,nt);
Tmax = zeros(5,nt);
for i = 1:5
    Mi = eval(sprintf('M%i;',i));
    for j = 1:nt
        Xj = Xhat(:,:,j);
        Xij = Xj(Mi);
        Tbar(i,j) = mean(Xij(:));
        Tmin(i,j) = min(Xij(:));
        Tmax(i,j) = max(Xij(:));
    end
end

% Smooth curves
w = 5; % Filter width
h = ones(1,w) / w;
for i = 1:5
    Tmin(i,:) = imfilter(Tmin(i,:),h,'replicate','same');
    Tmax(i,:) = imfilter(Tmax(i,:),h,'replicate','same');
end

%--------------------------------------------------------------------------
% Plot results #1
%--------------------------------------------------------------------------
figure('Color','w');

% Mean image
subplot(2,2,1);
imshow(XX);

% Colored mask image
subplot(2,2,2);
imshow(MM);

% Mean trajectories
subplot(2,2,3);
hold on;
for i = 1:5
    plot(1:nt,Tbar(i,:),'-o','Color',cm(i,:));
end
for j = 1:length(static_vals)
    plot(1:nt,static_vals(j) * ones(1,nt),'k:');
end
xlabel('Frame');
title('Mean Trajectories');
set(gca,'XLim',[1 nt]);
set(gca,'YLim',[0 max(Xhat(:))]);
box on;

% Max trajectories
subplot(2,2,4);
ph = zeros(1,5);
for i = 1:5
    ph(i) = plot(1:nt,Tmax(i,:),'-o','Color',cm(i,:));
    hold on;
end
for j = 1:length(static_vals)
    plot(1:nt,static_vals(j) * ones(1,nt),'k:');
end
xlabel('Frame');
title('Max Trajectories');
set(gca,'XLim',[1 nt]);
set(gca,'YLim',[0 max(Xhat(:))]);
box on;

% Save figure
set(gcf,'Position',[274 94 1375 1015]);
%export_fig -pdf -transparent otazo_traj1
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot results #2
%--------------------------------------------------------------------------
figure('Color','w');

% Images
subplot(3,5,1:10);
imshow(cell2mov({XX MM},64,255));

% Trajectories
Nbars = 5; inds = round(linspace(2,nt - 1,Nbars));
for i = 1:5
    subplot(3,5,10 + i);
    hold on;
    xvec = 1:nt;
    plot(xvec,Tbar(i,:),'-','Color',cm(i,:));
    plot(xvec,Tmax(i,:),':','Color',cm(i,:));
    plot(xvec,Tmin(i,:),':','Color',cm(i,:));
    Li = Tbar(i,inds) - Tmin(i,inds);
    Ui = Tmax(i,inds) - Tbar(i,inds);
    errorbar(xvec(inds),Tbar(i,inds),Li,Ui,'.','Color',cm(i,:));
    for j = 1:length(static_vals)
        plot(1:nt,static_vals(j) * ones(1,nt),'k:');
    end
    xlabel('Frame');
    title(sprintf('Region %i',i));
    set(gca,'XLim',[1 nt]);
    set(gca,'YLim',[0 max(Xhat(:))]);
    box on;
end

% Save figure
set(gcf,'Position',[317 260 1289 682]);
%export_fig -pdf -transparent otazo_traj
%--------------------------------------------------------------------------

%% Plot ROIs

% Knobs
inpath = './otazo_true.mat';
typeLabels = {'Whole image','Heart-only','Body-only'};

% Load data
load(inpath);

% Get images
Xt = abs(Xtrue);
Xm = mean(Xt,3);
lim = [min(Xt(:)), max(Xt(:))];
X1 = Xm; X1(~ROIs{1}) = lim(2);
X2 = Xm; X2(~ROIs{2}) = lim(2);

% Plot images
figure;
subplot(1,3,1); imshow(Xm,lim); title(typeLabels{1});
subplot(1,3,2); imshow(X1,lim); title(typeLabels{2});
subplot(1,3,3); imshow(X2,lim); title(typeLabels{3});

% Save figure
%export_fig -pdf -transparent masks

%% Generate NRMSE curves

% L + S soft T
fftpath = '/home/brimoor/Research/otazo/fft/data.mat';
svtpath = '/home/brimoor/Research/otazo/svt_LpS_soft_T/data.mat';
optpath = '/home/brimoor/Research/otazo/opt_LpS_soft_T/data.mat';
pcapath = '/home/brimoor/Research/otazo/pca_LpS_soft_T/data.mat';
outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LpS_soft_T.mat';
%outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LpS_soft_T_full.mat';
GenerateNRMSEcurves(fftpath,svtpath,optpath,pcapath,outpath);

% L + S soft
fftpath = '/home/brimoor/Research/otazo/fft/data.mat';
svtpath = '/home/brimoor/Research/otazo/svt_LpS_soft/data.mat';
optpath = '/home/brimoor/Research/otazo/opt_LpS_soft/data.mat';
pcapath = '/home/brimoor/Research/otazo/pca_LpS_soft/data.mat';
outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LpS_soft.mat';
%outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LpS_soft_full.mat';
GenerateNRMSEcurves(fftpath,svtpath,optpath,pcapath,outpath);

% L + S mixed21
fftpath = '/home/brimoor/Research/otazo/fft/data.mat';
svtpath = '/home/brimoor/Research/otazo/svt_LpS_mixed21/data.mat';
optpath = '/home/brimoor/Research/otazo/opt_LpS_mixed21/data.mat';
pcapath = '/home/brimoor/Research/otazo/pca_LpS_mixed21/data.mat';
outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LpS_mixed21.mat';
%outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LpS_mixed21_full.mat';
GenerateNRMSEcurves(fftpath,svtpath,optpath,pcapath,outpath);

% L & S soft T
fftpath = '/home/brimoor/Research/otazo/fft/data.mat';
svtpath = '/home/brimoor/Research/otazo/svt_LaS_soft_T/data.mat';
optpath = '/home/brimoor/Research/otazo/opt_LaS_soft_T/data.mat';
pcapath = '/home/brimoor/Research/otazo/pca_LaS_soft_T/data.mat';
outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LaS_soft_T.mat';
%outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LaS_soft_T_full.mat';
GenerateNRMSEcurves(fftpath,svtpath,optpath,pcapath,outpath);

% L & S soft
fftpath = '/home/brimoor/Research/otazo/fft/data.mat';
svtpath = '/home/brimoor/Research/otazo/svt_LaS_soft/data.mat';
optpath = '/home/brimoor/Research/otazo/opt_LaS_soft/data.mat';
pcapath = '/home/brimoor/Research/otazo/pca_LaS_soft/data.mat';
outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LaS_soft.mat';
%outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LaS_soft_full.mat';
GenerateNRMSEcurves(fftpath,svtpath,optpath,pcapath,outpath);

% L & S mixed21
fftpath = '/home/brimoor/Research/otazo/fft/data.mat';
svtpath = '/home/brimoor/Research/otazo/svt_LaS_mixed21/data.mat';
optpath = '/home/brimoor/Research/otazo/opt_LaS_mixed21/data.mat';
pcapath = '/home/brimoor/Research/otazo/pca_LaS_mixed21/data.mat';
outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LaS_mixed21.mat';
%outpath = '/home/brimoor/Research/otazo/nrmse/NRMSE_LaS_mixed21_full.mat';
GenerateNRMSEcurves(fftpath,svtpath,optpath,pcapath,outpath);

%% Plot NRMSE curves

% Knobs
%SNR = 40:5:60; % SNR
SNR = 10:10:60; % SNR
type = 1:3; % {1 = Whole image, 2 = Heart-only, 3 = Body-only}
typeLabels = {'Whole image','Heart-only','Body-only'};
NRMSEpath = './nrmse2/NRMSE_LpS_soft_T.mat';
%NRMSEpath = './nrmse2/NRMSE_LpS_soft.mat';
%NRMSEpath = './nrmse2/NRMSE_LpS_mixed21.mat';
%NRMSEpath = './nrmse2/NRMSE_LaS_soft_T.mat';
%NRMSEpath = './nrmse2/NRMSE_LaS_soft.mat';
%NRMSEpath = './nrmse2/NRMSE_LaS_mixed21.mat';
pIdx = ':'; % {':',1:4}

% Load NRMSE data
load(NRMSEpath);

% Generate plots
figure;
Ntype = length(type);
Nsnr = length(SNR);
cm = linspecer(Nsnr);
[ni nj] = bestSubplotShape(Ntype);
snrIdx = 1;
algoIdx = round(0.5 * Ntype);
phndlAlgo = zeros(1,4);
pstrAlgo = {'FFT','SVT','OptShrink','PCA'};
phndlSNR = zeros(1,Nsnr);
pstrSNR = cell(1,Nsnr);
for i = 1:Ntype
    subplot(ni,nj,i);
    for j = 1:Nsnr
        % FFT
        [~,idx_SNR_fft] = min(abs(labels_fft.SNR - SNR(j)));
        nrmse_fft = squeeze(NRMSE_fft{type(i)}(:,idx_SNR_fft));
        phndlAlgo(1) = plot(labels_fft.p(pIdx),nrmse_fft(pIdx),'-','Color',cm(j,:));
        hold on;
        
        % SVT
        [~,idx_SNR_svt] = min(abs(labels_svt.SNR - SNR(j)));
        nrmse_svt = squeeze(NRMSE_svt{type(i)}(:,idx_SNR_svt));
        phndlAlgo(2) = plot(labels_svt.p(pIdx),nrmse_svt(pIdx),'-x','Color',cm(j,:));
        hold on;
        
        % OptShrink
        [~,idx_SNR_opt] = min(abs(labels_opt.SNR - SNR(j)));
        nrmse_opt = squeeze(NRMSE_opt{type(i)}(:,idx_SNR_opt));
        phndlAlgo(3) = plot(labels_opt.p(pIdx),nrmse_opt(pIdx),'-o','Color',cm(j,:));
        hold on;
        
        % PCA
        [~,idx_SNR_pca] = min(abs(labels_pca.SNR - SNR(j)));
        nrmse_pca = squeeze(NRMSE_pca{type(i)}(:,idx_SNR_pca));
        phndlAlgo(4) = plot(labels_pca.p(pIdx),nrmse_pca(pIdx),'-*','Color',cm(j,:));
        hold on;
        
        % SNR labels
        nrmse_opt2 = squeeze(NRMSE_opt{type(i)}(:,idx_SNR_opt));
        pstrSNR{j} = sprintf('SNR = %.0fdB',SNR(j));
        phndlSNR(j) = plot(labels_opt.p(pIdx),nrmse_opt2(pIdx),'-','Color',cm(j,:));
    end
    xlabel('p');
    ylabel('NRMSE');
    title(typeLabels{i});
    if (i == snrIdx)
        legend(phndlSNR,pstrSNR{:},'Location','NorthEast');
    end
    if (i == algoIdx)
        legend(phndlAlgo,pstrAlgo{:},'Location','NorthEast');
    end
    box on;
    axis tight;
    padAxis();
end
set(gcf,'Position',[432 527 1097 293]);

% Save figure
[~,pdfName,~] = fileparts(NRMSEpath);
export_fig('-pdf','-transparent',pdfName);

%% Robustness to rank parameter

% Knobs
%SNR = 40:5:60; % SNR
SNR = 10:10:60; % SNR
%p = 96 / 128; % [8 12 16 32 48 64 96 128] / 128;
p = 4 / 128; % [4 8 12 16 20 24 28 32] / 128;
type = 1; % {1 = Whole image, 2 = Heart-only, 3 = Body-only}
typeLabels = {'Whole image','Heart-only','Body-only'};
NRMSEpath = './nrmse2/NRMSE_LpS_soft_T_full.mat';
%NRMSEpath = './nrmse2/NRMSE_LpS_soft_full.mat';
%NRMSEpath = './nrmse2/NRMSE_LpS_mixed21_full.mat';
%NRMSEpath = './nrmse2/NRMSE_LaS_soft_T_full.mat';
%NRMSEpath = './nrmse2/NRMSE_LaS_soft_full.mat';
%NRMSEpath = './nrmse2/NRMSE_LaS_mixed21_full.mat';
%lambdaS = 0.1;
lambdaL_max = 1;

% constants
asymptotes = true;
prox = 'soft';

% Parse inputs
nSNR = length(SNR);
load(NRMSEpath);
regionStr = regexprep(lower(typeLabels{type}),'\W','_');
[~,algoStr,~] = fileparts(NRMSEpath);
inds = regexp(algoStr,'_');
algoStr = algoStr((inds(1) + 1):(inds(end) - 1));
pStr = regexprep(num2str(p),'\.','p');

% Get lambdaS values
if exist('lambdaS','var')
    % Fixed lambdaS value
    lambdaS_opt = repmat(lambdaS,1,nSNR);
    lambdaS_svt = repmat(lambdaS,1,nSNR);
    
    % Title strings
    lStr = 'fixed';
    titlestr_opt = sprintf('OptShrink [p = %.3f, \\lambda_S = %01.03f]',p,lambdaS);
    titlestr_svt = sprintf('SVT [p = %.3f, \\lambda_S = %01.03f]',p,lambdaS);
else
    % mNRMSE lambdaS values
    data = load(regexprep(NRMSEpath,'_full','')); % NRMSE data
    lambdaS_opt = zeros(1,nSNR);
    lambdaS_svt = zeros(1,nSNR);
    for i = 1:nSNR
        % OptShrink
        [~,p_mmse_opt_idx] = min(abs(data.labels_opt.p - p));
        [~,SNR_mmse_opt_idx] = min(abs(data.labels_opt.SNR - SNR(i)));
        inds_mmse_opt = data.idxopt_opt{type}{p_mmse_opt_idx,SNR_mmse_opt_idx};
        lambdaS_opt(i) = data.labels_opt.lambdaS(inds_mmse_opt(2));
        
        % SVT
        [~,p_mmse_svt_idx] = min(abs(data.labels_svt.p - p));
        [~,SNR_mmse_svt_idx] = min(abs(data.labels_svt.SNR - SNR(i)));
        inds_mmse_svt = data.idxopt_svt{type}{p_mmse_svt_idx,SNR_mmse_svt_idx};
        lambdaS_svt(i) = data.labels_svt.lambdaS(inds_mmse_svt(2));
    end
    
    % Title strings
    lStr = 'mmse';
    titlestr_opt = sprintf('OptShrink [p = %.3f, \\lambda_S = mmse]',p);
    titlestr_svt = sprintf('SVT [p = %.3f, \\lambda_S = mmse]',p);
end

% Initialize figure
cm = linspecer(nSNR);
ylim = [inf -inf];
%figure('Position',[750 500 560 230]);
figure('Position',[750 500 685 280]);

% Plot rank parameter robustness curves
phndl = zeros(1,nSNR);
pstr = cell(1,nSNR);
for i = 1:nSNR
    % Get OptShrink NRMSEs
    [~,p_opt_idx] = min(abs(labels_opt.p - p));
    [~,SNR_opt_idx] = min(abs(labels_opt.SNR - SNR(i)));
    [~,lambdaS_opt_idx] = min(abs(labels_opt.lambdaS - lambdaS_opt(i)));
    NRMSEopt = squeeze(NRMSE_opt{type}(p_opt_idx,SNR_opt_idx,1,:,lambdaS_opt_idx));
    r = labels_opt.r;
    
    % Get SVT NRMSEs
    [~,p_svt_idx] = min(abs(labels_svt.p - p));
    [~,SNR_svt_idx] = min(abs(labels_svt.SNR - SNR(i)));
    [~,lambdaS_svt_idx] = min(abs(labels_svt.lambdaS - lambdaS_svt(i)));
    NRMSEsvt = squeeze(NRMSE_svt{type}(p_svt_idx,SNR_svt_idx,1,:,lambdaS_svt_idx));
    lambdaL = labels_svt.lambdaL;
    
    % Clip to lambdaL range
    indsL = (lambdaL < lambdaL_max);
    NRMSEsvt = NRMSEsvt(indsL);
    lambdaL = lambdaL(indsL);
    
    % Compute common y-limits
    ylim = [min([ylim(1); NRMSEopt(:); NRMSEsvt(:)]) ...
            max([ylim(2); NRMSEopt(:); NRMSEsvt(:)])];
    
    % OptShrink
    subplot(1,2,1);
    phndl(i) = plot(r,NRMSEopt,'-o','Color',cm(i,:));
    pstr{i} = sprintf('SNR = %.0fdB',SNR(i));
    hold on;
    
    % SVT
    subplot(1,2,2);
    semilogx(lambdaL,NRMSEsvt,'-o','Color',cm(i,:));
    hold on;
end

% Finalize OptShrink subplot
subplot(1,2,1);
xlabel('r');
ylabel(sprintf('NRMSE(%s)',typeLabels{type}));
title(titlestr_opt);
legend(phndl,pstr{:},'Location','Northwest');
ax1 = gca;
axis(ax1,'tight');
ylim1 = get(ax1,'ylim');

% Finalize SVT subplot
subplot(1,2,2);
xlabel('\lambda_L');
ylabel(sprintf('NRMSE(%s)',typeLabels{type}));
title(titlestr_svt);
ax2 = gca;
axis(ax2,'tight');
ylim2 = get(ax2,'ylim');

% Make y-axes consistent
ylim = [min(ylim1(1),ylim2(1)), max(ylim1(2),ylim2(2))];
set(ax1,'ylim',ylim);
set(ax2,'ylim',ylim);

% Add theory asymptotes
if exist('asymptotes','var') && asymptotes
    % Get regularizer data
    cdata = load('./otazo_true.mat');
    Xtrue = reshape(cdata.Xtrue,[cdata.nd cdata.nt]);
    CL = svd(Xtrue); CL = CL(2);
    if strcmpi(prox,'soft')
        CS = max(abs(Xtrue(:))); % soft
    elseif strcmpi(prox,'mixed21')
        CS = max(sqrt(sum(abs(Xtrue).^2,2))); % mixed21
    end
    n = cdata.nd;
    
    % Plot low-rank theory lines
    subplot(1,2,2);
    ylim = get(gca,'ylim');
    for i = 1:nSNR
        %i = round(nSNR / 2);
        lambdaL_theory = sqrt(n) * CS * lambdaS_svt(i) / CL;
        %plot(lambdaL_theory * [1 1],ylim,'k--');
        plot(lambdaL_theory * [1 1],ylim,'--','Color',cm(i,:));
    end
end

% Save figure
drawnow; pause(0.1);
figStr = sprintf('sensitivity_%s_%s_%s_%s',algoStr,regionStr,lStr,pStr);
export_fig('-pdf','-transparent',figStr);

%% Run algorithms

% Knobs
p = 16 / 128; % [8 12 16 32 48 64 96 128] / 128;
SNR = 60; % 40:5:60
seed = 1; % 1
algoStr = 'L+S PGM'; % {'L+S PGM','L&S ADMM'}
proxS = 'soft'; % {'soft','mixed21'}
T = TempFFT(2); % {1,TempFFT(2)}
r = 1; % 1:10
lambdaL = 0.1; % logspace(log10(0.001),log10(10),11)
lambdaS = 0.1; % logspace(log10(0.001),log10(10),11)
inpath = './otazo_true.mat';
proxpath = './prox/';

% Perform reconstructions
addpath(regexprep(proxpath,'\','/'));
[NRMSE_fft, recon_fft, Xtrue, ROIs] = run_fft_algo(p,SNR,seed,inpath);
[NRMSE_svt, recon_svt,     ~,    ~] = run_svt_algo(p,SNR,seed,inpath,lambdaL,proxS,lambdaS,algoStr,T);
[NRMSE_pca, recon_pca,     ~,    ~] = run_pca_algo(p,SNR,seed,inpath,r,proxS,lambdaS,algoStr,T);
[NRMSE_opt, recon_opt,     ~,    ~] = run_optshrink_algo(p,SNR,seed,inpath,r,proxS,lambdaS,algoStr,T);
[ny nx] = size(ROIs{1});
nt = size(Xtrue,2);
Xtrue = reshape(Xtrue,[ny nx nt]);
Xhat_fft = reshape(recon_fft.X,[ny nx nt]);
Xhat_svt = reshape(recon_svt.X,[ny nx nt]);
Xhat_pca = reshape(recon_pca.X,[ny nx nt]);
Xhat_opt = reshape(recon_opt.X,[ny nx nt]);

% Play movie
label = @(str,NRMSE) sprintf('%s (%.4f)',str,NRMSE.X(end));
movie.video = cat(1,cat(2,Xtrue,Xhat_fft,Xhat_svt), ...
                    cat(2,Xtrue,Xhat_pca,Xhat_opt));
opts.xlabels = {'truth',label('fft',NRMSE_fft),label('svt',NRMSE_svt);
                'truth',label('pca',NRMSE_pca),label('opt',NRMSE_opt)};
PlayMovie(movie,opts);
