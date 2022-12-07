%% [SPOT] drpca

% Knobs
inpath   = 'otazo_R8.mat';
initpath = 'otazo_R8_lps_visual.mat';
r        = nan;
lambdaL  = 0.5;
lambdaS  = 0.02;
lambdaB  = 0.03;
lambdaB0 = 0.03;
dr       = 1;
nIters   = 50;
np       = 0;
fixedD   = false;
%{
inpath   = 'invivo_R8.mat';
initpath = 'invivo_R8_lps_nrmse.mat';
r        = nan;
lambdaL  = 0.05;
lambdaS  = 0.01;
lambdaB  = 0.01;
lambdaB0 = 0.01;
dr       = 5;
nIters   = 50;
np       = 0;
fixedD   = false;
%}

% Load dependencies
addpath('./deps_lassi');

% Load undersampled data (Y, mask, Xfft, Xtrue)
data = load(inpath);
[ny, nx, nt] = size(data.Xtrue);
A  = Emat_xyt(data.mask,data.samp,[ny, nx, nt]);
%A = Afft(data.mask,[ny, nx, nt]);

% Load initialization
init = load(initpath);

% DRPCA
opts.A        = A;
opts.r        = r;
opts.sdim     = [ny, nx, nt];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 2];
opts.np       = np;
opts.type     = 'hard';
opts.dr       = dr;
opts.ddim     = [64, 5];
opts.fixedD   = fixedD;
opts.nIters   = nIters;
opts.nItersDB = 1;
opts.nItersLS = 5;
opts.L0       = reshape(init.Lhat,[],nt);
opts.S0       = reshape(init.Shat,[],nt);
opts.D0       = dctmtx(prod(opts.pdim));
opts.Xtrue    = reshape(data.Xtrue,[],nt);
opts.accel    = false;
opts.tau      = 1;
opts.flag     = 1;
lB = logspace(log10(lambdaB0),log10(lambdaB),nIters);
[Lhat, Shat, ~, ~,stats] = drpca(Y,lambdaL,lambdaS,lB,opts);
Lhat = reshape(Lhat,ny,nx,nt);
Shat = reshape(Shat,ny,nx,nt);

% Play movie
optsPM.xlabels = {'Xtrue','Xfft','Lhat + Shat','Lhat','Shat'};
PlayMovie(cat(2,data.Xtrue,data.Xfft,Lhat + Shat,Lhat,Shat),optsPM);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([950, 250]);

% Cost
subplot(1,3,1);
plot(1:nIters,stats.cost,'b-o');
xlabel('Iteration');
title('Cost');
axis tight; padAxis();

% NRMSE
subplot(1,3,2);
plot(1:nIters,stats.nrmse,'b-o');
xlabel('Iteration');
title('NRMSE');
axis tight; padAxis();

% Sparsity
subplot(1,3,3);
plot(1:nIters,stats.sparsity,'b-o');
xlabel('Iteration');
title('Sparsity');
axis tight; padAxis();
%--------------------------------------------------------------------------

%% [PLOT] LASSI Dfixed

% Knobs
%inpath = 'LASSI_R8_Dfixed_par1.mat';
inpath  = 'invivo_LASSI_R8_Dfixed_par1.mat';
SAVE_FIGURE = true;

% Load data (Lhat, Shat, vars, stats)
load(inpath);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([780, 424]);
ax  = [];
str = {'Learned','DCT'};

cax = zeros(1,2);
nax = zeros(1,2);
sax = zeros(1,2);
for i = 1:2
    % Cost trajectory
    cax(i) = subplot(2,3,1 + 3 * (i - 1));
    plot(1:vars.nIters,stats.cost(i,:),'b-');
    xlabel('Iteration');
    ylabel('Cost');
    title(sprintf('%s: cost',str{1 + vars.fixedD(i)}));
    axis tight; padAxis();
    
    % NRMSE trajectory
    nax(i) = subplot(2,3,2 + 3 * (i - 1));
    plot(1:vars.nIters,stats.nrmse(i,:),'b-');
    xlabel('Iteration');
    ylabel('NRMSE');
    title(sprintf('%s: NRMSE (%.2f%%)',str{1 + vars.fixedD(i)},100 * stats.nrmse(i,end)));
    axis tight; padAxis();
    
    % Sparsity trajectory
    sax(i) = subplot(2,3,3 + 3 * (i - 1));
    plot(1:vars.nIters,stats.sparsity(i,:),'b-');
    xlabel('Iteration');
    ylabel('Sparsity %');
    title(sprintf('%s: sparsity',str{1 + vars.fixedD(i)}));
    axis tight; padAxis();
end

% Match axes
matchAxes(cax,[],'y');
matchAxes(nax,[],'y');
matchAxes(sax,[],'y');

% Save figure
if exist('SAVE_FIGURE','var') && ~isempty(SAVE_FIGURE)
    [~, name] = fileparts(inpath);
    export_fig('-pdf','-transparent',[name, '.pdf']);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [PLOT] LASSI Diter

% Knobs
%inpath = 'LASSI_R8_Diter_par1.mat';
inpath  = 'invivo_LASSI_R8_Diter_par1.mat';
SAVE_FIGURE = true;

% Load data (Lhat, Shat, vars, stats)
load(inpath);
[m, n] = size(stats.nrmse);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([416, 344]);
cm = linspecer(m);

phndl = zeros(1,m);
pstr  = cell(1,m);
for i = 1:m
    phndl(i) = plot(1:n,stats.nrmse(i,:),'o','Color',cm(i,:)); hold on;
    pstr{i}  = sprintf('Outer: %02d,  iDB: %02d,  iLS: %02d', ...
                        vars.nIters(i),vars.nItersDB(i),vars.nItersLS(i));
end
xlabel('(D,B) updates');
ylabel('NRMSE');
title('NRMSE');
legend(phndl,pstr{:});
axis tight; padAxis();

% Save figure
if exist('SAVE_FIGURE','var') && ~isempty(SAVE_FIGURE)
    [~, name] = fileparts(inpath);
    export_fig('-pdf','-transparent',[name, '.pdf']);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [PLOT] LASSI Dover

% Knobs
%inpath = 'LASSI_R8_Dover_par1.mat';
%inpath = 'LASSI_R8_Dover_par2.mat';
%inpath = 'LASSI_R8_Dover_par3.mat';
%inpath = 'invivo_LASSI_R8_Dover_par1.mat';
%inpath = 'invivo_LASSI_R8_Dover_par2.mat';
inpath  = 'invivo_LASSI_R8_Dover_par3.mat';
SAVE_FIGURE = true;

% Load data (Lhat, Shat, vars, stats)
load(inpath);
m = numel(vars.np);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([416, 344]);
cm = linspecer(m);

phndl = zeros(1,m);
pstr  = cell(1,m);
for i = 1:m
    phndl(i) = plot(1:vars.nIters,stats.nrmse(i,:),'-','Color',cm(i,:)); hold on;
    pstr{i}  = sprintf('D = [DCT, %d patches]',vars.np(i));
end
xlabel('Iteration');
ylabel('NRMSE');
title('NRMSE');
legend(phndl,pstr{:});
axis tight; padAxis();

% Save figure
if exist('SAVE_FIGURE','var') && ~isempty(SAVE_FIGURE)
    [~, name] = fileparts(inpath);
    export_fig('-pdf','-transparent',[name, '.pdf']);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [PLOT] LASSI Dfixed

% Knobs
%inpath = 'LASSI_R8_LSinit_par1.mat';
inpath  = 'invivo_LASSI_R8_LSinit_par1.mat';
SAVE_FIGURE = true;

% Load data (Lhat, Shat, vars, stats)
load(inpath);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([780, 424]);
ax  = [];
str = {'L+S','k-t SLR'};

cax = zeros(1,2);
nax = zeros(1,2);
sax = zeros(1,2);
for i = 1:2
    % Cost trajectory
    cax(i) = subplot(2,3,1 + 3 * (i - 1));
    plot(1:vars.nIters,stats.cost(i,:),'b-');
    xlabel('Iteration');
    ylabel('Cost');
    title(sprintf('%s: cost',str{i}));
    axis tight; padAxis();
    
    % NRMSE trajectory
    nax(i) = subplot(2,3,2 + 3 * (i - 1));
    plot(1:vars.nIters,stats.nrmse(i,:),'b-');
    xlabel('Iteration');
    ylabel('NRMSE');
    title(sprintf('%s: NRMSE (%.2f%%)',str{i},100 * stats.nrmse(i,end)));
    axis tight; padAxis();
    
    % Sparsity trajectory
    sax(i) = subplot(2,3,3 + 3 * (i - 1));
    plot(1:vars.nIters,stats.sparsity(i,:),'b-');
    xlabel('Iteration');
    ylabel('Sparsity %');
    title(sprintf('%s: sparsity',str{i}));
    axis tight; padAxis();
end

% Match axes
matchAxes(cax,[],'y');
matchAxes(nax,[],'y');
matchAxes(sax,[],'y');

% Save figure
if exist('SAVE_FIGURE','var') && ~isempty(SAVE_FIGURE)
    [~, name] = fileparts(inpath);
    export_fig('-pdf','-transparent',[name, '.pdf']);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [PLOT] LASSI Sonly - tuning

% Knobs
%inpath = 'LASSI_R8_Sonly_par2.mat';
%inpath = 'LASSI_R8_Sonly_par2b.mat';
%inpath = 'LASSI_R8_Sonly_par2c.mat';
%inpath = 'LASSI_R8_Sonly_par3.mat';
%inpath = 'LASSI_R8_Sonly_par3b.mat';
%inpath = 'LASSI_R8_Sonly_par4.mat';
%inpath = 'LASSI_R8_Sonly_par5.mat';
%inpath = 'LASSI_R8_Sonly_par6.mat';
%inpath = 'LASSI_R8_Sonly_par6b.mat'; % incomplete, but sufficient
SAVE_FIGURE = true;

% Load data (Lhat, Shat, vars, stats)
load(inpath);

% Compute optimal parameters
nrmse = stats.nrmse(:,:,:,:,:,end);
[mnrmse, idx] = min(nrmse(:));
[ir, il, is, ib, id] = ind2sub(size(nrmse),idx);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([839, 463]);
ax = [];

% lambdaS
subplot(2,3,1);
if (numel(vars.lambdaS) > 1), ax(end + 1) = gca(); end
semilogx(vars.lambdaS,squeeze(nrmse(ir,il,:,ib,id)),'b-o');
xlabel('\lambda_S');
ylabel('NRMSE');
title(sprintf('\\lambda_S = %.3f',vars.lambdaS(is)));
axis tight; padAxis();

% lambdaB
subplot(2,3,2);
if (numel(vars.lambdaB) > 1), ax(end + 1) = gca(); end
semilogx(vars.lambdaB,squeeze(nrmse(ir,il,is,:,id)),'b-o');
xlabel('\lambda_B');
ylabel('NRMSE');
title(sprintf('\\lambda_B = %.3f',vars.lambdaB(ib)));
axis tight; padAxis();

% dr
subplot(2,3,3);
if (numel(vars.dr) > 1), ax(end + 1) = gca(); end
semilogx(vars.dr,squeeze(nrmse(ir,il,is,ib,:)),'b-o');
xlabel('r_d');
ylabel('NRMSE');
title(sprintf('r_d = %d',vars.dr(id)));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,3,4);
plot(1:vars.nIters,squeeze(stats.cost(ir,il,is,ib,id,:)),'b-');
xlabel('Iteration');
ylabel('cost');
title('Optimal cost');
axis tight; padAxis();

% Optimal MSE trajectory
subplot(2,3,5);
plot(1:vars.nIters,squeeze(stats.nrmse(ir,il,is,ib,id,:)),'b-');
xlabel('Iteration');
ylabel('NRMSE');
title(sprintf('Optimal NRMSE (%.2f%%)',100 * mnrmse));
axis tight; padAxis();

% Optimal sparsity trajectory
subplot(2,3,6);
plot(1:vars.nIters,squeeze(stats.sparsity(ir,il,is,ib,id,:)),'b-');
xlabel('Iteration');
ylabel('Sparsity %');
title('Optimal sparsity');
axis tight; padAxis();

% Match axes
if ~isempty(ax), matchAxes(ax,[],'y'); end

% Save figure
if SAVE_FIGURE
    [~, name] = fileparts(inpath);
    export_fig('-pdf','-transparent',[name, '.pdf']);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [PLOT] LASSI Sonly - dr

% Knobs
%inpath = 'LASSI_R8_Sonly_par1.mat';
%inpath = 'LASSI_R8_Sonly_par1b.mat';
%inpath = 'LASSI_R8_Sonly_par1c.mat';
%inpath = 'LASSI_R8_Sonly_par6c.mat';
inpath = 'LASSI_R8_Sonly_par6d.mat';
%inpath = 'invivo_LASSI_R8_Sonly_par1.mat';
SAVE_FIGURE = true;

% Load data (Lhat, Shat, vars, stats)
load(inpath);

% Exract data for the parameters of interest
nrmse = squeeze(stats.nrmse(:,:,:,:,:,end));
gain  = 20 * log10(nrmse(end) ./ nrmse);
m     = numel(nrmse);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([416, 344]);
cm = linspecer(m);

bar(vars.dr,gain);
xlabel('r_d');
ylabel('Gain (dB)');
title('NRMSE gain (\lambda_L = inf)');
axis tight; padAxis();

% Save figure
if exist('SAVE_FIGURE','var') && ~isempty(SAVE_FIGURE)
    [~, name] = fileparts(inpath);
    export_fig('-pdf','-transparent',[name, '.pdf']);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [PLOT] LASSI Sonly - sweep

%{
% Knobs: [8, 8, 5] patches
data(1).inpath = './RPCA_full_par1b.mat';
data(2).inpath = './RPCA_full_par1.mat';
data(3).inpath = './LASSI_full_Sonly_par1.mat';
data(4).inpath = './LASSI_full_Sonly_par2.mat';
data(5).inpath = './LASSI_full_Sonly_par3.mat';
data(6).inpath = './LASSI_full_Sonly_par4.mat';
data(7).inpath = './LASSI_full_Sonly_par5b.mat';
data(8).inpath = './LASSI_full_Sonly_par6.mat';
data(1).label = 'L + S [mmse]';
data(2).label = 'L + S [visual]';
data(3).label = 'L + S [visual] + SOUP-MRI [r = 1]';
data(4).label = 'L + S [mmse] + SOUP-MRI [r = 1]';
data(5).label = 'L + S [visual] + SOUP-MRI [r = 5]';
data(6).label = 'L + S [mmse] + SOUP-MRI [r = 5]';
data(7).label = 'L + S [visual] + SOUP-DCT';
data(8).label = 'L + S [mmse] + SOUP-DCT';
outpath = 'sweep_otazo885.pdf';
sortIdx = 2;
%}

% Knobs: [8, 8, 8] patches
data(1).inpath = './RPCA_full_par1b.mat';
data(2).inpath = './RPCA_full_par1.mat';
data(3).inpath = './LASSI_full_Sonly_par7a.mat';
data(4).inpath = './LASSI_full_Sonly_par7b.mat';
data(5).inpath = './LASSI_full_Sonly_par7c.mat';
data(6).inpath = './LASSI_full_Sonly_par7d.mat';
data(7).inpath = './LASSI_full_Sonly_par7e.mat';
data(8).inpath = './LASSI_full_Sonly_par7f.mat';
data(1).label = 'L + S [mmse]';
data(2).label = 'L + S [visual]';
data(3).label = 'L + S [visual] + SOUP-MRI [r = 1]';
data(4).label = 'L + S [mmse] + SOUP-MRI [r = 1]';
data(5).label = 'L + S [visual] + SOUP-MRI [r = 8]';
data(6).label = 'L + S [mmse] + SOUP-MRI [r = 8]';
data(7).label = 'L + S [visual] + SOUP-DCT';
data(8).label = 'L + S [mmse] + SOUP-DCT';
outpath = 'sweep_otazo888.pdf';
sortIdx = 2;

% Load data
nData = numel(data);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results
    data(i).ps    = datai.vars.ps;
    data(i).nrmse = nanmean(datai.stats.nrmse(:,:,end),2);
    
    % Compute gains
    data(i).gain  = 20 * log10(data(1).nrmse ./ data(i).nrmse);
    fprintf('Gain %s:\n',data(i).label);
    fprintf('%.2f, ',data(i).gain);
    fprintf('\n');
end

% Sort data, if requested
if exist('sortIdx','var') && ~isempty(sortIdx)
    data = data(argsort(cellfun(@(n)n(sortIdx),{data.nrmse}),'descend'));
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([437, 358]);
cm = flipud(linspecer(nData));

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    phndl(i) = semilogx(data(i).ps,data(i).nrmse,'o-','Color',cm(i,:)); hold on;
end
xlabel('ps');
ylabel('NRMSE');
title('otazo');
legend(phndl,data.label);
axis tight; padAxis();

% Save figure, if requested
if exist('outpath','var') && ~isempty(outpath)
    export_fig('-pdf','-transparent',outpath);
end
%--------------------------------------------------------------------------

%% [TABLE] LASSI Sonly - sweep

%{
% Knobs: [8, 8, 5] patches, visual
data(1).inpath = './RPCA_full_par1b.mat';
data(2).inpath = './LASSI_full_Sonly_par5b.mat'; % visual
data(3).inpath = './LASSI_full_Sonly_par3.mat';  % visual
data(4).inpath = './LASSI_full_Sonly_par1.mat';  % visual
data(1).label = 'L + S';
data(2).label = 'Fixed D';
data(3).label = 'r = 5';
data(4).label = 'r = 1';
tList = [4];
%}

% Knobs: [8, 8, 5] patches, mmse
data(1).inpath = './RPCA_full_par1b.mat';
data(2).inpath = './LASSI_full_Sonly_par6.mat'; % mmse
data(3).inpath = './LASSI_full_Sonly_par4.mat'; % mmse
data(4).inpath = './LASSI_full_Sonly_par2.mat'; % mmse
data(1).label = 'L + S';
data(2).label = 'Fixed D';
data(3).label = 'r = 5';
data(4).label = 'r = 1';
tList = [4];

%{
% Knobs: [8, 8, 8] patches, visual
data(1).inpath = './RPCA_full_par1b.mat';
data(2).inpath = './LASSI_full_Sonly_par7e.mat'; % visual
data(3).inpath = './LASSI_full_Sonly_par7c.mat'; % visual
data(4).inpath = './LASSI_full_Sonly_par7a.mat'; % visual
data(1).label = 'L + S';
data(2).label = 'Fixed D';
data(3).label = 'r = 8';
data(4).label = 'r = 1';
tList = [4];
%}

%{
% Knobs: [8, 8, 8] patches, mmse
data(1).inpath = './RPCA_full_par1b.mat';
data(2).inpath = './LASSI_full_Sonly_par7f.mat'; % mmse
data(3).inpath = './LASSI_full_Sonly_par7d.mat'; % mmse
data(4).inpath = './LASSI_full_Sonly_par7b.mat'; % mmse
data(1).label = 'L + S';
data(2).label = 'Fixed D';
data(3).label = 'r = 8';
data(4).label = 'r = 1';
tList = [4];
%}

% Compute gains
nData = numel(data);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results
    data(i).ps    = datai.vars.ps;
    data(i).nrmse = 100 * nanmean(datai.stats.nrmse(:,tList,end),2)';
    
    % Compute gains
    data(i).gain  = 20 * log10(data(1).nrmse ./ data(i).nrmse);
end
gainlr = 20 * log10(data(3).nrmse ./ data(4).nrmse);

% Display results
accel = 1 ./ data(1).ps;
line  = @() disp(repmat('-',1,75));
line();
fprintf('%25s ','Acceleration'); fprintf(' %2dx    ',accel); fprintf('\n');
line();
for i = 1:nData
    fprintf('%25s ',data(i).label); fprintf(' %.2f%% ',data(i).nrmse); fprintf('\n');
    if i > 1
        fprintf('%25s ','Gain over L + S (dB)'); fprintf(' % .2f  ',data(i).gain); fprintf('\n');
    end
    if i < nData
        line();
    end
end
fprintf('%25s ','Gain over r = 5 (dB)'); fprintf(' % .2f  ',gainlr); fprintf('\n');

%% [IMAGES] LASSI Sonly - R8 otazo

% Knobs
truthpath = 'otazo_R8.mat';                     % Ground truth
%souppath = 'recons/XXXXXXXXXXXXXXXXXXXX.mat';  % visual, DCT
%souppath = 'recons/XXXXXXXXXXXXXXXXXXXX.mat';  % mmse,   DCT
%souppath = 'LASSI_R8_Sonly_par2.mat';          % visual, r = 1
souppath  = 'LASSI_R8_Sonly_par2b.mat';         % mmse,   r = 1
%souppath = 'LASSI_R8_Sonly_par2c.mat';         % visual, r = 5
%souppath = 'recons/XXXXXXXXXXXXXXXXXXXX.mat';  % mmse,   r = 5
%lpspath  = 'otazo_R8_lps_visual.mat';          % visual
lpspath   = 'otazo_R8_lps_mse.mat';             % mmse

% Load ground truth (Xtrue)
load(truthpath);
NRMSEfcn = @(Xhat,Xref) 100 * norm(Xhat(:) - Xref(:)) / norm(Xref(:));
nt       = size(Xtrue,3);

% Load SOUP (Shat, Dhat)
load(souppath);
XSOUP     = reshape(Shat,[128, 128, 40]);
NRMSESOUP = NRMSEfcn(XSOUP,Xtrue) %#ok

% Load L + S
load(lpspath);
XLpS     = Lhat + Shat;
NRMSELpS = NRMSEfcn(XLpS,Xtrue) %#ok

% Movie
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.xlabels = {'Reference','SOUP-DILLA [r = 1]','L + S'};
PlayMovie(cat(2,Xtrue,XSOUP,XLpS),opts);
%}

% Per-frame errors
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Per-frame errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute errors
errLpS  = zeros(1,nt);
errSOUP = zeros(1,nt);
for i = 1:nt
    errLpS(i)  = NRMSEfcn(XLpS(:,:,i) ,Xtrue(:,:,i));
    errSOUP(i) = NRMSEfcn(XSOUP(:,:,i),Xtrue(:,:,i));
end

% Compute gains
gainDB   = 20 * log10(errLpS ./ errSOUP);
[~, idx] = sort(gainDB,'descend');
largestGains = idx(1:5) %#ok

% Plot error curves
cfigure([711, 263]);
cm = linspecer(3);

subplot(1,2,1);
plot(1:nt,errLpS ,'o-','Color',cm(1,:)); hold on;
plot(1:nt,errSOUP,'o-','Color',cm(2,:));
xlabel('Frame');
ylabel('NRMSE');
title('Per-frame error');
axis tight; padAxis();

subplot(1,2,2);
plot(1:nt,gainDB ,'o-','Color',cm(3,:)); hold on;
xlabel('Frame');
ylabel('Gain (dB)');
title('Per-frame gain');
axis tight; padAxis();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

% Error movie
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Knobs
f1    = 11; % [7, 13]
alpha = 1.25; % Alpha > 1 to clip large values
gamma = 1.25; % Gamma > 1 to emphasize large values

% Compute error maps
err1 = abs(abs(Xtrue) - abs(XLpS));
err2 = abs(abs(Xtrue) - abs(XSOUP));

% Gamma correction
err1 = sign(err1) .* abs(err1).^gamma;
err2 = sign(err2) .* abs(err2).^gamma;

% Scale data
L    = max(max(abs(err1(:))),max(abs(err2(:))));
err1 = min(1,alpha * err1 / L);
err2 = min(1,alpha * err2 / L);

% Error movie
E = cat(2,err1,err2);
opts.xlabels = {'L + S','SOUP-DILLA'};
opts.mag = 2;
PlayMovie(E,opts);
%}

% Error maps
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Knobs
f1    = 11; % [7, 13]
gamma = 1; % Gamma > 1 to emphasize large values
%cm   = spectral(64);
%cm   = jet(64);

% Compute error maps
err1 = abs(abs(Xtrue(:,:,f1)) - abs(XLpS(:,:,f1)));
err2 = abs(abs(Xtrue(:,:,f1)) - abs(XSOUP(:,:,f1)));

% Gamma correction
err1 = sign(err1) .* abs(err1).^gamma;
err2 = sign(err2) .* abs(err2).^gamma;

% Compute axis limits
L     = max(max(abs(err1(:))),max(abs(err2(:))));
%clim = [-L, L];
clim  = [0, L];

% Plot error maps
%cfigure([302, 161]);
figure;

subplot(1,2,1);
pcolor(flipud(err1));
axis equal; axis tight; axis off;
shading flat;
%colormap(cm);
caxis(clim);
colorbar;
title('L + S');

subplot(1,2,2);
pcolor(flipud(err2));
axis equal; axis tight; axis off;
shading flat;
%colormap(cm);
caxis(clim);
colorbar;
title('SOUP-DILLA');
%}

% Reconstruction frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Knobs
f1    = 7;
f2    = 13;
alpha = 1;
gamma = 0.75;

% Gamma correction
Xt = sign(Xtrue) .* abs(Xtrue).^gamma;
Xs = sign(XSOUP) .* abs(XSOUP).^gamma;
Xl = sign(XLpS)  .* abs(XLpS).^gamma;

% Clipping
L = max([Xt(:); Xs(:); Xl(:)]);
Xt = min(alpha * Xt / L,1);
Xs = min(alpha * Xs / L,1);
Xl = min(alpha * Xl / L,1);

%{
C = {Xt(:,:,f1), Xs(:,:,f1), Xl(:,:,f1);
     Xt(:,:,f2), Xs(:,:,f2), Xl(:,:,f2)};
%}
C = {Xt(:,:,f1), Xs(:,:,f1);
     Xt(:,:,f2), Xs(:,:,f2)};
I = cell2mov(C,1,1);

% Display image
%opts.xlabels = {'Reference','DINO-KAT MRI','L + S'};
opts.xlabels = {'Reference','DINO-KAT MRI'};
opts.ylabels = {sprintf('F%d',f2), sprintf('F%d',f1)};
opts.mag     = 1;
opts.repeat  = false;
PlayMovie(I,opts);

%{
% Save image
export_fig -pdf -transparent otazo_dinokat
%}

% ROIs
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Knobs
f1 = 13;
r  = 33:96; % Rows
c  = 33:96; % Columns

% Extract frames
L  = max(abs(Xtrue(:)));
Xt = min(abs(Xtrue) / L,1);
Xs = min(abs(XSOUP) / L,1);
Xl = min(abs(XLpS) / L ,1);
C  = {Xt(r,c,f1), Xs(r,c,f1), Xl(r,c,f1)};
I  = cell2mov(C,1,1);

% Display image
opts.xlabels = {'Reference','SOUP-DILLA','L + S'};
opts.ylabels = {sprintf('F%d',f1)};
opts.mag     = 1;
opts.repeat  = false;
PlayMovie(I,opts);
%}

%% [MOVIE] recons (by dataset)

flipFcn    = @(X) flipdim(flipdim(X,1),2);
flipInvivo = @(C) cellfun(flipFcn,C,'UniformOutput',false);

% Knobs (otazo 8x)
truthpath = 'otazo_full.mat';
dinopath  = 'data_TMI/recon_otazo_full_par2.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par3.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par4.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 0.90;  % Gamma < 1 to increase contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 1:128;
cols     = 1:128;
outpath  = './globalsip_otazo_8x_recon.gif';

%{
% Knobs (otazo 20x)
truthpath = 'otazo_full.mat';
dinopath  = 'data_TMI/recon_otazo_full_par10.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par11.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par12.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 0.90;  % Gamma < 1 to increase contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 1:128;
cols     = 1:128;
outpath  = './globalsip_otazo_20x_recon.gif';
%}

%{
% Knobs (invivo 8x)
truthpath = 'invivo_full.mat';
dinopath  = 'data_TMI/recon_invivo_full_par2.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par3.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par4.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 0.75;  % Gamma < 1 to increase contrast
mag      = 1.5;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 21:170;
cols     = 1:90;
outpath  = './globalsip_invivo_8x_recon.gif';
%}

%{
% Knobs (invivo 5x)
truthpath = 'invivo_full.mat';
dinopath  = 'data_TMI/recon_invivo_full_par7.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par8.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par9.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 0.75;  % Gamma < 1 to increase contrast
mag      = 1.5;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 21:170;
cols     = 1:90;
outpath  = './globalsip_invivo_5x_recon.gif';
%}

%{
% Knobs (pincat 9x)
truthpath = 'pincat_full.mat';
dinopath  = 'data_TMI/recon_pincat_full_par2.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par3.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par4.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 1.00;  % Gamma < 1 to increase contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 1:128;
cols     = 1:128;
outpath  = './globalsip_pincat_9x_recon.gif';
%}

%{
% Knobs (pincat 14x)
truthpath = 'pincat_full.mat';
dinopath  = 'data_TMI/recon_pincat_full_par6.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par7.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par8.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 1.00;  % Gamma < 1 to increase contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 1:128;
cols     = 1:128;
outpath  = './globalsip_pincat_14x_recon.gif';
%}

% Load data
truthdata = load(truthpath,'Xtrue');
dinodata  = load(dinopath,'Lhat','Shat');
rpcapath  = load(rpcapath,'L0','S0');
ktslrpath = load(ktslrpath,'X0');
Xref = truthdata.Xtrue;
Xdin = reshape(dinodata.Lhat + dinodata.Shat,size(Xref));
Xlps = reshape(rpcapath.L0 + rpcapath.S0,size(Xref));
Xkts = reshape(ktslrpath.X0,size(Xref));

% Check NRMSEs
NRMSEdin = 100 * norm(Xdin(:) - Xref(:)) / norm(Xref(:));
NRMSElps = 100 * norm(Xlps(:) - Xref(:)) / norm(Xref(:));
NRMSEkts = 100 * norm(Xkts(:) - Xref(:)) / norm(Xref(:));

% Correct images
LX    = max(abs(Xref(:)));
Xrefd = min(max(alpha * Xref / LX,0),1).^gamma;
Xdin2d = min(max(alpha * Xdin / LX,0),1).^gamma;
Xlpsd = min(max(alpha * Xlps / LX,0),1).^gamma;
Xktsd = min(max(alpha * Xkts / LX,0),1).^gamma;

% Clip images
Xrefd = Xrefd(rows,cols,:);
Xdin2d = Xdin2d(rows,cols,:);
Xlpsd = Xlpsd(rows,cols,:);
Xktsd = Xktsd(rows,cols,:);

% Generate movie
C = {Xrefd,Xdin2d,Xlpsd,Xktsd};
if ~isempty(strfind(truthpath,'invivo'))
    % Flip invivo images
    C = flipInvivo(C);
end
M = cell2mov(C,1,1);

% Play movie
lbl = @(s,n) sprintf('%s (%.1f%%)',s,n);
movie.video = M;
movie.Fv    = Fv;
opts.xlabels = {'Reference', ...
                lbl('DINO-KAT',NRMSEdin), ...
                lbl('L+S',NRMSElps), ...
                lbl('k-t SLR',NRMSEkts)};
opts.ylabels = {''};
opts.mag      = mag;
opts.fontSize = fontSize;
opts.gap      = gap;
P = PlayMovie(movie,opts);
P.SaveGIF(outpath);
P.Close();

%% [MOVIE] recons + error maps (by dataset)

flipFcn    = @(X) flipdim(flipdim(X,1),2);
flipInvivo = @(C) cellfun(flipFcn,C,'UniformOutput',false);

% Knobs (otazo 8x) [YES]
truthpath = 'otazo_full.mat';
dinopath  = 'data_TMI/recon_otazo_full_par2.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par3.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par4.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.90;  % Gamma < 1 to decrease contrast
alphaE = 0.40;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 1:128;
cols = 1:128;
outpath = './globalsip_errors_otazo_8x.gif';

%{
% Knobs (otazo 20x) [NO]
truthpath = 'otazo_full.mat';
dinopath  = 'data_TMI/recon_otazo_full_par10.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par11.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par12.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.90;  % Gamma < 1 to decrease contrast
alphaE = 0.75;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 1:128;
cols = 1:128;
outpath = './globalsip_errors_otazo_20x.gif';
%}

%{
% Knobs (invivo 8x) [YES]
truthpath = 'invivo_full.mat';
dinopath  = 'data_TMI/recon_invivo_full_par2.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par3.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par4.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.75;  % Gamma < 1 to decrease contrast
alphaE = 0.48;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1.5;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 31:160;
cols = 1:90;
outpath = './globalsip_errors_invivo_8x.gif';
%}

%{
% Knobs (invivo 5x) [NO]
truthpath = 'invivo_full.mat';
dinopath  = 'data_TMI/recon_invivo_full_par7.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par8.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par9.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.75;  % Gamma < 1 to decrease contrast
alphaE = 0.50;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1.5;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 31:160;
cols = 1:90;
outpath = './globalsip_errors_invivo_5x.gif';
%}

%{
% Knobs (pincat 9x) [YES]
truthpath = 'pincat_full.mat';
dinopath  = 'data_TMI/recon_pincat_full_par2.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par3.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par4.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 1.00;  % Gamma < 1 to decrease contrast
alphaE = 0.80;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 1:128;
cols = 1:128;
outpath = './globalsip_errors_pincat_9x.gif';
%}

%{
% Knobs (pincat 14x) [NO]
truthpath = 'pincat_full.mat';
dinopath  = 'data_TMI/recon_pincat_full_par6.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par7.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par8.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 1.00;  % Gamma < 1 to decrease contrast
alphaE = 0.80;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 1:128;
cols = 1:128;
outpath = './globalsip_errors_pincat_14x.gif';
%}

% Load data
truthdata = load(truthpath,'Xtrue');
dinodata  = load(dinopath,'Lhat','Shat');
rpcapath  = load(rpcapath,'L0','S0');
ktslrpath = load(ktslrpath,'X0');
Xref = truthdata.Xtrue;
Xdin = reshape(dinodata.Lhat + dinodata.Shat,size(Xref));
Xlps = reshape(rpcapath.L0 + rpcapath.S0,size(Xref));
Xkts = reshape(ktslrpath.X0,size(Xref));
[ny, nx, nt] = size(Xref);

% Check NRMSEs
NRMSEdin = 100 * norm(Xdin(:) - Xref(:)) / norm(Xref(:));
NRMSElps = 100 * norm(Xlps(:) - Xref(:)) / norm(Xref(:));
NRMSEkts = 100 * norm(Xkts(:) - Xref(:)) / norm(Xref(:));

% Format images for display
M      = max(abs(Xref(:)));
Xrefd = min(max(alphaX * abs(Xref) / M,0),1).^gammaX;
Xdin2d = min(max(alphaX * abs(Xdin) / M,0),1).^gammaX;

% Error maps
Edin = abs(Xref - Xdin);
Elps = abs(Xref - Xlps);
Ekts = abs(Xref - Xkts);

% Format error maps for display
E      = max([Edin(:);
              Elps(:);
              Ekts(:);]) * alphaE;
Edind = E * (min(max(Edin / E,0),1).^gammaE);
Elpsd = E * (min(max(Elps / E,0),1).^gammaE);
Ektsd = E * (min(max(Ekts / E,0),1).^gammaE);
clim   = [0, E];

% Clip data
Xrefd = Xrefd(rows,cols,:);
Xdin2d = Xdin2d(rows,cols,:);
Edind = Edind(rows,cols,:);
Elpsd = Elpsd(rows,cols,:);
Ektsd = Ektsd(rows,cols,:);

% Knobs
%cm  = jet(64);
cm  = parula(64);
%cm = hot(64);
%cm = spectral(64); cm = cm(end/2:end,:);
val  = [255, 255, 255];
gap0 = 1;
gap1 = 5;

% Construct color images
Xrefi = colors2im(Xrefd,gray(256),[0, 1]);
Xdini = colors2im(Xdin2d,gray(256),[0, 1]);
Edini = colors2im(Edind,cm,clim);
Elpsi = colors2im(Elpsd,cm,clim);
Ektsi = colors2im(Ektsd,cm,clim);

% Generate movie
CX = {Xrefi,Xdini};
CE = {Edini,Elpsi,Ektsi};
if ~isempty(strfind(truthpath,'invivo'))
    % Flip invivo images
    CX = flipInvivo(CX);
    CE = flipInvivo(CE);
end
MX = cell2mov(CX,gap0,val);
ME = cell2mov(CE,gap0,val);
M  = cell2mov({MX,ME},gap1,val);

% Play movie
movie.video = M;
movie.Fv    = Fv;
lbl = @(s,n) sprintf('%s (%.1f%%)',s,n);
opts.xlabels  = {'Reference','DINO-KAT', ...
                 lbl('DINO-KAT',NRMSEdin), ...
                 lbl('L+S',NRMSElps), ...
                 lbl('k-t SLR',NRMSEkts)};
%{
opts.xlabels = {'Reference', ...
                'DINO-KAT', ...
                'DINO-KAT error', ...
                'L+S error', ...
                'k-t SLR error'};
%}
opts.ylabels  = {''};
opts.mag      = mag;
opts.fontSize = fontSize;
opts.gap      = gap;
P = PlayMovie(movie,opts);
P.SaveGIF(outpath);
P.Close();

%% [IMAGE] otazo dictionary

% Knobs
dinopath = 'data_TMI/recon_otazo_full_par2.mat'; % DINO-KAT
fontSize  = 12;
sep       = [0, 0, 0];

% Load data
dinodata = load(dinopath);

% Initial dictionary
D0 = reshape(dctmtx(320),8,8,5,320);
D0 = reshape(squeeze(D0(:,:,1,:)),64,320);

% Learned dictioanry
addpath('deps_lassi');
D = reshape(dinodata.Dhat,8,8,5,320);
D = reshape(squeeze(D(:,:,1,:)),64,320);

%{
% Sort atoms
usage = sum((dinodata.Bhat ~= 0),2); % Count
%usage = sum(abs(dinodata.Bhat).^2,2); % Energy
[~, order] = sort(usage,'descend');
D = D(:,order);
%}

% Initial atoms
T0(:,:,1) = tilePatches(D0,[8, 8],[16, 20],1,sep(1),true);
T0(:,:,2) = tilePatches(D0,[8, 8],[16, 20],1,sep(2),true);
T0(:,:,3) = tilePatches(D0,[8, 8],[16, 20],1,sep(3),true);
T0 = im2uint8(T0);

% Atom magnitudes
Tm(:,:,1) = tilePatches(abs(D),[8, 8],[16, 20],1,sep(1),true);
Tm(:,:,2) = tilePatches(abs(D),[8, 8],[16, 20],1,sep(2),true);
Tm(:,:,3) = tilePatches(abs(D),[8, 8],[16, 20],1,sep(3),true);
Tm = im2uint8(Tm);

% Atom real-parts
Tr(:,:,1) = tilePatches(real(D),[8, 8],[16, 20],1,sep(1),true);
Tr(:,:,2) = tilePatches(real(D),[8, 8],[16, 20],1,sep(2),true);
Tr(:,:,3) = tilePatches(real(D),[8, 8],[16, 20],1,sep(3),true);
Tr = im2uint8(Tr);

% Atom imaginary-parts
Ti(:,:,1) = tilePatches(imag(D),[8, 8],[16, 20],1,sep(1),true);
Ti(:,:,2) = tilePatches(imag(D),[8, 8],[16, 20],1,sep(2),true);
Ti(:,:,3) = tilePatches(imag(D),[8, 8],[16, 20],1,sep(3),true);
Ti = im2uint8(Ti);

% Display initial dictionary
cfigure();
imshow(T0);
title('Initial atoms');
SetFigFontSize(fontSize);
%export_fig -pdf -transparent globalsip_atom_init
imwrite(T0,'globalsip_atom_init.png','png');

% Display atom magnitudes
cfigure();
imshow(Tm);
title('Atom magnitudes');
SetFigFontSize(fontSize);
%export_fig -pdf -transparent globalsip_atom_mag
imwrite(Tm,'globalsip_atom_mag.png','png');

% Display atom real-parts
cfigure();
imshow(Tr);
title('Atom real-parts');
SetFigFontSize(fontSize);
%export_fig -pdf -transparent globalsip_atom_real
imwrite(Tr,'globalsip_atom_real.png','png');

% Display atom imaginary-parts
figure;
imshow(Ti);
title('Atom imaginary-parts');
SetFigFontSize(fontSize);
%export_fig -pdf -transparent globalsip_atom_imag
imwrite(Ti,'globalsip_atom_imag.png','png');

%% [IMAGE] recons + error maps (by dataset)

flipFcn    = @(X) flipdim(flipdim(X,1),2);
flipInvivo = @(C) cellfun(flipFcn,C,'UniformOutput',false);

% Knobs (otazo 8x) [YES]
truthpath = 'otazo_full.mat';
dinopath  = 'data_TMI/recon_otazo_full_par2.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par3.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par4.mat';
f1 = 7;
f2 = 13;
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.90;  % Gamma < 1 to decrease contrast
alphaE = 0.40;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
gammaC = 1.75;  % Colormap gamma "correction"
fpos     = [923, 272];
fontSize = 13;
rows = 1:128;
cols = 1:128;
colobarFcn = @() cobar([],4,2,2);
Emanual = 0.12;

%{
% Knobs (invivo 8x) [YES]
truthpath = 'invivo_full.mat';
dinopath  = 'data_TMI/recon_invivo_full_par2.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par3.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par4.mat';
f1 = 13;
f2 = 45;
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.75;  % Gamma < 1 to decrease contrast
alphaE = 0.48;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
gammaC = 1.50;  % Colormap gamma "correction"
fpos = [923, 396];
fontSize = 14;
rows = 31:160;
cols = 1:90;
colobarFcn = @() cobar([],4,2,2);
Emanual = 0.12;
%}

%{
% Knobs (pincat 9x) [YES]
truthpath = 'pincat_full.mat';
dinopath  = 'data_TMI/recon_pincat_full_par2.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par3.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par4.mat';
f1 = 16;
f2 = 25;
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 1.00;  % Gamma < 1 to decrease contrast
alphaE = 0.80;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
gammaC = 1.20;  % Colormap gamma "correction"
fpos     = [923, 272];
fontSize = 13;
rows = 1:128;
cols = 1:128;
colobarFcn = @() cobar([],4,2,2);
Emanual = 0.30;
%}

% Load data
truthdata = load(truthpath,'Xtrue');
dinodata  = load(dinopath,'Lhat','Shat');
rpcapath  = load(rpcapath,'L0','S0');
ktslrpath = load(ktslrpath,'X0');
Xref = truthdata.Xtrue;
Xdin = reshape(dinodata.Lhat + dinodata.Shat,size(Xref));
Xlps = reshape(rpcapath.L0 + rpcapath.S0,size(Xref));
Xkts = reshape(ktslrpath.X0,size(Xref));
[ny, nx, nt] = size(Xref);

% Extract images
Xref1 = abs(Xref(:,:,f1));
Xdin1 = abs(Xdin(:,:,f1));
Xlps1 = abs(Xlps(:,:,f1));
Xkts1 = abs(Xkts(:,:,f1));
Xref2 = abs(Xref(:,:,f2));
Xdin2 = abs(Xdin(:,:,f2));
Xlps2 = abs(Xlps(:,:,f2));
Xkts2 = abs(Xkts(:,:,f2));

% Format images for display
M      = max(abs(Xref(:)));
Xref1d = min(max(alphaX * Xref1 / M,0),1).^gammaX;
Xref2d = min(max(alphaX * Xref2 / M,0),1).^gammaX;
Xdin1d = min(max(alphaX * Xdin1 / M,0),1).^gammaX;
Xdin2d = min(max(alphaX * Xdin2 / M,0),1).^gammaX;

% Error maps
Edin1 = abs(Xref1 - Xdin1);
Elps1 = abs(Xref1 - Xlps1);
Ekts1 = abs(Xref1 - Xkts1);
Edin2 = abs(Xref2 - Xdin2);
Elps2 = abs(Xref2 - Xlps2);
Ekts2 = abs(Xref2 - Xkts2);

% Format error maps for display
E      = max([Edin1(:); Edin2(:); ...
              Elps1(:); Elps2(:); ...
              Ekts1(:); Ekts2(:)]) * alphaE;
if exist('Emanual','var') && ~isempty(Emanual)
    E = Emanual;
end
Edin1d = E * (min(max(Edin1 / E,0),1).^gammaE);
Elps1d = E * (min(max(Elps1 / E,0),1).^gammaE);
Ekts1d = E * (min(max(Ekts1 / E,0),1).^gammaE);
Edin2d = E * (min(max(Edin2 / E,0),1).^gammaE);
Elps2d = E * (min(max(Elps2 / E,0),1).^gammaE);
Ekts2d = E * (min(max(Ekts2 / E,0),1).^gammaE);
clim = [0, E];

%--------------------------------------------------------------------------
% Display reconstructions
%--------------------------------------------------------------------------
% Knobs
%cm  = jet(64);
cm  = parula(64);
%cm = hot(64);
%cm = spectral(64); cm = cm(end/2:end,:);
val  = [255, 255, 255];
gap0 = 1;
gap1 = 5;

% Apply gamma "correction" to colormap
n = size(cm,1);
x = linspace(0,1,n);
cm0 = cm;
for i = 1:3
    cm(:,i) = pchip(x,cm0(:,i),x.^gammaC);
end

% Generate images
Xref1i = uint8(repmat(255 * Xref1d,[1 1 3]));
Xref2i = uint8(repmat(255 * Xref2d,[1 1 3]));
Xdin1i = uint8(repmat(255 * Xdin1d,[1 1 3]));
Xdin2i = uint8(repmat(255 * Xdin2d,[1 1 3]));
Edin1i = colors2im(Edin1d,cm,clim);
Edin2i = colors2im(Edin2d,cm,clim);
Elps1i = colors2im(Elps1d,cm,clim);
Elps2i = colors2im(Elps2d,cm,clim);
Ekts1i = colors2im(Ekts1d,cm,clim);
Ekts2i = colors2im(Ekts2d,cm,clim);

% Clip images
Xref1i = Xref1i(rows,cols,:);
Xdin1i = Xdin1i(rows,cols,:);
Edin1i = Edin1i(rows,cols,:);
Elps1i = Elps1i(rows,cols,:);
Ekts1i = Ekts1i(rows,cols,:);
Xref2i = Xref2i(rows,cols,:);
Xdin2i = Xdin2i(rows,cols,:);
Edin2i = Edin2i(rows,cols,:);
Elps2i = Elps2i(rows,cols,:);
Ekts2i = Ekts2i(rows,cols,:);

for i = 1:3
    % Construct images
    xlabels = {'Reference','DINO-KAT','DINO-KAT error','L+S error','k-t SLR error'};
    if i == 1
        CX      = {Xref1i,Xdin1i; Xref2i,Xdin2i}; % Both
        CE      = {Edin1i,Elps1i,Ekts1i; Edin2i,Elps2i,Ekts2i}; % Both
        ylabels = {sprintf('F%d',f1), sprintf('F%d',f2)}; % Both
    elseif i == 2
        CX      = {Xref1i,Xdin1i}; % First
        CE      = {Edin1i,Elps1i,Ekts1i}; % First
        ylabels = {sprintf('F%d',f1)}; % First
    elseif i == 3
        CX      = {Xref2i,Xdin2i}; % Second
        CE      = {Edin2i,Elps2i,Ekts2i}; % Second
        ylabels = {sprintf('F%d',f2)}; % Second
    end
    if ~isempty(strfind(truthpath,'invivo'))
        % Flip invivo images
        CX = flipInvivo(CX);
        CE = flipInvivo(CE);
    end
    MX = cell2mov(CX,gap0,val);
    ME = cell2mov(CE,gap0,val);
    M  = cell2mov({MX,ME},gap1,val);
    
    % Display results
    cfigure(fpos);
    imagesc(M);
    axis equal; axis tight; axis off;
    addXLabels(xlabels,[],0.01,'top' ,'FontSize',fontSize);
    addYLabels(ylabels,[],0.002,'left','FontSize',fontSize);
    caxis(clim);
    colormap(cm);
    colobarFcn();
    SetFigFontSize(fontSize);
    isDataset = @(t,r) strncmpi(t,r,numel(r));
    if i == 1
        if isDataset(truthpath,'otazo')
            export_fig -pdf -transparent globalsip_errors_otazo_8x_F7_F13
        elseif isDataset(truthpath,'invivo')
            export_fig -pdf -transparent globalsip_errors_invivo_8x_F13_F45
        elseif isDataset(truthpath,'pincat')
            export_fig -pdf -transparent globalsip_errors_pincat_9x_F16_F25
        end
    elseif i == 2
        if isDataset(truthpath,'otazo')
            export_fig -pdf -transparent globalsip_errors_otazo_8x_F7
        elseif isDataset(truthpath,'invivo')
            export_fig -pdf -transparent globalsip_errors_invivo_8x_F13
        elseif isDataset(truthpath,'pincat')
            export_fig -pdf -transparent globalsip_errors_pincat_9x_F16
        end
    elseif i == 3
        if isDataset(truthpath,'otazo')
            export_fig -pdf -transparent globalsip_errors_otazo_8x_F13
        elseif isDataset(truthpath,'invivo')
            export_fig -pdf -transparent globalsip_errors_invivo_8x_F45
        elseif isDataset(truthpath,'pincat')
            export_fig -pdf -transparent globalsip_errors_pincat_9x_F25
        end
    end
    close();
end
%--------------------------------------------------------------------------

%% [TABLE] sweep otazo

% Knobs
data(1).inpath = './LASSI_full_Sonly_par1.mat'; % L+S visual
%data(1).inpath = './LASSI_full_Sonly_par2.mat'; % L+S mmse
%data(2).inpath = './RPCA_full_par1.mat'; % L+S visual
data(2).inpath = './RPCA_full_par1b.mat'; % L+S mmse
data(3).inpath = './KTSLR_full_par2.mat'; % k-t SLR mmse
data(1).label = 'DINO-KAT';
data(2).label = 'L+S';
data(3).label = 'k-t SLR';
outpath = 'globalsip_sweep_otazo.pdf';
sortIdx = 4;
trials  = [3];

% Load data
nData = numel(data);
cm    = linspecer(nData);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results (otazo)
    data(i).ps    = datai.vars.ps; %#ok
    data(i).nrmse = nanmean(datai.stats.nrmse(:,trials,end),2); %#ok
end

% Compute effective accelerations
accel = round(1 ./ data(1).ps);
astr  = arrayfun(@(a)sprintf('%dx',a),accel,'UniformOutput',false);

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
labelFcn     = @(d) sprintf('%dx',d);
labels.strs  = [{'Acceleration'}, arrayfun(labelFcn,accel,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    rows(i).name = sprintf('NRMSE (%s) %%',data(i).label);
    rows(i).data = 100 * data(i).nrmse;
end
rows(nData).sep = true;
rows(nData + 1).name = 'Improvement over L+S (dB)';
rows(nData + 1).data = 20 * log10(data(2).nrmse ./ data(1).nrmse);
rows(nData + 2).name = 'Improvement over k-t SLR (dB)';
rows(nData + 2).data = 20 * log10(data(3).nrmse ./ data(1).nrmse);
[rows.fmt]    = deal('%.1f');
[rows.nAlign] = deal('center');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

% Sort data, if requested
if exist('sortIdx','var') && ~isempty(sortIdx)
    [~, idx] = sort(cellfun(@(n)n(sortIdx),{data.nrmse}),'descend');
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([282, 251]);

% Knobs
mk = 'oxd*';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = plot(accel,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Acceleration');
ylabel('NRMSE');
title('Cardiac perfusion');
legend(phndl,data(idx).label,'Location','NW');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
SetFigFontSize(fontSize);
export_fig('-pdf','-transparent',outpath);
%--------------------------------------------------------------------------

%% [TABLE] sweep invivo

% Knobs
%data(1).inpath = './invivo_LASSI_full_Sonly_par1.mat'; % r = 5, S = ktSLR
data(1).inpath = './invivo_LASSI_full_Sonly_par1b.mat'; % r = 1, S = ktSLR
data(2).inpath = './invivo_RPCA_full_par1b.mat';
data(3).inpath = './invivo_KTSLR_full_par1.mat';
data(1).label = 'DINO-KAT';
data(2).label = 'L+S';
data(3).label = 'k-t SLR';
outpath = 'globalsip_sweep_invivo.pdf';
sortIdx = 3;
trials  = [5];

% Load data
nData = numel(data);
cm    = linspecer(nData);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results (otazo)
    data(i).nLines = datai.vars.nLines; %#ok
    data(i).nrmse  = nanmean(datai.stats.nrmse(:,trials,end),2); %#ok
end

% Compute effective accelerations
%{
nCols = numel(data(1).nLines);
accel = zeros(1,nCols);
for i = 1:nCols
    Y = generateInvivoData(data(1).nLines(i),inf,trials(1),'invivo_full.mat');
    accel(i) = round(numel(Y) / nnz(Y));
end
accel %#ok
%}
accel = [23, 12, 8, 6, 5, 4];

% Flip the order of the data
accel = fliplr(accel);
for i = 1:nData
    data(i).nLines = fliplr(data(i).nLines(:)'); %#ok
    data(i).nrmse  = fliplr(data(i).nrmse(:)'); %#ok
end
astr = arrayfun(@(a)sprintf('%dx',a),accel,'UniformOutput',false);

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
%labelFcn    = @(d) sprintf('%d',d);
%labels.strs = [{'# Lines'}, arrayfun(labelFcn,data(1).nLines,'UniformOutput',false)];
labelFcn     = @(d) sprintf('%dx',d);
labels.strs  = [{'Acceleration'}, arrayfun(labelFcn,accel,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    rows(i).name = sprintf('NRMSE (%s) %%',data(i).label);
    rows(i).data = 100 * data(i).nrmse;
end
rows(nData).sep = true;
rows(nData + 1).name = 'Improvement over L+S (dB)';
rows(nData + 1).data = 20 * log10(data(2).nrmse ./ data(1).nrmse);
rows(nData + 2).name = 'Improvement over k-t SLR (dB)';
rows(nData + 2).data = 20 * log10(data(3).nrmse ./ data(1).nrmse);
[rows.fmt]    = deal('%.1f');
[rows.nAlign] = deal('center');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

% Sort data, if requested
if exist('sortIdx','var') && ~isempty(sortIdx)
    [~, idx] = sort(cellfun(@(n)n(sortIdx),{data.nrmse}),'descend');
end

%--------------------------------------------------------------------------
% Plot result
%--------------------------------------------------------------------------
cfigure([282, 251]);

% Knobs
mk = 'oxd*';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = semilogx(accel,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Acceleration');
ylabel('NRMSE');
title('Myocardial perfusion');
legend(phndl,data(idx).label,'Location','NW');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
SetFigFontSize(fontSize);
export_fig('-pdf','-transparent',outpath);
%--------------------------------------------------------------------------

%% [TABLE] sweep pincat

% Knobs
%data(1).inpath = './pincat_LASSI_full_Sonly_par1n.mat'; % r = 5, S = ktSLR
data(1).inpath = './pincat_LASSI_full_Sonly_par1nb.mat'; % r = 1, S = ktSLR
data(2).inpath = './pincat_RPCA_full_par1nb.mat';
data(3).inpath = './pincat_KTSLR_full_par1nb.mat';
data(1).label = 'DINO-KAT';
data(2).label = 'L+S';
data(3).label = 'k-t SLR';
outpath = 'globalsip_sweep_pincatn.pdf';
sortIdx = 3;
trials  = [4];

% Load data
nData = numel(data);
cm    = linspecer(nData);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results (otazo)
    data(i).nLines = datai.vars.nLines; %#ok
    data(i).nrmse  = nanmean(datai.stats.nrmse(:,trials,end),2); %#ok
end

% Compute effective accelerations
%{
nCols = numel(data(1).nLines);
accel = zeros(1,nCols);
for i = 1:nCols
    Y = generateInvivoData(data(1).nLines(i),45,trials(1),'pincat_full.mat');
    accel(i) = round(numel(Y) / nnz(Y));
end
accel %#ok
%}
accel = [27, 14, 9, 7, 6, 5];

% Flip the order of the data
accel = fliplr(accel);
for i = 1:nData
    data(i).nLines = fliplr(data(i).nLines(:)'); %#ok
    data(i).nrmse  = fliplr(data(i).nrmse(:)'); %#ok
end
astr = arrayfun(@(a)sprintf('%dx',a),accel,'UniformOutput',false);

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
%labelFcn    = @(d) sprintf('%d',d);
%labels.strs = [{'# Lines'}, arrayfun(labelFcn,data(1).nLines,'UniformOutput',false)];
labelFcn     = @(d) sprintf('%dx',d);
labels.strs  = [{'Acceleration'}, arrayfun(labelFcn,accel,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    rows(i).name = sprintf('NRMSE (%s) %%',data(i).label);
    rows(i).data = 100 * data(i).nrmse;
end
rows(nData).sep = true;
rows(nData + 1).name = 'Improvement over L+S (dB)';
rows(nData + 1).data = 20 * log10(data(2).nrmse ./ data(1).nrmse);
rows(nData + 2).name = 'Improvement over k-t SLR (dB)';
rows(nData + 2).data = 20 * log10(data(3).nrmse ./ data(1).nrmse);
[rows.fmt]    = deal('%.1f');
[rows.nAlign] = deal('center');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

% Sort data, if requested
if exist('sortIdx','var') && ~isempty(sortIdx)
    [~, idx] = sort(cellfun(@(n)n(sortIdx),{data.nrmse}),'descend');
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([282, 251]);

% Knobs
mk = 'oxd*';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = semilogx(accel,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Acceleration');
ylabel('NRMSE');
title('PINCAT phantom');
legend(phndl,data(idx).label,'Location','NW');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
SetFigFontSize(fontSize);
export_fig('-pdf','-transparent',outpath);
%--------------------------------------------------------------------------
