%% otazo data

% Knobs
ps = 1 / 8; % 1 ./ [4, 8, 12, 16, 20, 24]
SNR = inf;
seed = 1;
inpath = 'data/otazo_full.mat';

% Generete undersampled data
rng(seed);
[Y, Xtrue, Xfft, b1] = generateCardiacPerfData2(ps,SNR,inpath);
M = any(Y,4);

% Baseline NRMSE
NRMSE = norm(Xfft(:) - Xtrue(:)) / norm(Xtrue(:)) %#ok

%{
% Display baseline
PlayMovie(cat(2,Xtrue,Xfft));
%}

% Compare operators
frames = 1:5;
nf = numel(frames);
[ny, nx, nt] = size(Xtrue);
A  = Emat_xyt(M,b1,[ny, nx, nt]);
Af = Emat_xyt(M(:,:,frames),b1,[ny, nx, nf]);
y  = A * Xtrue;
yf = Af * Xtrue(:,:,frames);
x  = reshape(A' * y,[ny, nx, nt]);
xf = reshape(Af' * yf,[ny, nx, nf]);
erry = norm(yf(:) - vec(y(:,:,frames,:))) %#ok
errx = norm(vec(x(:,:,frames)) - xf(:)) %#ok

%% invivo data

% Knobs
nLines = 15; % [5, 10, 15, 20, 25, 30];
SNR = inf;
seed = 1;
inpath = 'data/invivo_full.mat';

% Generete undersampled data
rng(seed);
[Y, Xtrue, Xfft] = generateInvivoData2(nLines,SNR,inpath);
M = (abs(Y) ~= 0);

% Baseline NRMSE
NRMSE = norm(Xfft(:) - Xtrue(:)) / norm(Xtrue(:)) %#ok

%{
% Display baseline
PlayMovie(cat(2,Xtrue,Xfft));
%}

% Compare operators
frames = 1:5;
nf = numel(frames);
[ny, nx, nt] = size(Xtrue);
A  = Afft(M,[ny, nx, nt]);
Af = Afft(M(:,:,frames),[ny, nx, nf]);
y  = A * Xtrue;
yf = Af * Xtrue(:,:,frames);
x  = reshape(A' * y,[ny, nx, nt]);
xf = reshape(Af' * yf,[ny, nx, nf]);
erry = norm(yf(:) - vec(y(:,:,frames,:))) %#ok
errx = norm(vec(x(:,:,frames)) - xf(:)) %#ok

%% PINCAT data

% Knobs
nLines = 18; % [5, 10, 15, 20, 25, 30];
SNR = 45;
seed = 1;
inpath = 'data/pincat_full.mat';

% Generete undersampled data
rng(seed);
[Y, Xtrue, Xfft] = generateInvivoData2(nLines,SNR,inpath);
M = (abs(Y) ~= 0);

% Baseline NRMSE
NRMSE = norm(Xfft(:) - Xtrue(:)) / norm(Xtrue(:)) %#ok

%{
% Display baseline
PlayMovie(cat(2,Xtrue,Xfft));
%}

% Compare operators
frames = 1:5;
nf = numel(frames);
[ny, nx, nt] = size(Xtrue);
A  = Afft(M,[ny, nx, nt]);
Af = Afft(M(:,:,frames),[ny, nx, nf]);
y  = A * Xtrue;
yf = Af * Xtrue(:,:,frames);
x  = reshape(A' * y,[ny, nx, nt]);
xf = reshape(Af' * yf,[ny, nx, nf]);
erry = norm(yf(:) - vec(y(:,:,frames,:))) %#ok
errx = norm(vec(x(:,:,frames)) - xf(:)) %#ok

%% [PLOT] *_onlineDls_par*

% Knobs
%inpath = 'invivo_onlineDls_par1.mat';
%inpath = 'invivo_onlineDls_par1b.mat';
%inpath = 'invivo_onlineDls_par1c.mat';
%inpath = 'invivo_onlineDls_par1d.mat';
%inpath = 'invivo_onlineDls_par2.mat';
%inpath = 'invivo_onlineDls_par4.mat';
%inpath = 'pincat_onlineDls_par1.mat';
%inpath = 'otazo_onlineDls_par1.mat';
%inpath = 'otazo_onlineDls_par2.mat';
%inpath = 'otazo_onlineDls_par4.mat';
inpath = 'otazo_onlineDls_par5.mat';
SAVE_FIGURE = true;

% Load data (D, Xhat, err, stats, vars)
load(inpath);
[~, name, ~] = fileparts(inpath);
[vars, opts] = eval(sprintf('%s();',name));
IS_OTAZO = isfield(vars,'ps');

% Optimal parameters
NRMSE = nanmean(err.NRMSE,2); % average over seed
[mnrmse, idx] = min(NRMSE(:));
[ip, is, il, im, ir, ig] = ind2sub(size(NRMSE),idx);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([1190, 343]);
ax = [];

% ps/nLInes
subplot(2,5,1);
if IS_OTAZO
    % ps
    if (numel(vars.ps) > 1), ax(end + 1) = gca(); end
    plot(vars.ps,squeeze(NRMSE(:,is,il,im,ir,ig)),'b-o');
    xlabel('ps');
    ylabel('NRMSE');
    title(sprintf('ps = %.3f',vars.ps(ip)));
    axis tight; padAxis();
else
    % nLines
    if (numel(vars.nLines) > 1), ax(end + 1) = gca(); end
    plot(vars.nLines,squeeze(NRMSE(:,is,il,im,ir,ig)),'b-o');
    xlabel('nLines');
    ylabel('NRMSE');
    title(sprintf('nLines = %.0f',vars.nLines(ip)));
    axis tight; padAxis();
end

% lambda
subplot(2,5,2);
if (numel(vars.lambda) > 1), ax(end + 1) = gca(); end
semilogx(vars.lambda,squeeze(NRMSE(ip,is,:,im,ir,ig)),'b-o');
xlabel('lambda');
ylabel('NRMSE');
title(sprintf('lambda = %.3f',vars.lambda(il)));
axis tight; padAxis();

% mu
subplot(2,5,3);
if (numel(vars.mu) > 1), ax(end + 1) = gca(); end
semilogx(vars.mu,squeeze(NRMSE(ip,is,il,:,ir,ig)),'b-o');
xlabel('mu');
ylabel('NRMSE');
title(sprintf('mu = %.3f',vars.mu(im)));
axis tight; padAxis();

% dr
subplot(2,5,4);
if (numel(vars.dr) > 1), ax(end + 1) = gca(); end
plot(vars.dr,squeeze(NRMSE(ip,is,il,im,:,ig)),'b-o');
xlabel('dr');
ylabel('NRMSE');
title(sprintf('dr = %d',vars.dr(ir)));
axis tight; padAxis();

% gamma
subplot(2,5,5);
%if (numel(vars.gamma) > 1), ax(end + 1) = gca(); end
plot(vars.gamma,squeeze(NRMSE(ip,is,il,im,ir,:)),'b-o'); hold on;
plot(vars.gamma2,squeeze(NRMSE(ip,is,il,im,ir,:)),'r-x'); hold on;
xlabel('gamma, gamma2');
ylabel('NRMSE');
title(sprintf('gamma = %.3f, gamma2 = %.3f',vars.gamma(ig),vars.gamma2(ig)));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,5,6);
plotOnlineDlsMetric2(squeeze(stats.cost(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,vars.nReps);
xlabel('time + iteration');
ylabel('cost');
title(sprintf('cost, (T, dt, np) = (%d, %d, %d)',vars.T,vars.dt,vars.np));
axis tight; padAxis();

% Optimal NRMSE trajectory
subplot(2,5,7);
plotOnlineDlsMetric2(squeeze(stats.nrmse(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,vars.nReps);
xlabel('time + iteration');
ylabel('NRMSE');
title(sprintf('NRMSE = %.5f',mnrmse));
axis tight; padAxis();

% Optimal sparsity trajectory
subplot(2,5,8);
plotOnlineDlsMetric2(squeeze(stats.sparsity(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,vars.nReps);
xlabel('time + iteration');
ylabel('sparsity (%)');
title('sparsity');
axis tight; padAxis();

% Optimal deltaX trajectory
subplot(2,5,9);
plotOnlineDlsMetric2(squeeze(stats.deltaX(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,vars.nReps);
xlabel('time + iteration');
ylabel('deltaX');
title(sprintf('deltaX, pgap = [%d, %d, %d]',opts.pgap));
set(gca,'YScale','log');
axis tight; padAxis();

% Optimal deltaD trajectory
subplot(2,5,10);
plotOnlineDlsMetric2(squeeze(stats.deltaD(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,vars.nReps);
xlabel('time + iteration');
ylabel('deltaD');
title(sprintf('deltaD, pdim = [%d, %d, %d]',opts.pdim));
set(gca,'YScale','log');
axis tight; padAxis();

% Match axes
if ~isempty(ax), matchAxes(ax,[],'y'); end

% Save figure, if requested
if SAVE_FIGURE
    [~, name, ~] = fileparts(inpath);
    export_fig('-pdf','-transparent',name);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [TABLE] *_onlineDls_par*

%{
% Otazo #1 (baseline init)
inpaths = {'otazo_onlineDls_par3.mat' , 'Online [nReps = 1]';
           'otazo_onlineDls_par3a.mat', 'Online [nReps = 2]';
           'otazo_onlineDls_par3b.mat', 'Online [DCT]'};
%}
%{
% Otazo #2 (baseline init)
inpaths = {'otazo_onlineDls_par3e.mat', 'Online [nReps = 1]';
           'otazo_onlineDls_par3f.mat', 'Online [nReps = 2]';
           'otazo_onlineDls_par3g.mat', 'Online [DCT]'};
%}
%{
% Otazo #3 (baseline init)
inpaths = {'otazo_onlineDls_par6.mat' , 'Online [nReps = 1]';
           'otazo_onlineDls_par6a.mat', 'Online [nReps = 2]';
           'otazo_onlineDls_par6b.mat', 'Online [DCT]'};
%}
%{
% Invivo (baseline init)
inpaths = {'invivo_onlineDls_par3.mat', 'Online [nReps = 3]'};
%}
%{
% Invivo (L+S init)
inpaths = {'invivo_onlineDls_par5.mat' , 'Online [nReps = 1]';
           'invivo_onlineDls_par5b.mat', 'Online [DCT]'};
%}
%{
% PINCAT (baseline init)
inpaths = {'pincat_onlineDls_par3.mat' , 'Online [nReps = 3]';
           'pincat_onlineDls_par3b.mat', 'Online [nReps = 1]';
           'pincat_onlineDls_par3c.mat', 'Online [DCT]'};
%}

% Load data
nData = size(inpaths,1);
data = struct();
for k = 1:nData
    datak = load(inpaths{k,1});
    nrmse = squeeze(datak.stats.nrmse); % 6 x ni
    data(k).label = inpaths{k,2};
    data(k).nrmse = 100 * nrmse(:,end)';
end
IS_OTAZO = isfield(datak.vars,'ps');
if IS_OTAZO
    ps = datak.vars.ps;
else
    nLines = datak.vars.nLines;
end

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
if IS_OTAZO
    labelFcn = @(ps) sprintf('1/%d',round(1 / ps));
    labels.strs = ['ps', arrayfun(labelFcn,ps,'UniformOutput',false)];
else
    labelFcn = @(nLines) sprintf('%d',nLines);
    labels.strs = ['nLines', arrayfun(labelFcn,nLines,'UniformOutput',false)];
end

% Row data
rows = struct();
for k = 1:nData
    rows(k).name = data(k).label;
    rows(k).data = data(k).nrmse;
end
[rows.fmt]    = deal('%.2f%%');
[rows.nAlign] = deal('right');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

%% 0/5 [PLOT, IMAGE, VIDEOS] recon
%  Otazo: accel = [4x ,  8x, 12x, 16x, 20x, 24x]
% Invivo: accel = [23x, 12x,  8x,  6x,  5x,  4x]
% PINCAT: accel = [27x, 14x,  9x,  7x,  6x,  5x]

% Universal constants
FRAME_RATE = 12;
FONT_SIZE = 14;
CMX = gray(256);
CME = parula(64);
GAP = 1;
VAL = 255;


% Otazo
BASE = '/Users/Brian/Desktop/Online DINO-KAT/ISMRM';
inpaths = {'batch_otazo_rpca'     , 'L+S';
           'otazo_onlineDls_par3' , 'Online DINO-KAT';      % 1 pass
%          'otazo_onlineDls_par3a', 'Online DINO-KAT';      % 2 passes
%          'batch_otazo_fft'      , 'Baseline';
           'otazo_onlineDls_par3b', 'Online (DCT)'};
%IDX = 5; % ps = 1 ./ [4, 8, 12, 16, 20, 24]
IDX = 3; % ps = 1 ./ [4, 8, 12, 16, 20, 24]
GAMMAX = 0.9;
GAMMAE = 1.5;
ALPHAE = 0.9;
%FRAMES = [3, 31];
FRAMES = [4, 11];


%{
% PINCAT
BASE = '/Users/Brian/Desktop/Online DINO-KAT/ISMRM';
inpaths = {'batch_pincat_rpca'     , 'L+S';
           'pincat_onlineDls_par3b', 'Online DINO-KAT';     % 1 pass
%          'pincat_onlineDls_par3' , 'Online DINO-KAT';     % 3 passes
%          'batch_pincat_fft'      , 'Baseline';
           'pincat_onlineDls_par3c', 'Online (DCT)'};
%IDX = 5; % Invivo: nLines   = [5, 10, 15, 20, 25, 30]
IDX = 5; % Invivo: nLines   = [5, 10, 15, 20, 25, 30]
GAMMAX = 1.0;
GAMMAE = 1.3;
ALPHAE = 0.77;
%FRAMES = [5, 15];
FRAMES = [3, 29];
%}

% Load ground truth
vars = eval(sprintf('%s();',inpaths{end,1}));
datat = load(vars.inpath);
Xtrue = datat.Xtrue;
trueMax = max(abs(Xtrue(:)));

% Load recons
nData = size(inpaths,1);
Xhat = cell(1,nData);
for k = 1:nData
    datak = load(sprintf('%s/%s/data%d.mat',BASE,inpaths{k,1},IDX));
    Xhat{k} = datak.Xhat;
end

% All data
X = [{Xtrue}, Xhat{:}];
xlabels = [{'Reference'}, inpaths{:,2}];

% Compute per-frame errors
clear err;
for k = 1:nData
    err(k) = computeErrorMetrics(Xhat{k},Xtrue); %#ok
end

%% 1/5 Recon movie

% Format data
CMXg = gammaCorrectColormap(CMX,GAMMAX);
VX = abs(cell2mat(X));
VX = min(VX / trueMax,1);
VX = colors2im(VX,CMXg,[0, 1]);

% Play movie
movie = struct();
movie.video = VX;
movie.Fv = FRAME_RATE;
opts = struct();
opts.xlabels = xlabels;
opts.fontSize = FONT_SIZE;
PlayMovie(movie,opts);

%{
% Save movie
./otazo_recons.gif
./pincat_recons.gif
%}

%% 2/5 Recon + error movie

% Compute error videos
E = cell(1,nData + 1);
E{1} = zeros(size(Xtrue));
for k = 1:nData
    E{k + 1} = abs(Xhat{k} - Xtrue);
end

% Format X videos
CMXg = gammaCorrectColormap(CMX,GAMMAX);
VX = abs(cell2mat(X));
VX = min(VX / trueMax,1);
VX = colors2im(VX,CMXg,[0, max(VX(:))]);

% Format E videos
CMEg = gammaCorrectColormap(CME,GAMMAE);
VE = cell2mat(E);
VE = colors2im(VE,CMEg,[0, ALPHAE * max(VE(:))]);

% Play movie
movie = struct();
movie.video = cat(1,VX,VE);
movie.Fv = FRAME_RATE;
opts = struct();
opts.xlabels  = xlabels;
opts.ylabels  = {'Errors','Recons'};
opts.fontSize = FONT_SIZE;
PlayMovie(movie,opts);

%{
% Save movie
./otazo_errors.gif
./pincat_errors.gif
%}

%% 3/5 Single-frame recons

% Find best frames
lpsIdx = 1; % L+S method
onlineIdx = 2; % best online method
[~, idx] = sort(err(lpsIdx).pNRMSE - err(onlineIdx).pNRMSE,'descend');
best_frames = idx(1:10) %#ok

% Generate frames
nFrames = numel(FRAMES);
IX = cell(nFrames,nData + 1);
formatFrame = @(X,f) min(abs(X(:,:,f)) / trueMax,1);
for j = 1:nFrames
    IX(j,:) = cellfun(@(Xk) formatFrame(Xk,FRAMES(j)),X,'UniformOutput',false);
end

% Format frames
CMXg = gammaCorrectColormap(CMX,GAMMAX);
IX = cellfun(@(X) colors2im(X,CMXg,[0, 1]),IX,'UniformOutput',false);
IX = cell2mov(IX,GAP,VAL);

% Plot frames
imshow(IX,[0, 1]);
ylabels = arrayfun(@(f) sprintf('Frame %d',f),FRAMES,'UniformOutput',false);
addXLabels(xlabels,[],0.015,'top');
addYLabels(ylabels,[],0.005,'left');
SetFigFontSize(FONT_SIZE);

%{
% Save figure
export_fig('-pdf','-transparent','otazo_recons');
export_fig('-pdf','-transparent','pincat_recons');
%}

%% 4/5 Single-frame recons + errors

% Find best frames
lpsIdx = 1; % L+S method
onlineIdx = 2; % best online method
[~, idx] = sort(err(lpsIdx).pNRMSE - err(onlineIdx).pNRMSE,'descend');
best_frames = idx(1:10) %#ok

% Generate frames
nFrames = numel(FRAMES);
IX = cell(nFrames,nData + 1);
formatFrame = @(X,f) min(abs(X(:,:,f)) / trueMax,1);
for j = 1:nFrames
    IX(j,:) = cellfun(@(Xk) formatFrame(Xk,FRAMES(j)),X,'UniformOutput',false);
end

% Generate errors
IE = cell(nFrames,nData + 1);
for j = 1:nFrames
    f = FRAMES(j);
    Xtruef = Xtrue(:,:,f);
    IE{j,1} = zeros(size(Xtruef));
    for k = 1:nData
        IE{j,k + 1} = abs(Xhat{k}(:,:,f) - Xtruef);
    end
end

% Extract relevant data
IX = IX(:,[1, 3]); % Truth and best online
IE = IE(:,2:end);
XLABELS = [xlabels([1, 3]), xlabels(2:end)];
YLABELS = arrayfun(@(f) sprintf('Frame %d',f),FRAMES,'UniformOutput',false);

% Format frames
CMXg = gammaCorrectColormap(CMX,GAMMAX);
IX = cellfun(@(X) colors2im(X,CMXg,[0, 1]),IX,'UniformOutput',false);
IX = cell2mov(IX,GAP,VAL);

% Format errors
CMEg = gammaCorrectColormap(CME,GAMMAE);
Elim = [0, ALPHAE * max(vec(cellfun(@(X) max(X(:)),IE)))];
IE = cellfun(@(E) colors2im(E,CMEg,Elim),IE,'UniformOutput',false);
IE = cell2mov(IE,GAP,VAL);

% Plot frames + errors
imshow(cell2mov({IX, IE},GAP,VAL),[0, 1]);
addXLabels(XLABELS,[],0.015,'top');
addYLabels(YLABELS,[],0.005,'left');

% Add colorbar #1
ax1 = gca();
ax2 = axes('Position',get(ax1,'Position'),'Visible','off');
set(colorbar('Peer',ax1),'Visible','off');
colorbar('Peer',ax2);
colormap(ax2,CMEg);
caxis(ax2,Elim);

%{
% Add colorbar #2
c = colorbar();
colormap(CMEg);
yTicks = get(c,'YTick');
yTicks2 = (Elim(2) / size(CMEg,1)) * yTicks;
yTickLabels = arrayfun(@(t) sprintf('%.2f',t),yTicks2,'UniformOutput',false);
set(c,'YTickLabel',yTickLabels);
%}

% Set font size
SetFigFontSize(FONT_SIZE);

%{
% Save figure
export_fig('-pdf','-transparent','otazo_errors');
export_fig('-pdf','-transparent','pincat_errors');
%}

%% 5/5 Per-frame errors

% Setup figure
CM = linspecer(nData);
cfigure([468, 273]);
nt = size(Xtrue,3);
st = {'-','.-','*-'};
ms = [6, 12, 6];
lw = [1, 1, 1];

% Sort frames
[~, idx] = sort([err.NRMSE],'descend');

% Plot pNRMSEs
phndl = zeros(1,nData);
for k = 1:nData
    k1 = idx(k);
    k2 = nData + 1 - k;
    phndl(k) = plot(1:nt,100 * err(k1).pNRMSE,st{k2},'Color',CM(k2,:),'Linewidth',lw(k2),'Markersize',ms(k2));
    hold on;
end
xlabel('Frame');
ylabel('NRMSE (%)');
title('Per-frame NRMSEs');
legend(phndl,inpaths{idx,2},'Location','NE');
axis tight; padAxis();

% Adjust x-ticks
xTicks = get(gca,'XTick'); xTicks(1) = 1;
set(gca,'XTick',xTicks);

% Set font size
SetFigFontSize(FONT_SIZE);

%{
% Save figure
export_fig('-pdf','-transparent','otazo_pnrmse');
export_fig('-pdf','-transparent','pincat_pnrmse');
%}
