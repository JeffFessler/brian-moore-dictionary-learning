%% [RAW DATA] videos

BASE = '/Users/Brian/Archive/MATLAB/video_data/BM3D/';

% In order of priority
%path = [BASE, 'gstennis.mat']; % 150
%path = [BASE, 'coastguard.mat']; % 300
%path = [BASE, 'coastguard2.mat']; % 300
%path = [BASE, 'gbus.mat']; % 150
%path = [BASE, 'gflower.mat']; % 150
%path = [BASE, 'gforeman.mat'];
%path = [BASE, 'gsalesman.mat'];
%path = [BASE, 'gbicycle.mat'];
%path = [BASE, 'gmissa.mat'];

idim = [nan, nan, 150]; % first 
T = 5;
dt = 1;

X = loadVideo(path,idim,T,dt);
size(X)

PlayMovie(X);

%% [TABLE,PLOT] interp_par*

%{
% gstennis
inpaths = {
    'interp_par1a.mat', 'cubic (2D)';
    'interp_par1b.mat', 'cubic (3D)';
};
name = 'gstennis_interp.pdf';
%}

%{
% coastguard2
inpaths = {
    'interp_par3a.mat', 'cubic (2D)';
    'interp_par3b.mat', 'cubic (3D)';
};
name = 'coastguard2_interp.pdf';
%}

%{
% gbus
inpaths = {
    'interp_par5a.mat', 'cubic (2D)';
    'interp_par5b.mat', 'cubic (3D)';
};
name = 'gbus_interp.pdf';
%}

%{
% gflower
inpaths = {
    'interp_par7a.mat', 'cubic (2D)';
    'interp_par7b.mat', 'cubic (3D)';
};
name = 'gflower_interp.pdf';
%}

pSortIdx = 3;
metric = 'pSNR'; % NRMSE, pSNR
SAVE_FIGURE = true;

% Load data
data = struct();
for i = 1:size(inpaths,1)
    datai = load(inpaths{i,1},'vars','err');
    data(i).p = squeeze(datai.vars.p);
    data(i).(metric) = squeeze(datai.err.(metric));
    data(i).label = inpaths{i,2};
end
nData = numel(data);

% Sort by metric
val = @(data) data.(metric)(pSortIdx);
[~, idx] = sort(arrayfun(val,data),'descend');

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
labelFcn = @(p) sprintf('%d%%',100 * p);
labels.strs  = [metric, arrayfun(labelFcn,data(1).p,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    ii = idx(i);
    rows(i).name = data(ii).label;
    rows(i).data = data(ii).(metric);
end
[rows.fmt]    = deal('%.2f');
[rows.nAlign] = deal('right');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([281, 231]);
cm = linspecer(nData);

% Plot metric
phndl = zeros(1,nData);
pstr = cell(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = plot(data(ii).p,data(ii).(metric),'-o','Color',cm(i,:));
    pstr{i} = data(ii).label;
    hold on;
end
xlabel('p');
ylabel(metric);
legend(phndl,pstr{:},'Location','SW');
axis tight; padAxis();

% Save figure, if requested
if SAVE_FIGURE
    export_fig('-pdf','-transparent',name);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [PLOT] onlineDls_par*

% Knobs (noiseless)
%inpath = 'onlineDls_par1a.mat';
%inpath = 'onlineDls_par1a2.mat';
%inpath = 'onlineDls_par1a3.mat';
%inpath = 'onlineDls_par1b.mat';
%inpath = 'onlineDls_par1b2.mat';
%inpath = 'onlineDls_par1c.mat';
%inpath = 'onlineDls_par3a.mat';
%inpath = 'onlineDls_par3a2.mat';
%inpath = 'onlineDls_par3b.mat';
%inpath = 'onlineDls_par3c.mat';
%inpath = 'onlineDls_par3g.mat';
%inpath = 'onlineDls_par5a.mat';
%inpath = 'onlineDls_par5a2.mat';
%inpath = 'onlineDls_par5b.mat';
%inpath = 'onlineDls_par5g.mat';
%inpath = 'onlineDls_par7a.mat';
%inpath = 'onlineDls_par7a2.mat';
%inpath = 'onlineDls_par7b.mat';
%inpath = 'onlineDls_par7g.mat';
% Knobs (noisy)
%inpath = 'onlineDls_npar3a.mat';
%inpath = 'onlineDls_npar3b.mat';
%inpath = 'onlineDls_npar3c.mat';
%inpath = 'onlineDls_npar3g.mat';
SAVE_FIGURE = true;

% Load data (err, stats, vars)
load(inpath);
[~, name, ~] = fileparts(inpath);
[vars, opts] = eval(sprintf('%s();',name));

% Optimal parameters
NRMSE = nanmean(err.NRMSE,2); % average over seed
[mnrmse, idx] = min(NRMSE(:));
[ip, is, il, im, ir, ig] = ind2sub(size(NRMSE),idx);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([1190, 343]);
ax = [];

% p
subplot(2,5,1);
if (numel(vars.p) > 1), ax(end + 1) = gca(); end
plot(vars.p,squeeze(NRMSE(:,is,il,im,ir,ig)),'b-o');
xlabel('p');
ylabel('NRMSE');
title(sprintf('p = %.3f',vars.p(ip)));
axis tight; padAxis();

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
if (numel(vars.gamma) > 1), ax(end + 1) = gca(); end
plot(vars.gamma,squeeze(NRMSE(ip,is,il,im,ir,:)),'b-o'); hold on;
plot(vars.gamma2,squeeze(NRMSE(ip,is,il,im,ir,:)),'r-x'); hold on;
xlabel('gamma, gamma2');
ylabel('NRMSE');
title(sprintf('gamma = %.3f, gamma2 = %.3f',vars.gamma(ig),vars.gamma2(ig)));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,5,6);
plotOnlineDlsMetric(squeeze(stats.cost(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,'b-');
xlabel('time + iteration');
ylabel('cost');
title(sprintf('cost, (T, dt, np) = (%d, %d, %d)',vars.T,vars.dt,vars.np));
axis tight; padAxis();

% Optimal NRMSE trajectory
subplot(2,5,7);
plotOnlineDlsMetric(squeeze(stats.nrmse(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,'b-');
xlabel('time + iteration');
ylabel('NRMSE');
title(sprintf('NRMSE = %.5f',mnrmse));
axis tight; padAxis();

% Optimal sparsity trajectory
subplot(2,5,8);
plotOnlineDlsMetric(squeeze(stats.sparsity(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,'b-');
xlabel('time + iteration');
ylabel('sparsity (%)');
title(sprintf('sparsity, idim = [%d, %d, %d]',vars.idim));
axis tight; padAxis();

% Optimal deltaX trajectory
subplot(2,5,9);
plotOnlineDlsMetric(squeeze(stats.deltaX(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,'b-');
xlabel('time + iteration');
ylabel('deltaX');
title(sprintf('deltaX, pgap = [%d, %d, %d]',opts.pgap));
set(gca,'YScale','log');
axis tight; padAxis();

% Optimal deltaD trajectory
subplot(2,5,10);
plotOnlineDlsMetric(squeeze(stats.deltaD(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,'b-');
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

%% [PLOT] dls_par*

% Knobs (noiseless)
%inpath = 'dls_par1a.mat';
%inpath = 'dls_par1b.mat';
%inpath = 'dls_par3a.mat';
%inpath = 'dls_par3b.mat';
%inpath = 'dls_par5b.mat';
%inpath = 'dls_par7b.mat';
% Knobs (noiseless)
%inpath = 'dls_npar3a.mat';
SAVE_FIGURE = true;

% Load data (err, stats, vars)
load(inpath);
[~, name, ~] = fileparts(inpath);
[vars, opts] = eval(sprintf('%s();',name));

% Optimal parameters
NRMSE = nanmean(err.NRMSE,2); % average over seed
[mnrmse, idx] = min(NRMSE(:));
[ip, is, il, im, ir] = ind2sub(size(NRMSE),idx);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([1190, 343]);
ax = [];

% p
subplot(2,5,1);
if (numel(vars.p) > 1), ax(end + 1) = gca(); end
plot(vars.p,squeeze(NRMSE(:,is,il,im,ir)),'b-o');
xlabel('p');
ylabel('NRMSE');
title(sprintf('p = %.3f',vars.p(ip)));
axis tight; padAxis();

% lambda
subplot(2,5,2);
if (numel(vars.lambda) > 1), ax(end + 1) = gca(); end
semilogx(vars.lambda,squeeze(NRMSE(ip,is,:,im,ir)),'b-o');
xlabel('lambda');
ylabel('NRMSE');
title(sprintf('lambda = %.3f',vars.lambda(il)));
axis tight; padAxis();

% mu
subplot(2,5,3);
if (numel(vars.mu) > 1), ax(end + 1) = gca(); end
semilogx(vars.mu,squeeze(NRMSE(ip,is,il,:,ir)),'b-o');
xlabel('mu');
ylabel('NRMSE');
title(sprintf('mu = %.3f',vars.mu(im)));
axis tight; padAxis();

% dr
subplot(2,5,4);
if (numel(vars.dr) > 1), ax(end + 1) = gca(); end
plot(vars.dr,squeeze(NRMSE(ip,is,il,im,:)),'b-o');
xlabel('dr');
ylabel('NRMSE');
title(sprintf('dr = %d',vars.dr(ir)));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,5,6);
plot(1:vars.nIters,squeeze(stats.cost(ip,is,il,im,ir,:)),'b-');
xlabel('iteration');
ylabel('cost');
title(sprintf('cost, (np) = (%d)',vars.np));
axis tight; padAxis();

% Optimal NRMSE trajectory
subplot(2,5,7);
plot(1:vars.nIters,squeeze(stats.nrmse(ip,is,il,im,ir,:)),'b-');
xlabel('iteration');
ylabel('NRMSE');
title(sprintf('NRMSE = %.5f',mnrmse));
axis tight; padAxis();

% Optimal sparsity trajectory
subplot(2,5,8);
plot(1:vars.nIters,squeeze(stats.sparsity(ip,is,il,im,ir,:)),'b-');
xlabel('iteration');
ylabel('sparsity (%)');
title(sprintf('sparsity, idim = [%d, %d, %d]',vars.idim));
axis tight; padAxis();

% Optimal deltaX trajectory
subplot(2,5,9);
plot(1:vars.nIters,squeeze(stats.deltaX(ip,is,il,im,ir,:)),'b-');
xlabel('iteration');
ylabel('deltaX');
title(sprintf('deltaX, pgap = [%d, %d, %d]',opts.pgap));
set(gca,'YScale','log');
axis tight; padAxis();

% Optimal deltaD trajectory
subplot(2,5,10);
plot(1:vars.nIters,squeeze(stats.deltaD(ip,is,il,im,ir,:)),'b-');
xlabel('iteration');
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

%% [TABLE,PLOT] p sweeps
% Plots (p, dr) slices of data

%{
% gstennis
ipaths = {'interp_par1a.mat', 'cubic (2D)';
          'interp_par1b.mat', 'cubic (3D)'};
dpaths = {'onlineDls_par2a.mat', 'Online DINO-KAT'; % r = 5
          'onlineDls_par2b.mat', 'Online DINO-KAT [unitary 1x1]';
          'onlineDls_par2c.mat', 'Online DINO-KAT [fixed D]';
          'onlineDls_par2d.mat', 'Online DINO-KAT'; % r = 1
%         'onlineDls_par2e.mat', 'Online DINO-KAT [nIters0 = 5]'; % r = 5
%         'onlineDls_par2f.mat', 'Online DINO-KAT [nIters0 = 5]'; % r = 1
          'dls_par2a.mat', 'Batch DINO-KAT'; % r = 1
          'dls_par2b.mat', 'Batch DINO-KAT'}; % r = 5
name = 'gstennis_psweep.pdf';
%}

%{
% coastguard2
ipaths = {'interp_par3a.mat', 'cubic (2D)';
          'interp_par3b.mat', 'cubic (3D)'};
dpaths = {'onlineDls_par4a.mat', 'Online DINO-KAT'; % r = 5
          'onlineDls_par4b.mat', 'Online DINO-KAT [unitary 1x1]';
          'onlineDls_par4c.mat', 'Online DINO-KAT [fixed D]';
          'onlineDls_par4d.mat', 'Online DINO-KAT'; % r = 1
%         'onlineDls_par4e.mat', 'Online DINO-KAT [nIters0 = 5]'; % r = 5
%         'onlineDls_par4f.mat', 'Online DINO-KAT [nIters0 = 5]'; % r = 1
          'onlineDls_par4g.mat', 'Online DINO-KAT [unitary 2x2]';
          'dls_par4a.mat', 'Batch DINO-KAT'; % r = 1
          'dls_par4b.mat', 'Batch DINO-KAT'}; % r = 5
name = 'coastguard2_psweep';
%}

%{
% coastguard2 (noisy)
ipaths = {'interp_npar3a.mat', 'cubic (2D)';
          'interp_npar3b.mat', 'cubic (3D)'};
dpaths = {'onlineDls_npar4a.mat', 'Online DINO-KAT'; % r = 5
          'onlineDls_npar4b.mat', 'Online DINO-KAT [unitary 1x1]';
          'onlineDls_npar4c.mat', 'Online DINO-KAT [fixed D]';
          'onlineDls_npar4d.mat', 'Online DINO-KAT'; % r = 1
          'onlineDls_npar4g.mat', 'Online DINO-KAT [unitary 2x2]';
          'dls_npar4a.mat', 'Batch DINO-KAT'; % r = 5
          'dls_npar4b.mat', 'Batch DINO-KAT'}; % r = 1
name = 'coastguard2_noisy_psweep';
%}

%{
% gbus
ipaths = {'interp_par5a.mat', 'cubic (2D)';
          'interp_par5b.mat', 'cubic (3D)'};
dpaths = {'onlineDls_par6a.mat', 'Online DINO-KAT'; % r = 5
          'onlineDls_par6b.mat', 'Online DINO-KAT [unitary 1x1]';
          'onlineDls_par6c.mat', 'Online DINO-KAT [fixed D]';
          'onlineDls_par6d.mat', 'Online DINO-KAT'; % r = 1
%         'onlineDls_par6e.mat', 'Online DINO-KAT [nIters0 = 5]'; % r = 5
%         'onlineDls_par6f.mat', 'Online DINO-KAT [nIters0 = 5]'; % r = 1
          'onlineDls_par6g.mat', 'Online DINO-KAT [unitary 2x2]';
          'dls_par6a.mat', 'Batch DINO-KAT'; % r = 1
          'dls_par6b.mat', 'Batch DINO-KAT'}; % r = 5
name = 'gbus_psweep';
%}

%{
% gflower
ipaths = {'interp_par7a.mat', 'cubic (2D)';
          'interp_par7b.mat', 'cubic (3D)'};
dpaths = {'onlineDls_par8a.mat', 'Online DINO-KAT'; % r = 5
          'onlineDls_par8b.mat', 'Online DINO-KAT [unitary 1x1]';
          'onlineDls_par8c.mat', 'Online DINO-KAT [fixed D]';
          'onlineDls_par8d.mat', 'Online DINO-KAT'; % r = 1
%         'onlineDls_par8e.mat', 'Online DINO-KAT [nIters0 = 5]'; % r = 5
%         'onlineDls_par8f.mat', 'Online DINO-KAT [nIters0 = 5]'; % r = 1
          'onlineDls_par8g.mat', 'Online DINO-KAT [unitary 2x2]';
          'dls_par8a.mat', 'Batch DINO-KAT'; % r = 1
          'dls_par8b.mat', 'Batch DINO-KAT'}; % r = 5
name = 'gflower_psweep';
%}

pSortIdx = 3;
metric = 'pSNR'; % NRMSE, pSNR
SAVE_FIGURE = true;

% Load data
data = struct();
for i = 1:size(ipaths,1)
    % Interp data
    datai = load(ipaths{i,1},'vars','err');
    data(i).p = squeeze(datai.vars.p);
    data(i).(metric) = squeeze(datai.err.(metric));
    data(i).label = ipaths{i,2};
end
for j = 1:size(dpaths,1)
    % DINO-KAT data
    dataj = load(dpaths{j,1},'vars','err');
    errData = dataj.err.(metric); % 7D
    pData = squeeze(dataj.vars.p);
    drData = dataj.vars.dr;
    for k = 1:size(errData,5)
        data(end + 1).p = pData; %#ok
        data(end).(metric) = squeeze(errData(:,:,:,:,k,:,:));
        if isnan(drData(k))
            % No dr
            data(end).label = sprintf('%s',dpaths{j,2});
        else
            % Add dr
            data(end).label = sprintf('%s [r = %d]',dpaths{j,2},drData(k));
        end
    end
end
nData = numel(data);

% Sort by metric
val = @(data) data.(metric)(pSortIdx);
[~, idx] = sort(arrayfun(val,data),'descend');

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
labelFcn = @(p) sprintf('%d%%',100 * p);
labels.strs  = [metric, arrayfun(labelFcn,data(1).p,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    rows(i).name = data(i).label;
    rows(i).data = data(i).(metric);
end
[rows.fmt]    = deal('%.2f');
[rows.nAlign] = deal('right');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows,pSortIdx,3,true);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([492, 231]);
cm = linspecer(nData);

% Plot metric
phndl = zeros(1,nData);
pstr = cell(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = plot(data(ii).p,data(ii).(metric),'-o','Color',cm(i,:));
    pstr{i} = data(ii).label;
    hold on;
end
xlabel('p');
ylabel(metric);
legend(phndl,pstr{:},'Location','EastOutside');
axis tight; padAxis();

% Save figure, if requested
if SAVE_FIGURE
    savePDF(name);
end
%--------------------------------------------------------------------------

%% [PLOT] dictionaries

BASE_PATH = '/Users/Brian/Desktop/Online DINO-KAT/TCI-inpaint';

%{
% coastguard2
% p = [0.5, 0.6, 0.7, 0.8, 0.9];
inpath = 'onlineDls_par4a'; % r = 5
outpath = 'tci_coastguard2_dict_learned';
idx = 1;
%}

%{
% coastguard2 noisy
% p = [0.5, 0.6, 0.7, 0.8, 0.9];
inpath = 'onlineDls_npar4a'; % r = 5
outpath = 'tci_coastguard2_noisy_dict_learned';
idx = 1;
%}

%{
% gbus
% p = [0.5, 0.6, 0.7, 0.8, 0.9];
inpath = 'onlineDls_par6a'; % r = 5
outpath = 'tci_gbus_dict_learned';
idx = 1;
%}

%{
% gflower
% p = [0.5, 0.6, 0.7, 0.8, 0.9];
inpath = 'onlineDls_par8a'; % r = 5
outpath = 'tci_gflower_dict_learned';
idx = 1;
%}

gap = 1; val = 1;
col = 4;

% Load learned dictionary
data = load(sprintf('%s/%s/data%d.mat',BASE_PATH,inpath,idx));
Dhat = data.D;

% Construct initial dictionary
D0 = dctmtx(320)';

% Display dictionaries
D0 = reshape(D0,[8, 8, 5, 320]);
Dhat = reshape(Dhat,[8, 8, 5, 320]);
C0a = cell(16,20);
C0b = cell(16,20);
Ca = cell(16,20);
Cb = cell(16,20);
for k = 1:320
    %D0k = reshape(D0(:,:,:,k),[64, 5]);
    %[U0k, S0k, V0k] = svd(D0k);
    %C0a{k} = formatDict(reshape(U0k(:,1),[8, 8]));
    C0a{k} = formatDict(D0(:,:,1,k));
    
    C0b{k} = formatDict(squeeze(D0(:,col,:,k)));
    
    %Dk = reshape(Dhat(:,:,:,k),[64, 5]);
    %[Uk, Sk, Vk] = svd(Dk);
    %Ca{k} = formatDict(reshape(Uk(:,1),[8, 8]));
    Ca{k} = formatDict(Dhat(:,:,1,k));
    
    Cb{k} = formatDict(squeeze(Dhat(:,col,:,k)));
end
X0a = cell2mov(C0a,gap,val);
X0b = cell2mov(C0b,gap,val);
Xa = cell2mov(Ca,gap,val);
Xb = cell2mov(Cb,gap,val);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([464, 461]);
dy = 0.05;
dx = 0.03;
h = 0.5 * (1 - 3 * dy);
w1 = 8 / 13 * (1 - 3 * dx);
w2 = 5 / 13 * (1 - 3 * dx);

axes('Position',[dx, h + 2 * dy, w1, h]);
imshow(X0a,[0, 1]);
%title('Initial (1st PC)','FontWeight','normal');
title('Initial (1st slice)','FontWeight','normal');

axes('Position',[w1 + 2 * dx, h + 2 * dy, w2, h]);
imshow(X0b,[0, 1]);
title('Initial (y-t)','FontWeight','normal');

axes('Position',[dx, dy, w1, h]);
imshow(Xa,[0, 1]);
%title('Learned (1st PC)','FontWeight','normal');
title('Learned (1st slice)','FontWeight','normal');

axes('Position',[w1 + 2 * dx, dy, w2, h]);
imshow(Xb,[0, 1]);
title('Learned (y-t)','FontWeight','normal');

SetFigFontSize(14);
savePDF(outpath);
%--------------------------------------------------------------------------

%% 0/4 [DATA] load recons

%dr = [5]; p = [0.5];  % gbus per-frame psnrs
dr = [5]; p = [0.7]; % coastguard2 per-frame psnrs
%dr = [5]; p = [0.8];
BASE = '/Users/Brian/Desktop/Online DINO-KAT/TCI-inpaint';
SAVE_FIGURE = true;
FONT_SIZE = 14;
LINEWIDTH = 1.5;


% coastguard2
ONLINE_DLS = 'onlineDls_par4a';
%UNITARY_DLS = 'onlineDls_par4b'; % 1x1
UNITARY_DLS = 'onlineDls_par4g'; % 2x2
FIXED_DLS = 'onlineDls_par4c';
BATCH_DLS = 'dls_par4b';
%CUBIC_2D = 'interp_par3a';
CUBIC_3D = 'interp_par3b';
BASE_OUT_NAME = 'coastguard2';


%{
% coastguard2 (noisy)
ONLINE_DLS = 'onlineDls_npar4a';
%UNITARY_DLS = 'onlineDls_npar4b'; % 1x1
%UNITARY_DLS = 'onlineDls_npar4g'; % 2x2
FIXED_DLS = 'onlineDls_npar4c';
BATCH_DLS = 'dls_npar4b';
%CUBIC_2D = 'interp_npar3a';
CUBIC_3D = 'interp_npar3b';
BASE_OUT_NAME = 'coastguard2_noisy';
%}

%{
% gbus
ONLINE_DLS = 'onlineDls_par6a';
%UNITARY_DLS = 'onlineDls_par6b'; % 1x1
UNITARY_DLS = 'onlineDls_par6g'; % 2x2
FIXED_DLS = 'onlineDls_par6c';
BATCH_DLS = 'dls_par6b';
%CUBIC_2D = 'interp_par5a';
CUBIC_3D = 'interp_par5b';
BASE_OUT_NAME = 'gbus';
%}

%{
% gflower
ONLINE_DLS = 'onlineDls_par8a';
%UNITARY_DLS = 'onlineDls_par8b'; % 1x1
%UNITARY_DLS = 'onlineDls_par8g'; % 2x2
FIXED_DLS = 'onlineDls_par8c';
BATCH_DLS = 'dls_par8b';
%CUBIC_2D = 'interp_par7a';
CUBIC_3D = 'interp_par7b';
BASE_OUT_NAME = 'gflower';
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np = numel(p);
nr = numel(dr);
data = cell(np,0);
xlabels = cell(1,0);
ylabels = cell(1,0);
genPath = @(path) sprintf('%s/%s',BASE,path);

% Load vars
[vars, ~] = eval(sprintf('%s();',ONLINE_DLS));
if ~isfield(vars,'SNR'), vars.SNR = inf; end
na = 0;

% Load (Xture, Y)
for i = 1:np
    [Xtruei, Yi, ~] = loadData(vars.inpath,vars.idim,vars.T,vars.dt,p(i),vars.seed,vars.SNR);
    data{i,na + 1} = Xtruei;
    data{i,na + 2} = Yi;
end
xlabels{na + 1} = 'Truth';
xlabels{na + 2} = 'Corrupted';
na = na + 2;

% Load Online DINO-KAT
for j = 1:nr
    drIdxj = find(ismember(vars.dr,dr(j)));
    for i = 1:np
        pIdxi = find(ismember(vars.p,p(i)));
        idxi = pIdxi + numel(vars.p) * (drIdxj - 1);
        tmp = load(genPath(sprintf('%s/data%d.mat',ONLINE_DLS,idxi)));
        data{i,na + 1} = tmp.Xhat;
    end
    %xlabels{na + 1} = sprintf('Online (r = %d)',dr(j));
    xlabels{na + 1} = 'OnAIR-FD';  % TCI
    na = na + 1;
end

% Load Online DINO-KAT unitary
if exist('UNITARY_DLS','var') && ~isempty(UNITARY_DLS)
    for i = 1:np
        pIdxi = find(ismember(vars.p,p(i)));
        idxi = pIdxi;
        tmp = load(genPath(sprintf('%s/data%d.mat',UNITARY_DLS,idxi)));
        data{i,na + 1} = tmp.Xhat;
    end
    %xlabels{na + 1} = 'Online (unitary)';
    xlabels{na + 1} = 'OnAIR-UD';  % TCI
    na = na + 1;
end

% Load Online DINO-KAT fixed D
if exist('FIXED_DLS','var') && ~isempty(FIXED_DLS)
    for i = 1:np
        pIdxi = find(ismember(vars.p,p(i)));
        idxi = pIdxi;
        tmp = load(genPath(sprintf('%s/data%d.mat',FIXED_DLS,idxi)));
        data{i,na + 1} = tmp.Xhat;
    end
    xlabels{na + 1} = 'Online (DCT)';
    na = na + 1;
end

% Load DINO-KAT
if exist('BATCH_DLS','var') && ~isempty(BATCH_DLS)
    for j = 1:nr
        drIdxj = find(ismember(vars.dr,dr(j)));
        for i = 1:np
            pIdxi = find(ismember(vars.p,p(i)));
            idxi = pIdxi + numel(vars.p) * (drIdxj - 1);
            tmp = load(genPath(sprintf('%s/data%d.mat',BATCH_DLS,idxi)));
            data{i,na + 1} = tmp.Xhat;
        end
        %xlabels{na + 1} = sprintf('Batch (r = %d)',dr(j));
        xlabels{na + 1} = 'Batch';  % TCI
        na = na + 1;
    end
end

% Load interp
if exist('CUBIC_3D','var') && ~isempty(CUBIC_3D)
    for i = 1:np
        pIdxi = find(ismember(vars.p,p(i)));
        idxi = pIdxi;
        tmp = load(genPath(sprintf('%s/data%d.mat',CUBIC_3D,idxi)));
        data{i,na + 1} = tmp.Xhat;
    end
    xlabels{na + 1} = 'Interp (3D)';
    na = na + 1;
end
if exist('CUBIC_2D','var') && ~isempty(CUBIC_2D)
    for i = 1:np
        pIdxi = find(ismember(vars.p,p(i)));
        idxi = pIdxi;
        tmp = load(genPath(sprintf('%s/data%d.mat',CUBIC_2D,idxi)));
        data{i,na + 1} = tmp.Xhat;
    end
    xlabels{na + 1} = 'Interp (2D)';
    na = na + 1;
end

% y-labels
for i = 1:np
    ylabels{i} = sprintf('%.0f%%',100 * p(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute errors
pSNR = nan(np,na);
ppSNR = cell(np,na);
for i = 1:np
    for j = 2:na
        errij = computeErrorMetrics(data{i,j},data{i,1});
        pSNR(i,j) = errij.pSNR;
        ppSNR{i,j} = errij.ppSNR;
    end
end

%% 1/4 [IMAGE] single-frame images

% Knobs
gap = 1;
val = 1;
gamma = 0.9;
%topFrameInds = [1, 3]; % p70, coastguard2
%topFrameInds = [1, 2]; % p80, coastguard2
%topFrameInds = [1, 3]; % p70, coastguard2_noisy
%topFrameInds = [1, 3]; % p80, coastguard2_noisy
%topFrameInds = [2, 4]; % p70, gbus
%topFrameInds = [2, 8]; % p80, gbus
%topFrameInds = [1, 4]; % p70, gflower
%topFrameInds = [1, 5]; % p80, gflower
compIdx = size(ppSNR,2) - 2;

% thesis slides
%topFrameInds = [1, 4];
%compIdx = size(ppSNR,2) - 1;
%FONT_SIZE = 14;

% Find best frames
[~, order] = sort(ppSNR{1,3} - ppSNR{1,compIdx},'descend');
%[~, order] = sort(ppSNR{1,3} + (ppSNR{1,3} - ppSNR{1,compIdx}),'descend');
TOP_5_FRAMES = order(1:10) %#ok

% Generate images
ni = size(data,2);
nf = numel(topFrameInds);
frames = sort(order(topFrameInds));
C = cell(nf,ni);
yLabels = cell(1,nf);
for kk = 1:nf
    fIdx = frames(kk);

    % Generate image
    C(kk,:) = cellfun(@(D) D(:,:,fIdx),data(1,:),'UniformOutput',false);
    improvf = ppSNR{1,3}(fIdx) - ppSNR{1,compIdx}(fIdx) %#ok
    yLabels{kk} = sprintf('Frame %d',fIdx);
end

% Aggregate image
CC = cell2mov(C,gap,val);

% Gamma correction
if gamma ~= 1
    cmg = gammaCorrectColormap(gray(256),gamma);
    CC = colors2im(CC,cmg,[0, 1]);
end

% Plot image
figure();
imshow(CC,[0, 1]);
xLabels = xlabels;
xLabels{2} = sprintf('Corrupted (%.0f%%)',100 * p(1));
addXLabels(xLabels,[],0.013,'top');
addYLabels(yLabels,[],0.003,'left');
SetFigFontSize(FONT_SIZE);

% Resize figure
if ni == 6
    cfigure([863, 312],gcf()); % 4 algos
elseif ni == 7
    cfigure([1028, 312],gcf()); % 5 algos
end

% Save figure
if SAVE_FIGURE
    fStr = sprintf('f%d_',frames);
    name = sprintf('%s_%sp%.0f',BASE_OUT_NAME,fStr,100 * p(1));
    savePDF(name);
end

%% 2/4 [MOVIE] recon + error movies

% Knobs
cm = parula(64);
alpha = 1.0;  % Alpha < 1 to clip large values
gamma = 1.3;  % Gamma < 1 to decrease contrast

% Compute error videos
E = cell(np,na - 2);
for i = 1:np
    for j = 3:na
        E{i,j - 2} = abs(data{i,j} - data{i,1});
    end
end

% Format recon videos
M = cell2mat(data);
M = colors2im(M,gray(256),[0, 1]);

% Format error videos
cmg = gammaCorrectColormap(cm,gamma);
E = cell2mat(E);
E = colors2im(E,cmg,[0, alpha * max(E(:))]);

% For ground truth-only
M1 = cell2mat(data(:,1:2));
M1 = colors2im(M1,gray(256),[0, 1]);

video = cat(1,M,cat(2,M1,E));

% Visualize results
movie = struct(); opts = struct();
movie.video = video;
opts.xlabels = xlabels;
opts.ylabels = ylabels(end:-1:1);
opts.fontSize = FONT_SIZE;
%opts.mag = 0.5; % gbus and gflower
P = PlayMovie(movie,opts);

% Name
name = sprintf('./%s_p%.0f_errors.gif',BASE_OUT_NAME,100 * p(1));
P.SaveGIF(name);
P.Close();

%% 3/4 [MOVIE] recon movies

%{
% thesis
dd = cell(size(data));
for k = 1:numel(dd)
    dd{k} = max(0,min(data{k},1));
end
X = cell2mov(dd,1,1);
%}

% Visualize results
X = max(0,min(cell2mat(data),1));
movie = struct(); opts = struct();
movie.video = X;
opts.xlabels = xlabels;
opts.ylabels = ylabels(end:-1:1);
opts.fontSize = FONT_SIZE;
%opts.mag = 0.5; % gbus and gflower
%P = PlayMovie(movie,opts);

%{
% Name
name = sprintf('./%s_p%.0f_recons.gif',BASE_OUT_NAME,100 * p(1)) %#ok
P.SaveGIF(name);
P.Close();
%}

%{
% Save frames
name = sprintf('%s_p%.0f_recons/frame.png',BASE_OUT_NAME,100 * p(1)) %#ok
saveFrames(X(:,:,1:2:100),name);
%}

%% 4/4 [PLOT] per-frame PSNRs

nc = na - 2;
cm = linspecer(nc);
nt = size(data{1},3);
%ls = {'-','-','-','-','-'};
%ls = {'-','--','-.',':','.-'};
ls = {'x-','v-','^-','d-','s-'};
%ls = {'x','v','^','d','s'};
mn = {11,11,11,11,11}; % number of markers
ms = {6,6,6,6,6}; % markersize
lw = 1.5;

pstr = cell(1,nc);
for i = 1:np
    %cfigure([703, 296]);
    %cfigure([519, 274]);
    cfigure([597, 232]);
    pi = 100 * p(i);
    
    % order
    [~, idx] = sort(cellfun(@mean,ppSNR(i,3:end)),'descend');
    
    % Plot pSNRs
    phndl = nan(1,nc);
    for j = 1:nc
        jj = idx(j);
        %phndl(j) = plot(1:nt,ppSNR{i,2+jj},ls{jj},'Color',cm(jj,:),'Linewidth',LINEWIDTH);
        phndl(j) = line_fewer_markers(1:nt,ppSNR{i,2+jj},mn{jj},ls{jj},'Color',cm(jj,:),'Linewidth',lw,'MarkerSize',ms{jj},'spacing','curve');
        pstr{j} = xlabels{2+jj};
        hold on;
    end
    xlabel('Frame');
    ylabel('PSNR (dB)');
    title('Per-frame PSNR','FontWeight','normal');
    %title(sprintf('Per-frame PSNR [p = %.0f%%]',pi));
    legend(phndl,pstr{:},'Location','EastOutside');
    axis tight; padAxis(); box on;
    SetFigFontSize(FONT_SIZE);
    
    % HACKs
    xTick = get(gca,'XTick'); xTick(1) = 1; set(gca,'XTick',xTick);
    if exist('YLIM','var'), set(gca,'YLim',YLIM); end % HACK
    if exist('XTICK','var'), set(gca,'XTick',XTICK); end % HACK
    
    % Save results
    if SAVE_FIGURE
        name = sprintf('%s_psnrs_p%.0f',BASE_OUT_NAME,pi);
        savePDF(name);
    end
end
