%% [TEST]  Undersampling

% Knobs
p      = 1 / 8;
SNR    = inf;
seed   = 1;
inpath = 'otazo_full.mat';

% Generate data
[Y, ~, ~, Xtrue, Xfft] = generateCardiacPerfData(p,SNR,seed,inpath);

% Play results
PlayMovie(cat(2,Xtrue,Xfft));

%% [DATA] FFT undersampling sweep

% Knobs
inpath   = 'otazo_full.mat';
%outpath = 'otazo_full_fft.mat';
%outpath = 'otazo_full2_fft.mat';
%p       = [1/64, 1/32, 1/16, 1/8, 1/4, 1/2];
p        = [1/20, 1/16, 1/12, 1/8, 1/4];
seed     = 1;

% Load full data
data = load(inpath);

% Reconstructions
nP    = numel(p);
nS    = numel(seed);
mse   = nan(nP,nS);
count = 0;
for jj = 1:nS
for ii = 1:nP
    % Generate data
    count = count + 1;
    fprintf('\nSimulation %i/%i\n\n',count,numel(mse));
    [~, ~, ~, Xtrue, Xfft] = generateCardiacPerfData(p(ii),inf,seed(jj),data);
    
    % Save results
    mse(ii,jj) = norm(Xtrue(:) - Xfft(:));
end
end

% Save results
if exist('outpath','var') && ~isempty(outpath)
    save(outpath,'mse','p','seed');
end

%% [CHECK] L + S R8

% Knobs
inpath   = 'otazo_R8.mat';
outpath  = '';
%lambdaL = 1.1955;
%lambdaS = 0.01;
r        = 1;
lambdaS  = 0.012;
nIters   = 250;

% Load undersampled data (Y, mask, samp, Xfft, Xtrue)
load(inpath);
A = Emat_xyt(mask,samp);

% L + S
%[Lhat,Shat] = runLpS(A,Y,Xtrue,Xfft,zeros(size(Xfft)),{'svt',lambdaL},lambdaS,nIters);
[Lhat,Shat]  = runLpS(A,Y,Xtrue,Xfft,zeros(size(Xfft)),{'opt',r},lambdaS,nIters);

% Play results
PlayMovie(cat(2,Lhat + Shat,Lhat,Shat));

% Save results
if exist('outpath','var') && ~isempty(outpath)
    save(outpath,'Lhat','Shat','lambdaL','lambdaS','nIters');
end

%% [DATA]  L + S R8 sweep
% mse:    (lambdaL, lambdaS) = (0.525, 0.01)
% tuned:  (lambdaL, lambdaS) = (0.800, 0.01)
% visual: (lambdaL, lambdaS) = (1.196, 0.01)

rng(1);

% Knobs
inpath  = 'otazo_R8.mat';
outpath = 'otazo_R8_lps_mse.mat';
lambdaL = 0.70 * logspace(-0.5,0.5,9);
lambdaS = 0.01 * logspace(-0.5,0.5,9);
nIters  = 250;

% Load undersampled data
load(inpath);
A = Emat_xyt(mask,samp);

% Reconstructions
nL = numel(lambdaL);
nS = numel(lambdaS);
time   = nan(nL,nS,nIters);
cost   = nan(nL,nS,nIters);
mse    = nan(nL,nS,nIters);
delta  = nan(nL,nS,nIters);
mmse   = inf;
Lhat   = nan(size(Xtrue));
Shat   = nan(size(Xtrue));
count  = 0;
for ii = 1:nL
for jj = 1:nS
    % L + S
    count = count + 1;
    fprintf('\nSimulation %i/%i\n\n',count,numel(mse) / nIters);
    [Lh,Sh,it,ti,co,ms,de] = runLpS(A,Y,Xtrue,Xfft,zeros(size(Xfft)),lambdaL(ii),lambdaS(jj),nIters);
    
    % Save results
    time(ii,jj,1:it)  = ti;
    cost(ii,jj,1:it)  = co;
    mse(ii,jj,1:it)   = ms;
    delta(ii,jj,1:it) = de;
    if ms(it) < mmse
        % Update optimal solution
        mmse = ms(it);
        Lhat = Lh;
        Shat = Sh;
    end
end
end
save(outpath,'Lhat','Shat','time','cost','mse','delta','lambdaL','lambdaS','nIters');
load(outpath);

% Optimal parameters
mseFinal     = mse(:,:,nIters);
[~,idx]      = nanmin(mseFinal(:));
[ihat, jhat] = ind2sub([nL nS],idx);
lambdaLopt   = lambdaL(ihat);
lambdaSopt   = lambdaS(jhat);

% Play results
PlayMovie(cat(2,Lhat + Shat,Lhat,Shat));

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
figure;

% lambdaL sensitivity
subplot(2,2,1);
plot(lambdaL,squeeze(mse(:,jhat,nIters)),'b-o');
xlabel('\lambda_L');
ylabel('MSE');
title(sprintf('\\lambda_L = %.3f',lambdaLopt));
axis tight; padAxis();

% lambdaS sensitivity
subplot(2,2,2);
plot(lambdaS,squeeze(mse(ihat,:,nIters)),'b-o');
xlabel('\lambda_S');
ylabel('MSE');
title(sprintf('\\lambda_S = %.3f',lambdaSopt));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,2,3);
semilogy(1:nIters,squeeze(cost(ihat,jhat,:)),'b-');
xlabel('Iteration');
ylabel('MSE');
title('Optimal cost trajectory');
axis tight; padAxis();

% Optimal MSE trajectory
subplot(2,2,4);
semilogy(1:nIters,squeeze(mse(ihat,jhat,:)),'b-');
xlabel('Iteration');
ylabel('MSE');
title('Optimal MSE trajectory');
axis tight; padAxis();

%{
export_fig -pdf -transparent lps_R8
%}
%--------------------------------------------------------------------------

%% [DATA]  L + S undersampling sweep
% mse:    (lambdaL, lambdaS) = (0.525, 0.01)
% tuned:  (lambdaL, lambdaS) = (0.800, 0.01)
% visual: (lambdaL, lambdaS) = (1.196, 0.01)

% Knobs
inpath   = 'otazo_full.mat';
%outpath = 'otazo_full_lps_mse.mat';
%outpath = 'otazo_full_lps_tuned.mat';
%outpath = 'otazo_full_lps_visual.mat';
outpath  = 'otazo_full2_lps_mse.mat';
%p       = [1/64, 1/32, 1/16, 1/8, 1/4, 1/2];
p        = [1/20, 1/16, 1/12, 1/8, 1/4];
seed     = 1;
lambdaL  = 0.525;
lambdaS  = 0.01;
nIters   = 250;

% Load full data
data = load(inpath);

% Reconstructions
nP     = numel(p);
nS     = numel(seed);
time   = nan(nP,nIters,nS);
cost   = nan(nP,nIters,nS);
mse    = nan(nP,nIters,nS);
delta  = nan(nP,nIters,nS);
count  = 0;
for jj = 1:nS
for ii = 1:nP
    % Generate data
    [Y, A, ~, Xtrue, Xfft] = generateCardiacPerfData(p(ii),inf,seed(jj),data);
    
    % L + S
    count = count + 1;
    fprintf('\nSimulation %i/%i\n\n',count,numel(mse) / nIters);
    [~,~,it,ti,co,ms,de] = runLpS(A,Y,Xtrue,Xfft,zeros(size(Xfft)),lambdaL,lambdaS,nIters);
    
    % Save results
    time(ii,1:it,jj)  = ti;
    cost(ii,1:it,jj)  = co;
    mse(ii,1:it,jj)   = ms;
    delta(ii,1:it,jj) = de;
end
end
save(outpath,'time','cost','mse','delta','p','seed','lambdaL','lambdaS','nIters');

%% [PLOT]  L + S undersampling sweep [P1]

% Load data
data_lps_mse_p1  = load('otazo_full2_lps_mse.mat');
PSDENOM          = [20, 16, 12, 8, 4];

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
figure('Position',[705 391 513 421]);
cm = linspecer(1);

phndl(1) = semilogx(data_lps_mse_p1.p,squeeze(nanmean(data_lps_mse_p1.mse(:,data_lps_mse_p1.nIters,:),3)),'-o','Color',cm(1,:)); hold on;
xlabel('Undersampling');
ylabel('MSE');
legend(phndl,sprintf('\\lambda_L = %.3f, \\lambda_S = %.3f',data_lps_mse_p1.lambdaL,data_lps_mse_p1.lambdaS));
axis tight; padAxis();
title('L + S');
set(gca,'XTick',1 ./ PSDENOM);
set(gca,'XTickLabel',arrayfun(@(d)sprintf('1/%d',d),PSDENOM,'UniformOutput',false));

%{
export_fig -pdf -transparent lps_full_p1
%}
%--------------------------------------------------------------------------

%% [PLOT]  L + S undersampling sweep [P3]

% Load data
data_lps_mse    = load('otazo_full_lps_mse.mat');
data_lps_tuned  = load('otazo_full_lps_tuned.mat');
data_lps_visual = load('otazo_full_lps_visual.mat');
PSDENOM         = [64, 32, 16, 8, 4, 2];

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
figure('Position',[705 391 513 421]);
cm = linspecer(3);

phndl(1) = semilogx(data_lps_mse.p,squeeze(nanmean(data_lps_mse.mse(:,data_lps_mse.nIters,:),3)),'-o','Color',cm(1,:)); hold on;
phndl(2) = semilogx(data_lps_tuned.p,squeeze(nanmean(data_lps_tuned.mse(:,data_lps_tuned.nIters,:),3)),'-o','Color',cm(2,:)); hold on;
phndl(3) = semilogx(data_lps_visual.p,squeeze(nanmean(data_lps_visual.mse(:,data_lps_visual.nIters,:),3)),'-o','Color',cm(3,:)); hold on;
xlabel('Undersampling');
ylabel('MSE');
legend(phndl,sprintf('\\lambda_L = %.3f, \\lambda_S = %.3f',data_lps_mse.lambdaL,data_lps_mse.lambdaS), ...
             sprintf('\\lambda_L = %.3f, \\lambda_S = %.3f',data_lps_tuned.lambdaL,data_lps_tuned.lambdaS), ...
             sprintf('\\lambda_L = %.3f, \\lambda_S = %.3f',data_lps_visual.lambdaL,data_lps_visual.lambdaS));
axis tight; padAxis();
title('L + S');
set(gca,'XTick',1 ./ PSDENOM);
set(gca,'XTickLabel',arrayfun(@(d)sprintf('1/%d',d),PSDENOM,'UniformOutput',false));

%{
export_fig -pdf -transparent lps_full
%}
%--------------------------------------------------------------------------

%% [CHECK] LADS R8

rng(1);

% Knobs
inpath    = 'otazo_R8.mat';
%initpath = 'otazo_R8_lps_mse.mat';
%initpath = 'otazo_R8_lps_tuned.mat';
initpath  = 'otazo_R8_lps_visual.mat';
%initpath = 'otazo_R8_lps_visual20.mat';
outpath   = 'otazo_R8_lads_p1_mse250.mat';
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = 0.03;
lambdaB0  = 0.03;
p         = 1;
nIters    = 50;

% Load undersampled data
load(inpath);
A = Emat_xyt(mask,samp);

% Load L + S initialization
init = load(initpath);
L0   = init.Lhat;
S0   = init.Shat;

% LADS
lB = logspace(log10(lambdaB0),log10(lambdaB),nIters);
[Lhat,Shat,Dhat,~,ite,time,cost,mse,delta,sparsity] = runLADS(A,Y,Xtrue,L0,S0,lambdaL,lambdaS,lB,p,nIters); %#ok

% Save results
save(outpath,'Lhat','Shat','Dhat','time','cost','mse','delta','sparsity', ...
             'lambdaL','lambdaS','lambdaB','lambdaB0','p','nIters');
fprintf('DONE\n');

%% [PLOT]  LADS R8 sweep

% Knobs
inpath  = 'otazo_R8_lads_mse250.mat';
%inpath = 'otazo_R8_lads_mse250a.mat';
%inpath = 'otazo_R8_lads_mse250b.mat';
%inpath = 'otazo_R8_lads_mse20.mat';

% Load data
load(inpath);

% Compute optimal parameters
nL = numel(lambdaL);
nS = numel(lambdaS);
nB = numel(lambdaB);
mseFinal   = mse(:,:,:,nIters);
[mmse,idx] = nanmin(mseFinal(:));
[ihat, jhat, khat] = ind2sub([nL nS nB],idx);
lambdaLopt = lambdaL(ihat);
lambdaSopt = lambdaS(jhat);
lambdaBopt = lambdaB(khat);

%{
% Play results
PlayMovie(cat(2,Lhat + Shat,Lhat,Shat));
%}

%--------------------------------------------------------------------------
% Plot sensitivities
%--------------------------------------------------------------------------
figure('Position',[489 335 945 532]);

% lambdaL
subplot(2,3,1);
plot(lambdaL,squeeze(mse(:,jhat,khat,nIters)),'b-o');
xlabel('\lambda_L');
ylabel('MSE');
title(sprintf('\\lambda_L = %.3f',lambdaLopt));
axis tight; padAxis();

% lambdaS
subplot(2,3,2);
plot(lambdaS,squeeze(mse(ihat,:,khat,nIters)),'b-o');
xlabel('\lambda_S');
ylabel('MSE');
title(sprintf('\\lambda_S = %.3f',lambdaSopt));
axis tight; padAxis();

% lambdaB
subplot(2,3,3);
plot(lambdaB,squeeze(mse(ihat,jhat,:,nIters)),'b-o');
xlabel('\lambda_B');
ylabel('MSE');
title(sprintf('\\lambda_B = %.3f',lambdaBopt));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,3,4);
semilogy(1:nIters,squeeze(cost(ihat,jhat,khat,:)),'b-');
xlabel('Iteration');
ylabel('cost');
title('Optimal cost trajectory');
axis tight; padAxis();

% Optimal MSE trajectory
subplot(2,3,5);
semilogy(1:nIters,squeeze(mse(ihat,jhat,khat,:)),'b-');
xlabel('Iteration');
ylabel('MSE');
title(sprintf('Optimal MSE (%.3f) trajectory',mmse));
axis tight; padAxis();

% Optimal sparsity trajectory
subplot(2,3,6);
semilogy(1:nIters,squeeze(sparsity(ihat,jhat,khat,:)),'b-');
xlabel('Iteration');
ylabel('Sparsity %');
title('Optimal sparsity trajectory');
axis tight; padAxis();

%{
export_fig -pdf -transparent lads_R8_250
export_fig -pdf -transparent lads_R8_250a
export_fig -pdf -transparent lads_R8_250b
export_fig -pdf -transparent lads_R8_20
%}
%--------------------------------------------------------------------------

%% [CHECK] LADS undersampling sweep
% 250:    (lambdaL, lambdaS, lambdaB, lambdaB0, p) = (0.5, 0.01, 0.03, 0.03, 3)
% 250b:   (lambdaL, lambdaS, lambdaB, lambdaB0, p) = (0.5, 0.01, 0.04, 0.04, 3)
% 20:     (lambdaL, lambdaS, lambdaB, lambdaB0, p) = (0.5, 0.01, 0.03, 0.06, 3)

% Knobs
inpath    = 'otazo_full.mat';
initpath  = 'otazo_R8_lps_visual.mat';      % 250
outpath   = 'otazo_full_lads250.mat';       % 250
%initpath = 'otazo_R8_lps_mse.mat';         % 250b
%outpath  = 'otazo_full_lads250b.mat';      % 250b
%initpath = 'otazo_R8_lps_visual20.mat';    % 20
%outpath  = 'otazo_full_lads20.mat';        % 20
ps        = [1/64, 1/32, 1/16, 1/8, 1/4, 1/2];
seed      = 1;
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = 0.03;
lambdaB0  = 0.03;
p         = 3;
nIters    = 50;

% Load full data
data = load(inpath);

% Load L + S initialization (Lhat, Shat, lambdaL, lambdaS, nIters)
init = load(initpath);
L0   = init.Lhat;
S0   = init.Shat;

% Reconstructions
nP       = numel(ps);
nS       = numel(seed);
time     = nan(nP,nIters,nS);
cost     = nan(nP,nIters,nS);
mse      = nan(nP,nIters,nS);
delta    = nan(nP,nIters,nS);
sparsity = nan(nP,nIters,nS);
count    = 0;
for jj = 1:nS
for ii = 1:nP
    % Generate data
    [Y, A, ~, Xtrue, Xfft] = generateCardiacPerfData(ps(ii),inf,seed(jj),data);
    
    % LADS
    count = count + 1;
    fprintf('\nSimulation %i/%i\n\n',count,numel(mse) / nIters);
    lB = logspace(log10(lambdaB0),log10(lambdaB),nIters);
    [Lh,Sh,~,~,it,ti,co,ms,de,sp] = runLADS(A,Y,Xtrue,L0,S0,lambdaL,lambdaS,lB,p,nIters);
    
    % Save results
    time(ii,1:it,jj)     = ti;
    cost(ii,1:it,jj)     = co;
    mse(ii,1:it,jj)      = ms;
    delta(ii,1:it,jj)    = de;
    sparsity(ii,1:it,jj) = sp;
end
end
save(outpath,'time','cost','mse','delta','sparsity','ps','seed', ...
             'lambdaL','lambdaS','lambdaB','lambdaB0','p','nIters');

%% [PLOT]  LADS undersampling sweep [P1]

% Load data
data_lads_p1_250 = load('otazo_full_lads_p1_250.mat');
PSDENOM          = [20, 16, 12, 8, 4];

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
figure('Position',[705 391 513 421]);
cm = linspecer(1);

phndl(1) = semilogx(data_lads_p1_250.ps,squeeze(nanmean(data_lads_p1_250.mse(:,data_lads_p1_250.nIters,:),3)),'-o','Color',cm(1,:)); hold on;
xlabel('Undersampling');
ylabel('MSE');
legend(phndl,sprintf('\\lambda_L = %.3f, \\lambda_S = %.3f, \\lambda_B = %.3f, p = %d',data_lads_p1_250.lambdaL,data_lads_p1_250.lambdaS,data_lads_p1_250.lambdaB,data_lads_p1_250.p));
axis tight; padAxis();
title('LADS');
set(gca,'XTick',1 ./ PSDENOM);
set(gca,'XTickLabel',arrayfun(@(d)sprintf('1/%d',d),PSDENOM,'UniformOutput',false));

%{
export_fig -pdf -transparent lads_full_p1
%}
%--------------------------------------------------------------------------

%% [PLOT]  LADS undersampling sweep [P3]

% Load data
%data_lads250  = load('otazo_full_lads250.mat');
%data_lads250b = load('otazo_full_lads250b.mat');
%data_lads20   = load('otazo_full_lads20.mat');
data_lads250   = load('otazo_full_lads_2_250.mat');
data_lads250b  = load('otazo_full_lads_2_250b.mat');
data_lads20    = load('otazo_full_lads_2_20.mat');
PSDENOM        = [64, 32, 16, 8, 4, 2];

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
figure('Position',[705 391 513 421]);
cm = linspecer(3);

phndl(1) = semilogx(data_lads250.ps,squeeze(nanmean(data_lads250.mse(:,data_lads250.nIters,:),3)),'-o','Color',cm(1,:)); hold on;
phndl(2) = semilogx(data_lads250b.ps,squeeze(nanmean(data_lads250b.mse(:,data_lads250b.nIters,:),3)),'-o','Color',cm(2,:)); hold on;
phndl(3) = semilogx(data_lads20.ps,squeeze(nanmean(data_lads20.mse(:,data_lads20.nIters,:),3)),'-o','Color',cm(3,:)); hold on;
xlabel('Undersampling');
ylabel('MSE');
legend(phndl,sprintf('\\lambda_L = %.3f, \\lambda_S = %.3f, \\lambda_B = %.3f, p = %d',data_lads250.lambdaL,data_lads250.lambdaS,data_lads250.lambdaB,data_lads250.p), ...
             sprintf('\\lambda_L = %.3f, \\lambda_S = %.3f, \\lambda_B = %.3f, p = %d',data_lads250b.lambdaL,data_lads250b.lambdaS,data_lads250b.lambdaB,data_lads250b.p), ...
             sprintf('\\lambda_L = %.3f, \\lambda_S = %.3f, \\lambda_B = %.3f, p = %d',data_lads20.lambdaL,data_lads20.lambdaS,data_lads20.lambdaB,data_lads20.p));
axis tight; padAxis();
title('LADS');
set(gca,'XTick',1 ./ PSDENOM);
set(gca,'XTickLabel',arrayfun(@(d)sprintf('1/%d',d),PSDENOM,'UniformOutput',false));

%{
export_fig -pdf -transparent lads_full
export_fig -pdf -transparent lads_full_2
%}
%--------------------------------------------------------------------------

%% [CHECK] LADS p sweep
% 250:    (lambdaL, lambdaS, lambdaB, lambdaB0, p) = (0.5, 0.01, 0.03, 0.03, 3)
% 250b:   (lambdaL, lambdaS, lambdaB, lambdaB0, p) = (0.5, 0.01, 0.04, 0.04, 3)
% 20:     (lambdaL, lambdaS, lambdaB, lambdaB0, p) = (0.5, 0.01, 0.03, 0.06, 3)

rng(1);

% Knobs
inpath    = 'otazo_R8.mat';
initpath  = 'otazo_R8_lps_visual.mat';      % 250
outpath   = 'otazo_p_lads250.mat';          % 250
%initpath = 'otazo_R8_lps_mse.mat';         % 250b
%outpath  = 'otazo_p_lads250b.mat';         % 250b
%initpath = 'otazo_R8_lps_visual20.mat';    % 20
%outpath  = 'otazo_p_lads20.mat';           % 20
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = 0.03;
lambdaB0  = 0.03;
p         = [1, 2, 3, 4, 5];
nIters    = 50;

% Load undersampled data
load(inpath);
A = Emat_xyt(mask,samp);

% Load L + S initialization (Lhat, Shat, lambdaL, lambdaS, nIters)
init = load(initpath);
L0   = init.Lhat;
S0   = init.Shat;

% Reconstructions
nP       = numel(p);
time     = nan(nP,nIters);
cost     = nan(nP,nIters);
mse      = nan(nP,nIters);
delta    = nan(nP,nIters);
sparsity = nan(nP,nIters);
count  = 0;
for ii = 1:nP
    % LADS
    count = count + 1;
    fprintf('\nSimulation %i/%i\n\n',count,numel(mse) / nIters);
    lB = logspace(log10(lambdaB0),log10(lambdaB),nIters);
    [Lh,Sh,~,~,it,ti,co,ms,de,sp] = runLADS(A,Y,Xtrue,L0,S0,lambdaL,lambdaS,lB,p(ii),nIters);
    
    % Save results
    time(ii,1:it)     = ti;
    cost(ii,1:it)     = co;
    mse(ii,1:it)      = ms;
    delta(ii,1:it)    = de;
    sparsity(ii,1:it) = sp;
end
save(outpath,'time','cost','mse','delta','sparsity', ...
             'lambdaL','lambdaS','lambdaB','lambdaB0','p','nIters');

%% [PLOT]  LADS p sweep

% Load data
data_lads250  = load('otazo_p_lads250.mat');
data_lads250b = load('otazo_p_lads250b.mat');
data_lads20   = load('otazo_p_lads20.mat');

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
figure('Position',[742 437 439 329]);
cm = linspecer(3);

phndl(1) = plot(data_lads250.p,squeeze(data_lads250.mse(:,data_lads250.nIters)),'-o','Color',cm(1,:)); hold on;
phndl(2) = plot(data_lads250b.p,squeeze(data_lads250b.mse(:,data_lads250b.nIters)),'-o','Color',cm(2,:)); hold on;
phndl(3) = plot(data_lads20.p,squeeze(data_lads20.mse(:,data_lads20.nIters)),'-o','Color',cm(3,:)); hold on;
xlabel('Dictionary atom rank');
ylabel('MSE');
legend(phndl,sprintf('\\lambda_L = %.3f, \\lambda_S = %.3f, \\lambda_B = %.3f',data_lads250.lambdaL,data_lads250.lambdaS,data_lads250.lambdaB), ...
             sprintf('\\lambda_L = %.3f, \\lambda_S = %.3f, \\lambda_B = %.3f',data_lads250b.lambdaL,data_lads250b.lambdaS,data_lads250b.lambdaB), ...
             sprintf('\\lambda_L = %.3f, \\lambda_S = %.3f, \\lambda_B = %.3f',data_lads20.lambdaL,data_lads20.lambdaS,data_lads20.lambdaB));
axis tight; padAxis();
title('LADS');

%{
export_fig -pdf -transparent lads_p
%}
%--------------------------------------------------------------------------

%% [FIGS1] UNDERSAMPLING [P1]
% ps = [1/20, 1/16, 1/12, 1/8, 1/4]

% Load data
data_gt   = load('otazo_R8.mat');
data_fft  = load('otazo_full2_fft.mat');
data_lps  = load('otazo_full2_lps_mse.mat');
data_lads = load('otazo_full_lads_p1_250.mat');

% Extract MSEs
mse_fft  = nanmean(data_fft.mse,2);
mse_lps  = squeeze(nanmean(data_lps.mse(:,data_lps.nIters,:),3));
mse_lads = squeeze(nanmean(data_lads.mse(:,data_lads.nIters,:),3));

% Compute NRMSEs
ps       = data_fft.p(:)' %#ok
%ps      = data_lps.p(:)' %#ok
%ps      = data_lads.ps(:)' %#ok
normXtrue  = norm(data_gt.Xtrue(:));
nrmse_fft  = 100 * (mse_fft(:)' / normXtrue) %#ok
nrmse_lps  = 100 * (mse_lps(:)' / normXtrue) %#ok
nrmse_lads = 100 * (mse_lads(:)' / normXtrue) %#ok
gainDB     = 20 * log10(nrmse_lps ./ nrmse_lads) %#ok

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
figure('Position',[786 454 351 295]);
cm = linspecer(2);

% L + S
phndl(1) = semilogx(ps,mse_lps,'-o','Color',cm(1,:)); hold on;

% LADS
phndl(2) = semilogx(ps,mse_lads,'-o','Color',cm(2,:)); hold on;

xlabel('Undersampling');
ylabel('MSE');
legend(phndl,'L + S','LADS');
axis tight; padAxis();
set(gca,'XTick',ps);
set(gca,'XTickLabel',arrayfun(@(ps)sprintf('1/%d',round(1 ./ ps)),ps,'UniformOutput',false));

%{
export_fig -pdf -transparent undersampling_p1
%}
%--------------------------------------------------------------------------

%% [FIGS1] UNDERSAMPLING [P3]
% ps = [1/64, 1/32, 1/16, 1/8, 1/4, 1/2]

% Load data
data_gt    = load('otazo_R8.mat');
data_fft   = load('otazo_full_fft.mat');
data_lps   = load('otazo_full_lps_mse.mat');
%data_lps  = load('otazo_full_lps_tuned.mat');
%data_lps  = load('otazo_full_lps_visual.mat');
%data_lads = load('otazo_full_lads250.mat');
%data_lads = load('otazo_full_lads250b.mat');
%data_lads = load('otazo_full_lads20.mat');
data_lads  = load('otazo_full_lads_2_250.mat');
%data_lads = load('otazo_full_lads_2_250b.mat');
%data_lads = load('otazo_full_lads_2_20.mat');

% Extract MSEs
mse_fft  = nanmean(data_fft.mse,2);
mse_lps  = squeeze(nanmean(data_lps.mse(:,data_lps.nIters,:),3));
mse_lads = squeeze(nanmean(data_lads.mse(:,data_lads.nIters,:),3));

% Compute NRMSEs
ps       = data_fft.p(:)' %#ok
%ps      = data_lps.p(:)' %#ok
%ps      = data_lads.ps(:)' %#ok
normXtrue  = norm(data_gt.Xtrue(:));
nrmse_fft  = 100 * (mse_fft(:)' / normXtrue) %#ok
nrmse_lps  = 100 * (mse_lps(:)' / normXtrue) %#ok
nrmse_lads = 100 * (mse_lads(:)' / normXtrue) %#ok
gainDB     = 20 * log10(nrmse_lps ./ nrmse_lads) %#ok

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
figure('Position',[786 454 351 295]);
cm = linspecer(2);

% L + S
phndl(1) = semilogx(ps,mse_lps,'-o','Color',cm(1,:)); hold on;

% LADS
phndl(2) = semilogx(ps,mse_lads,'-o','Color',cm(2,:)); hold on;

xlabel('Undersampling');
ylabel('MSE');
legend(phndl,'L + S','LADS');
axis tight; padAxis();
set(gca,'XTick',ps);
set(gca,'XTickLabel',arrayfun(@(ps)sprintf('1/%d',round(1 ./ ps)),ps,'UniformOutput',false));

%{
export_fig -pdf -transparent undersampling
export_fig -pdf -transparent undersampling_2
%}
%--------------------------------------------------------------------------

%% [FIGS2] PICTURES 

% Knobs
gtpath    = 'otazo_R8.mat';                 % Xtrue
lpspath   = 'otazo_R8_lps_mse.mat';         % L + S MMSE
%lpspath  = 'otazo_R8_lps_visual.mat';      % L + S visual
%ladspath = 'otazo_R8_lads_mse250.mat';     % LADS MMSE (L+S: 250 visual)
%ladspath = 'otazo_R8_lads_mse20.mat';      % LADS MMSE (L+S: 20 visual)
ladspath  = 'otazo_R8_lads_p1_mse250.mat';  % LADS MMSE p = 1 (L+S: 250 visual)
f1        = 7;
f2        = 13;

% Load data
gtdata   = load(gtpath);
lpsdata  = load(lpspath);
ladsdata = load(ladspath);
Xtrue    = gtdata.Xtrue;
Llps     = lpsdata.Lhat;
Slps     = lpsdata.Shat;
Xlps     = lpsdata.Lhat + lpsdata.Shat;
Llads    = ladsdata.Lhat;
Slads    = ladsdata.Shat;
Xlads    = ladsdata.Lhat + ladsdata.Shat;

%{
% Check abs() vs. non-abs() MSEs
mse_lps  = [norm(vec(Xtrue - Xlps)),  norm(vec(abs(Xtrue) - abs(Xlps)))] %#ok
mse_lads = [norm(vec(Xtrue - Xlads)), norm(vec(abs(Xtrue) - abs(Xlads)))] %#ok
gain_db  = 20 * log10(mse_lps ./ mse_lads) %#ok
%}

% Normalize data
M = max(abs(Xtrue(:)));
Xtrue = abs(Xtrue) / M; Xtrue(Xtrue > 1) = 1;
Xlps  = abs(Xlps)  / M;  Xlps(Xlps > 1)  = 1;
Llps  = abs(Llps)  / M;  Llps(Llps > 1)  = 1;
Slps  = abs(Slps)  / M;  Slps(Slps > 1)  = 1;
Xlads = abs(Xlads) / M; Xlads(Xlads > 1) = 1;
Llads = abs(Llads) / M; Llads(Llads > 1) = 1;
Slads = abs(Slads) / M; Slads(Slads > 1) = 1;

% Reconstruction video
Xt = cat(1,zeros(64,128,40),Xtrue,zeros(64,128,40));
M  = cat(1,cat(2,Xlps,Llps,Slps), ...
           cat(2,Xlads,Llads,Slads));
%PlayMovie(cat(2,Xt,M));

% Error video
Elps  = abs(abs(Xtrue) - abs(Xlps));
Elads = abs(abs(Xtrue) - abs(Xlads));
%PlayMovie(cat(2,Elps,Elads));

%{
%--------------------------------------------------------------------------
% Per-frame errors
%-------------------------------------------------------------------------
% Compute stats
errlps  = nanmean(reshape(Elps,[],size(Elps,3)).^2,1);
errlads = nanmean(reshape(Elads,[],size(Elads,3)).^2,1);
improve = (errlps - errlads);
[~,bestFrames] = sort(improve,'descend') %#ok

figure;

subplot(1,2,1);
phndl(1) = plot(errlps,'b-o'); hold on;
phndl(2) = plot(errlads,'r-o'); hold on;
legend(phndl,'L + S','LADS');
xlabel('frame');
ylabel('MSE');

subplot(1,2,2);
plot(improve,'b-o');
xlabel('frame');
ylabel('Improvement');
%--------------------------------------------------------------------------
%}

%--------------------------------------------------------------------------
% Pictures
%--------------------------------------------------------------------------
xlabels  = {'Reference','','L','S'};
ylabels  = {sprintf('F%d',f1),sprintf('F%d',f2)};
gap      = 2;
fontSize = 12;

% L + S
C1 = {Xtrue(:,:,f1), Xlps(:,:,f1), Llps(:,:,f1), Slps(:,:,f1);
      Xtrue(:,:,f2), Xlps(:,:,f2), Llps(:,:,f2), Slps(:,:,f2)};
xlabels{2} = 'L + S  [11]';
labelImage(cell2mov(C1,1),xlabels,ylabels,gap,fontSize);
set(gcf,'Position',[622 430 679 343]);
%{
export_fig -pdf -transparent lps_frames_msec
export_fig -pdf -transparent lps_frames_visual
%}

% LADS
C2 = {Xtrue(:,:,f1), Xlads(:,:,f1), Llads(:,:,f1), Slads(:,:,f1);
      Xtrue(:,:,f2), Xlads(:,:,f2), Llads(:,:,f2), Slads(:,:,f2)};
xlabels{2} = 'LASSI';
labelImage(cell2mov(C2,1),xlabels,ylabels,gap,fontSize);
set(gcf,'Position',[622 430 679 343]);
%{
export_fig -pdf -transparent lads_frames250
export_fig -pdf -transparent lads_frames20
export_fig -pdf -transparent lads_frames_p1_250c
%}
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Zoomed pictures
% NOTE: must add arrows manually
%--------------------------------------------------------------------------
xlabels  = {'Reference','L + S  [11]','LASSI'};
ylabels  = {sprintf('F%d',f2)};
rows     = 33:96; % Rows
cols     = 33:96; % Columns
gap      = 2;
fontSize = 12;

% Zoomed images
C3 = {Xtrue(rows,cols,f2), Xlps(rows,cols,f2), Xlads(rows,cols,f2)};
labelImage(cell2mov(C3,1),xlabels,ylabels,gap,fontSize);
set(gcf,'Position',[615 430 692 343]);
%{
export_fig -pdf -transparent zoom_frames_mse
export_fig -pdf -transparent zoom_frames_visual
%}
%--------------------------------------------------------------------------

%% [PLOT]  L + S [OPT] 

% Knobs
inpath  = 'otazo_R8_lps_opt.mat';

% Load data
load(inpath);

% Compute optimal parameters
nR = numel(r);
nS = numel(lambdaS);
mseFinal     = mse(:,:,nIters);
[mmse,idx]   = nanmin(mseFinal(:));
[ihat, jhat] = ind2sub([nR nS],idx);
ropt       = r(ihat);
lambdaSopt = lambdaS(jhat);

% Play results
PlayMovie(cat(2,Lhat + Shat,Lhat,Shat));

%--------------------------------------------------------------------------
% Plot sensitivities
%--------------------------------------------------------------------------
figure;

% lambdaL
subplot(1,3,1);
plot(r,squeeze(mse(:,jhat,nIters)),'b-o');
xlabel('r');
ylabel('MSE');
title(sprintf('r = %d',ropt));
axis tight; padAxis();

% lambdaS
subplot(1,3,2);
plot(lambdaS,squeeze(mse(ihat,:,nIters)),'b-o');
xlabel('\lambda_S');
ylabel('MSE');
title(sprintf('\\lambda_S = %.3f',lambdaSopt));
axis tight; padAxis();

% Optimal MSE trajectory
subplot(1,3,3);
semilogy(1:nIters,squeeze(mse(ihat,jhat,:)),'b-');
xlabel('Iteration');
ylabel('MSE');
title(sprintf('Optimal MSE (%.3f) trajectory',mmse));
axis tight; padAxis();

%{
export_fig -pdf -transparent lps_R8_opt
%}
%--------------------------------------------------------------------------

%% [CHECK] LPS_IST vs. robustPCA

rng(1);

% Knobs
inpath  = 'otazo_R8.mat';
lambdaL = 1.1955;
lambdaS = 0.01;
nIters  = 20;

% Load undersampled data
load(inpath);
[ny, nx, nt] = size(Xtrue);

% LPS
param.E        = Emat_xyt(mask,samp);
param.d        = Y;
param.T        = TempFFT(3);
param.L0       = Xfft;
param.S0       = zeros(size(Xfft));
param.Lfull    = Xtrue;
param.nite     = nIters;
param.tol      = -1;
param.lambda_L = lambdaL;
param.lambda_S = lambdaS;
[~,~,~,~,cost1,mse1,~] = lps_ist2(param);
nrmse1 = mse1 / norm(Xtrue(:));

% Robust PCA
opts.A      = Emat_xyt(mask,samp,[ny, nx, nt]);
opts.T      = TempFFT(3,[ny, nx, nt]);
opts.nIters = nIters;
opts.L0     = reshape(Xfft,[],nt);
opts.S0     = zeros(ny * nx,nt);
opts.Xtrue  = reshape(Xtrue,[],nt);
opts.accel  = false;
opts.tau    = 1;
opts.flag   = 1;
[~,~,stats2] = robustPCA(Y,lambdaL,lambdaS,opts);
cost2  = stats2.cost;
nrmse2 = stats2.nrmse;

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure();
cm = linspecer(2);

% Cost
subplot(1,2,1);
plot(1:nIters,cost1,'o','Color',cm(1,:)); hold on;
plot(1:nIters,cost2,'-','Color',cm(2,:));
xlabel('Iteration');
title('Cost');
axis tight; padAxis();

% NRMSE
subplot(1,2,2);
phndl = zeros(1,2);
phndl(1) = plot(1:nIters,nrmse1,'o','Color',cm(1,:)); hold on;
phndl(2) = plot(1:nIters,nrmse2,'-','Color',cm(2,:));
xlabel('Iteration');
title('NRMSE');
legend(phndl,'L + S','Robust PCA');
axis tight; padAxis();
%--------------------------------------------------------------------------

%% [CHECK] LADS vs. DRPCA

rng(1);

% Knobs
inpath    = 'otazo_R8.mat';
initpath  = 'otazo_R8_lps_visual.mat';
nIters    = 3;
lambdaL   = 0.5;
lambdaS   = 0.01;
lambdaB   = logspace(log10(0.06),log10(0.03),nIters);
p         = 1;

% Load undersampled data
load(inpath);
[ny, nx, nt] = size(Xtrue);

% Load L + S initialization
init = load(initpath);
L0   = init.Lhat;
S0   = init.Shat;

% LADS
param.E         = Emat_xyt(mask,samp);
param.d         = Y;
param.L0        = L0;
param.S0        = S0;
param.Lfull     = Xtrue;
param.nite      = nIters;
param.niteDL    = 1;
param.niteLS    = 5;
param.tol       = -1;
param.patchsize = [8, 8, 5];
param.slidedis  = [2, 2];
param.D         = dctmtx(prod(param.patchsize));
param.lambda_L  = lambdaL;
param.lambda_S  = lambdaS;
param.lambda_B  = lambdaB;
param.p         = p;
[~,~,~,~,~,time1,cost1,mse1,delta1,sparsity1] = lads_ist2(param);
nrmse1 = mse1 / norm(Xtrue(:));

% DRPCA
opts.A        = Emat_xyt(mask,samp,[ny, nx, nt]);
opts.sdim     = [ny, nx, nt];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 2];
opts.type     = 'hard';
opts.dr       = p;
opts.ddim     = [64, 5];
opts.nIters   = nIters;
opts.nItersDB = 1;
opts.nItersLS = 5;
opts.L0       = reshape(L0,[],nt);
opts.S0       = reshape(S0,[],nt);
opts.D0       = dctmtx(prod(opts.pdim));
opts.Xtrue    = reshape(Xtrue,[],nt);
opts.accel    = false;
opts.tau      = 1;
opts.flag     = 1;
[~,~,~,~,stats2] = drpca(Y,lambdaL,2 * lambdaS,lambdaB,opts);
time2     = stats2.time;
cost2     = stats2.cost;
nrmse2    = stats2.nrmse;
sparsity2 = stats2.sparsity;

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([950, 250]);
cm = linspecer(2);

% Cost
subplot(1,3,1);
plot(1:nIters,cost1,'o','Color',cm(1,:)); hold on;
plot(1:nIters,cost2,'-','Color',cm(2,:));
xlabel('Iteration');
title('Cost');
axis tight; padAxis();

% NRMSE
subplot(1,3,2);
phndl = zeros(1,2);
phndl(1) = plot(1:nIters,nrmse1,'o','Color',cm(1,:)); hold on;
phndl(2) = plot(1:nIters,nrmse2,'-','Color',cm(2,:));
xlabel('Iteration');
title('NRMSE');
legend(phndl,'LASSI','DRPCA');
axis tight; padAxis();

% Sparsity
subplot(1,3,3);
plot(1:nIters,sparsity1,'o','Color',cm(1,:)); hold on;
plot(1:nIters,sparsity2,'-','Color',cm(2,:));
xlabel('Iteration');
title('Sparsity');
axis tight; padAxis();
%--------------------------------------------------------------------------
