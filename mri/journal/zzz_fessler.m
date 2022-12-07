%% Run FFT simulations on Fessler phantom

% Path
SVTpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\FFT';

% MRI parameters
nc = 8;
SNR = [50];
nt = [2 4 6 8 10 12 16 20 24 30 32 40];

% Loop over SNR
NSNR = length(SNR);
for j = 1:NSNR
    % Perform simulations
    rng(1);
    otimer = tic;
    RunFFTAlgos(nc,nt,SNR(j),SVTpath);
    
    % Display progress
    fprintf('*** Overall trial %i/%i complete - [Time = %.2f] ***\n',j,NSNR,toc(otimer));
end

%% Run SVT simulations on Fessler phantom

% L + S
SVTpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LpS_SVT';
algoStr = 'L+S PGM';
T = 1;

% L & S
%basepath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LaS_SVT';
%algoStr = 'L&S ADMM';
%T = TempFFT(2);

% MRI parameters
nc = 8;
SNR = [50];
nt = [2 4 6 8 10 12 16 20 24 30 32 40];
normA = [22.985 21.143 20.805 20.638 20.543 20.526 20.402 20.378 20.317 20.299 20.344 20.274];

% Regularization parameters
minLambdaL = 0.001;
maxLambdaL = 1;
NlambdaL = 10;
minLambdaS = 0.001;
maxLambdaS = 0.05;
NlambdaS = 10;
%nb = [24 27];

% Loop over SNR
NSNR = length(SNR);
Nt = length(nt);
for j = 1:NSNR
    % Loop over frame counts
    for i = 1:Nt
        % Perform simulations
        rng(1);
        otimer = tic;
        RunSVTAlgos(nc,nt(i),SNR(j),normA(i),T,SVTpath,algoStr,minLambdaL,maxLambdaL,NlambdaL,minLambdaS,maxLambdaS,NlambdaS);
        %RunSVTAlgos(nc,nt(i),SNR(j),normA(i),T,basepath,algoStr,minLambdaL,maxLambdaL,NlambdaL,minLambdaS,maxLambdaS,NlambdaS,nb);
        
        % Display progress
        fprintf('*** Overall trial (%i,%i)/(%i,%i) complete - [Time = %.2f] ***\n',j,i,NSNR,Nt,toc(otimer));
    end
end

%% Run OptShrink simulations on Fessler phantom

% L + S
SVTpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LpS_OptShrink';
algoStr = 'L+S PGM';
T = 1;

% L & S
%basepath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LaS_OptShrink';
%algoStr = 'L&S ADMM';
%T = TempFFT(2);

% MRI parameters
nc = 8;
SNR = [50];
nt = [2 4 6 8 10 12 16 20 24 30 32 40];
normA = [22.985 21.143 20.805 20.638 20.543 20.526 20.402 20.378 20.317 20.299 20.344 20.274];

% Regularization parameters
Nr = 10;
minLambdaS = 0.001;
maxLambdaS = 0.05;
NlambdaS = 10;
%nb = [24 27];

% Loop over SNR
NSNR = length(SNR);
Nt = length(nt);
for j = 1:length(SNR)
    % Loop over frame counts
    for i = 1:length(nt)
        % Perform simulations
        rng(1);
        otimer = tic;
        RunOptShrinkAlgos(nc,nt(i),SNR(j),normA(i),T,SVTpath,algoStr,Nr,minLambdaS,maxLambdaS,NlambdaS);
        %RunOptShrinkAlgos(nc,nt(i),SNR(j),normA(i),T,basepath,algoStr,Nr,minLambdaS,maxLambdaS,NlambdaS,nb);
        
        % Display progress
        fprintf('*** Overall trial (%i,%i)/(%i,%i) complete - [Time = %.2f] ***\n',j,i,NSNR,Nt,toc(otimer));
    end
end

%% Plot mNRMSE curves

% Knobs
FFTpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\FFT_nc_8_SNR_50';
SVTpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LpS_SVT_nc_8_SNR_50';
OPTpath = 'C:\Users\brimoor\Documents\Brian''s Stuff\OptShrink4MRI\LpS_OptShrink_nc_8_SNR_50';
str = {'X','ROI','ROI1','ROI2','ROI3'};
%str = {'X','ROI'};
Nstr = length(str);

% Locate files
basefft = regexprep(fileparts(FFTpath),'\','/');
basesvt = regexprep(fileparts(SVTpath),'\','/');
baseopt = regexprep(fileparts(OPTpath),'\','/');
Dfft = dir(sprintf('%s*_NRMSE.mat',FFTpath));
Dsvt = dir(sprintf('%s*_NRMSE.mat',SVTpath));
Dopt = dir(sprintf('%s*_NRMSE.mat',OPTpath));

% Get NRMSE for FFT
data = load(sprintf('%s/%s',basefft,Dfft(1).name));
nt_fft = data.nt;
NRMSE_fft = cell(Nstr,1);
for j = 1:Nstr
    NRMSE_fft{j} = [data.NRMSE.(str{j})];
end

% Get mNRMSEs for SVT
Nsvt = length(Dsvt);
NRMSE_svt = cell(Nstr,1);
idx_svt = cell(Nstr,1);
nt_svt = zeros(1,Nsvt);
for i = 1:Nsvt
    % Load data
    matObj = matfile(sprintf('%s/%s',basesvt,Dsvt(i).name));
    nt_svt(i) = matObj.nt;
end
[nt_svt inds_svt] = sort(nt_svt);
Dsvt = Dsvt(inds_svt);
for i = 1:Nsvt
    data = load(sprintf('%s/%s',basesvt,Dsvt(i).name));
    for j = 1:Nstr
        [NRMSE_svt{j}(i) idx] = min([data.NRMSE.(str{j})]);
        [idx_svt{j}(1,i) idx_svt{j}(2,i)] = ind2sub(size(data.NRMSE),idx);
    end
end

% Get mNRMSE for OptShrink
Nopt = length(Dopt);
NRMSE_opt = cell(Nstr,1);
idx_opt = cell(Nstr,1);
nt_opt = zeros(1,Nopt);
for i = 1:Nopt
    % Load data
    matObj = matfile(sprintf('%s/%s',baseopt,Dopt(i).name));
    nt_opt(i) = matObj.nt;
end
[nt_opt inds_opt] = sort(nt_opt);
Dopt = Dopt(inds_opt);
for i = 1:Nopt
    % Load data
    data = load(sprintf('%s/%s',baseopt,Dopt(i).name));
    for j = 1:Nstr
        [NRMSE_opt{j}(i) idx] = min([data.NRMSE.(str{j})]);
        [idx_opt{j}(1,i) idx_opt{j}(2,i)] = ind2sub(size(data.NRMSE),idx);
    end
end

% Plot NRMSEs
figure;
[ni nj] = bestSubplotShape(Nstr);
for j = 1:Nstr
    subplot(ni,nj,j);
    p(1) = plot(nt_fft,NRMSE_fft{j},'g-o');
    hold on;
    p(2) = plot(nt_svt,NRMSE_svt{j},'b-o');
    p(3) = plot(nt_opt,NRMSE_opt{j},'r-o');
    xlabel('Number of frames');
    ylabel('NRMSE');
    if (j == 1)
        legend(p,'FFT','SVT','OptShrink');
    end
    title(sprintf('NRMSE [%s]',str{j}));
end
