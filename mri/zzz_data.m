%
% Demo the performance of ground truth data and its reconstructions via
% IST [Otazo] and SLRN [Proposed]
%

%% 1a. Synthetic PINCAT phantom (not low-rank...)

% Knobs
fps = 8; % frames per second
res = 512; % # pixels width
title = 'PINCAT Phantom';

% Show video
load('phantom.mat');
M = abs(reshape(Mtrue,[nx ny nt]));
DataViewer(M,fps,res,title);

%% 1b. (from kt-slr) Synthetic PINCAT phantom (not low-rank...)

% Knobs
fps = 8; % frames per second
res = 512; % # pixels width
title = 'PINCAT Phantom';

% Show video
load('./ktslr/aperiodic_pincat.mat');
DataViewer(new,fps,res,title);

%% 1c. (from kt-slr) In vivo perfusion (full data scan)

% Knobs
fps = 8; % frames per second
res = 570; % # pixels width
title = 'In Vivo Perfusion [Full Scan]';

% Show video
load('./ktslr/invivo_perfusion.mat');
DataViewer(x,fps,res,title);

%% 2a. Cardiac perfusion #1 demo (very low-rank)

% Knobs
fps = 8; % frames per second
res = 500; % # pixels width
title = 'Cardiac Perfusion #1 [Proposed]';

% Show video
load('./figures/LS_cardiac_SLRN.mat');
M = abs(Ltrue + Strue);
DataViewer(M,fps,res,title);

%% 2b. Cardiac perfusion #1 demo (very low-rank)

% Knobs
fps = 8; % frames per second
res = 500; % # pixels width
title = {'Cardiac Perfusion #1 [Proposed]';
         'Cardiac Perfusion #1 [IST]'};

% Show video
load('./figures/LS_cardiac_SLRN.mat');
M1 = abs(Ltrue + Strue);
load('./figures/LS_cardiac_ist.mat');
M2 = abs(Ltrue + Strue);
DataViewer({M1;M2},fps,res,title);

%% 3a. Cardiac perfusion #2 demo (modestly low-rank)

% Knobs
fps = 8; % frames per second
res = 500; % # pixels width
title = 'Cardiac Perfusion #2 [Proposed]';

% Show video
load('./figures/LS_invivo_SLRN.mat');
M = abs(Ltrue + Strue);
DataViewer(M,fps,res,title);

%% 3b. Cardiac perfusion #2 demo (modestly low-rank)

% Knobs
fps = 8; % frames per second
res = 500; % # pixels width
title = {'Cardiac Perfusion #2 [Proposed]';
         'Cardiac Perfusion #2 [IST]'};

% Show video
load('./figures/LS_invivo_SLRN.mat');
M1 = abs(Ltrue + Strue);
load('./figures/LS_invivo_ist.mat');
M2 = abs(Ltrue + Strue);
DataViewer({M1;M2},fps,res,title);

%% 4a. Simuated Block-M demo

% Knobs
fps = 8; % frames per second
res = 500; % # pixels width
title = 'Synthetic Block M [Ground Truth]';

% Show video
load('./Block M sims/BlockM_GT.mat');
M = Mtrue;
DataViewer(M,fps,res,title);

%% 4b. Simuated Block-M demo

% Knobs
fps = 8; % frames per second
res = 500; % # pixels width
title = {'Synthetic Block M [Proposed]';
         'Synthetic Block M [IST]'};

% Show video
load('./Block M sims/BlockM_SLRN.mat'); % SLRN r = 1
%load('./Block M sims/BlockM_SLRN2.mat'); % SLRN r = 2
M1 = abs(Ltrue + Strue);
load('./Block M sims/BlockM_ist.mat');
M2 = abs(Ltrue + Strue);
DataViewer({M1;M2},fps,res,title);
