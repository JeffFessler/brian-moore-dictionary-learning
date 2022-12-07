%%
% Multicoil cardiac perfusion MRI reconstruction algorithms
%

clear;
rng(1);

%--------------------------------------------------------------------------
% Knobs
%--------------------------------------------------------------------------
MinLambdaL = 0.001;
MaxLambdaL = 0.1;
NLambdaL = 11;

MinLambdaS = 0.001;
MaxLambdaS = 0.1;
NLambdaS = 11;

sigma = 1; % normalized noise variance

FullySampled = false; % {true,false}

tol = 0.0025;
MaxIters = 50;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Load data
%--------------------------------------------------------------------------
load('phantom.mat');
if (FullySampled == true)
    Mask = true(size(Mask)); % fully sampled!
end
E = Emat_xyt(Mask,CoilSensitivities);
T = TempFFT(3);
ktrue = E * reshape(Mtrue,[nx ny nt]); % true subsampled k-space data

% Add synthetic noise
idx = (ktrue ~= 0);
noise = (sigma / sqrt(2 * nd)) * (randn(size(ktrue)) + 1i * randn(size(ktrue)));
d = zeros(size(ktrue));
d(idx) = ktrue(idx) + noise(idx);
if (FullySampled == true)
    M0 = E' * d;
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Simulation lists
%--------------------------------------------------------------------------
lambdaLlist = logspace(log10(MinLambdaL),log10(MaxLambdaL),NLambdaL);
lambdaSlist = logspace(log10(MinLambdaS),log10(MaxLambdaS),NLambdaS);
rlist = (0:(NLambdaL - 1));
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Data storage
%--------------------------------------------------------------------------
IST.lambdaL = nan(NLambdaL,NLambdaS);
IST.lambdaS = nan(NLambdaL,NLambdaS);
IST.L = cell(NLambdaL,NLambdaS);
IST.S = cell(NLambdaL,NLambdaS);
IST.cost = cell(NLambdaL,NLambdaS);
IST.deltaM = cell(NLambdaL,NLambdaS);
IST.time = cell(NLambdaL,NLambdaS);
IST.its = nan(NLambdaL,NLambdaS);
IST.NRMSE = cell(NLambdaL,NLambdaS);
IST.NRMSEend = nan(NLambdaL,NLambdaS);
IST.SNR = cell(NLambdaL,NLambdaS);
IST.SNRend = nan(NLambdaL,NLambdaS);

SLRN.r = nan(NLambdaL,NLambdaS);
SLRN.lambdaS = nan(NLambdaL,NLambdaS);
SLRN.L = cell(NLambdaL,NLambdaS);
SLRN.S = cell(NLambdaL,NLambdaS);
SLRN.cost = cell(NLambdaL,NLambdaS);
SLRN.deltaM = cell(NLambdaL,NLambdaS);
SLRN.time = cell(NLambdaL,NLambdaS);
SLRN.its = nan(NLambdaL,NLambdaS);
SLRN.NRMSE = cell(NLambdaL,NLambdaS);
SLRN.NRMSEend = nan(NLambdaL,NLambdaS);
SLRN.SNR = cell(NLambdaL,NLambdaS);
SLRN.SNRend = nan(NLambdaL,NLambdaS);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Perform simulations
%--------------------------------------------------------------------------
C = max(svd(reshape(E' * d,[nd nt]))); % scaling constant for Otazo's SVT
for i = 1:NLambdaL
    for j = 1:NLambdaS
        % Display progress
        fprintf('\nTrial %i/%i\n',(i - 1) * NLambdaS + j,NLambdaL * NLambdaS);
        
        % Set parameters
        r = rlist(i);
        lambdaL = lambdaLlist(i);
        lambdaS = lambdaSlist(j);
        
        % Otazo's L + S algo
        if (FullySampled == true)
            % Fully-sampled version
            % Note: deltaM actually contains deltaL values...
            [L S cost deltaM time its NRMSE SNR] = lps_ist_FS(E,T,d,C * lambdaL,lambdaS,tol,MaxIters,Mtrue,M0);
        else
            % Undersampled version
            [L S cost deltaM time its NRMSE SNR] = lps_ist(E,T,d,C * lambdaL,lambdaS,tol,MaxIters,Mtrue);
        end
        IST.lambdaL(i,j) = lambdaL;
        IST.lambdaS(i,j) = lambdaS;
        IST.L{i,j} = L;
        IST.S{i,j} = S;
        IST.cost{i,j} = cost;
        IST.deltaM{i,j} = deltaM;
        IST.time{i,j} = time;
        IST.its(i,j) = its;
        IST.NRMSE{i,j} = NRMSE;
        IST.NRMSEend(i,j) = NRMSE(end);
        IST.SNR{i,j} = SNR;
        IST.SNRend(i,j) = SNR(end);
        
        % Our SLRN algo
        if (FullySampled == true)
            % Fully-sampled version
            % Note: deltaM actually contains deltaL values...
            [L S cost deltaM time its NRMSE SNR] = SLRN_MRI_FS(E,T,d,r,lambdaS,tol,MaxIters,Mtrue,M0);
        else
            % Undersampled version
            [L S cost deltaM time its NRMSE SNR] = SLRN_MRI(E,T,d,r,lambdaS,tol,MaxIters,Mtrue);
        end
        SLRN.r(i,j) = r;
        SLRN.lambdaS(i,j) = lambdaS;
        SLRN.L{i,j} = L;
        SLRN.S{i,j} = S;
        SLRN.cost{i,j} = cost;
        SLRN.deltaM{i,j} = deltaM;
        SLRN.time{i,j} = time;
        SLRN.its(i,j) = its;
        SLRN.NRMSE{i,j} = NRMSE;
        SLRN.NRMSEend(i,j) = NRMSE(end);
        SLRN.SNR{i,j} = SNR;
        SLRN.SNRend(i,j) = SNR(end);
    end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Save data
%--------------------------------------------------------------------------
save(['C:\Users\brimoor\Desktop\PhantomComparison_sigma_' int2str(round(sigma)) '.mat'],'IST','SLRN','-v7.3');
%save(['PhantomComparison_sigma_' int2str(round(sigma)) '.mat'],'IST','SLRN','-v7.3');
%--------------------------------------------------------------------------

%% Extract best NRMSE/SNR

% Knobs
sigma = 1; % normalized noise variance

[nx ny] = size(IST.NRMSE);
NRMSE_IST = inf;
SNR_IST = -inf;
idx_IST = [nan nan];
NRMSE_SLRN = inf;
SNR_SLRN = -inf;
idx_SLRN = [nan nan];
for i = 1:nx
    for j = 1:ny
        % IST
        NRMSE_IST = min(NRMSE_IST,min(IST.NRMSE{i,j}));
        temp_IST = max(IST.SNR{i,j});
        if (temp_IST > SNR_IST)
            SNR_IST = temp_IST;
            idx_IST = [i j];
        end
        
        % SLRN
        NRMSE_SLRN = min(NRMSE_SLRN,min(SLRN.NRMSE{i,j}));
        temp_SLRN = max(SLRN.SNR{i,j});
        if (temp_SLRN > SNR_SLRN)
            SNR_SLRN = temp_SLRN;
            idx_SLRN = [i j];
        end
    end
end
NRMSE_IST2 = min(min(IST.NRMSEend));
SNR_IST2 = max(max(IST.SNRend));
NRMSE_SLRN2 = min(min(SLRN.NRMSEend));
SNR_SLRN2 = max(max(SLRN.SNRend));

%--------------------------------------------------------------------------
% Display results
%--------------------------------------------------------------------------
% NRMSE
NRMSE_IST = NRMSE_IST %#ok
NRMSE_SLRN = NRMSE_SLRN %#ok
NRMSE_IST2 = NRMSE_IST2 %#ok
NRMSE_SLRN2 = NRMSE_SLRN2 %#ok

% SNR
SNR_IST = SNR_IST %#ok
SNR_SLRN = SNR_SLRN %#ok
SNR_IST2 = SNR_IST2 %#ok
SNR_SLRN2 = SNR_SLRN2 %#ok

% Maximizing indices
idx_IST = idx_IST %#ok
idx_SLRN = idx_SLRN %#ok
%--------------------------------------------------------------------------

% Save NRMSEs and SNRs
save(['PhantomComparison_sigma_' int2str(round(sigma)) '_Results.mat'],'NRMSE_IST','NRMSE_IST2','NRMSE_SLRN','NRMSE_SLRN2','SNR_IST','SNR_IST2','SNR_SLRN','SNR_SLRN2','idx_IST','idx_SLRN');
