%%
% Multicoil MRI applied to synthetic block M data
%
% NOTE: Run this script from the /MRI tests folder
%

%--------------------------------------------------------------------------
% Knobs
%--------------------------------------------------------------------------
% IST [Otazo]
lambdaL_ist = 0.01;
lambdaS_ist = 0.01;

% SLRN [Proposed]
r = 1; % 0 ==> RankDetect
lambdaS_slrn = 0.01;

% Normalized noise variance
sigma = 1;

% Algo params
tol = 0.0025;
MaxIters = 50;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Load data
%--------------------------------------------------------------------------
load('./Block M sims/BlockM.mat');
E = Emat_xyt(Mask,CoilSensitivities);
T = TempFFT(3);
ktrue = E * reshape(Mtrue,[nx ny nt]); % true subsampled k-space data

% Add synthetic noise
idx = (ktrue ~= 0);
noise = (sigma / sqrt(2 * nd)) * (randn(size(ktrue)) + 1i * randn(size(ktrue)));
d = zeros(size(ktrue));
d(idx) = ktrue(idx) + noise(idx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Perform simulations
%--------------------------------------------------------------------------
% IST [Otazo]
C = max(svd(reshape(E' * d,[nd nt]))); % scaling constant for Otazo's SVT
[L S cost deltaM time its NRMSE SNR] = lps_ist(E,T,d,C * lambdaL_ist,lambdaS_ist,tol,MaxIters,Mtrue);
IST.lambdaL = lambdaL_ist;
IST.lambdaS = lambdaS_ist;
IST.L = L;
IST.S = S;
IST.cost = cost;
IST.deltaM = deltaM;
IST.time = time;
IST.its = its;
IST.NRMSE = NRMSE;
IST.NRMSEend = NRMSE(end);
IST.SNR = SNR;
IST.SNRend = SNR(end);

% SLRN [Proposed]
[L S cost deltaM time its NRMSE SNR] = SLRN_MRI(E,T,d,r,lambdaS_slrn,tol,MaxIters,Mtrue);
SLRN.r = r;
SLRN.lambdaS = lambdaS_slrn;
SLRN.L = L;
SLRN.S = S;
SLRN.cost = cost;
SLRN.deltaM = deltaM;
SLRN.time = time;
SLRN.its = its;
SLRN.NRMSE = NRMSE;
SLRN.NRMSEend = NRMSE(end);
SLRN.SNR = SNR;
SLRN.SNRend = SNR(end);
%--------------------------------------------------------------------------

%% Save data

%{
% IST
Ltrue = IST.L;
Strue = IST.S;
lambdaLtrue = IST.lambdaL;
lambdaStrue = IST.lambdaS;
save('./Block M sims/BlockM_ist.mat','Ltrue','Strue','lambdaLtrue','lambdaStrue');
%}

% SLRN
Ltrue = SLRN.L;
Strue = SLRN.S;
rtrue = SLRN.r;
lambdaStrue = SLRN.lambdaS;
save('./Block M sims/BlockM_SLRN.mat','Ltrue','Strue','rtrue','lambdaStrue');
