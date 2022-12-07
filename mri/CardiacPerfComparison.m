%
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

tol = 0.0025;
MaxIters = 50;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Load data
%--------------------------------------------------------------------------
load('cardiac_perf_R8.mat');
[nx ny nt nc] = size(kdata);
nd = nx * ny;
E = Emat_xyt(kdata(:,:,:,1) ~= 0,b1); % Aquisition operator
T = TempFFT(3); % Sparsifying transformation : 3D temporal DFT
d = kdata; % Sub-sampled k-space data
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

SLRN.r = nan(NLambdaL,NLambdaS);
SLRN.lambdaS = nan(NLambdaL,NLambdaS);
SLRN.L = cell(NLambdaL,NLambdaS);
SLRN.S = cell(NLambdaL,NLambdaS);
SLRN.cost = cell(NLambdaL,NLambdaS);
SLRN.deltaM = cell(NLambdaL,NLambdaS);
SLRN.time = cell(NLambdaL,NLambdaS);
SLRN.its = nan(NLambdaL,NLambdaS);
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
        [L S cost deltaM time its] = lps_ist(E,T,d,C * lambdaL,lambdaS,tol,MaxIters);
        IST.lambdaL(i,j) = lambdaL;
        IST.lambdaS(i,j) = lambdaS;
        IST.L{i,j} = L;
        IST.S{i,j} = S;
        IST.cost{i,j} = cost;
        IST.deltaM{i,j} = deltaM;
        IST.time{i,j} = time;
        IST.its(i,j) = its;
        
        % Our SLRN algo
        [L S cost deltaM time its] = SLRN_MRI(E,T,d,r,lambdaS,tol,MaxIters);
        SLRN.r(i,j) = r;
        SLRN.lambdaS(i,j) = lambdaS;
        SLRN.L{i,j} = L;
        SLRN.S{i,j} = S;
        SLRN.cost{i,j} = cost;
        SLRN.deltaM{i,j} = deltaM;
        SLRN.time{i,j} = time;
        SLRN.its(i,j) = its;
    end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Save data
%--------------------------------------------------------------------------
%save('C:\Users\brimoor\Desktop\CardiacPerfComparison.mat','IST','SLRN','-v7.3');
save('CardiacPerfComparison.mat','IST','SLRN','-v7.3');
%--------------------------------------------------------------------------
