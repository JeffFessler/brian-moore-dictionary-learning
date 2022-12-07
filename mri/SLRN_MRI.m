function [L S cost deltaM time its NRMSE SNR] = SLRN_MRI(E,T,d,r,lambdaS,tol,MaxIters,Mtrue)
% 
% [L S cost deltaM time its] = SLRN_MRI(E,T,d,r,lambdaS,tol,MaxIters);
% [L S cost deltaM time its NRMSE SNR] = SLRN_MRI(E,T,d,r,lambdaS,tol,MaxIters,Mtrue);
% 
% d: undersampled k-t data (nx,ny,nt,nc)
% E: data acquisition operator
% T: sparsifying transform
% r: desired rank of low-rank component (r = 0 employs RankDetect)
% lambdaS: l1-norm weight
% MaxIters: max number of allowed iterations
% tol: stopping tolerance
% [OPTIONAL] Mtrue: true target image
% 
% New sparse + low-rank decomposition algorithm with optimal low-rank
% shrinkage
% 
% Brian Moore (2013)
% 

% Check to see if we should compute NRMSEs and SNRs
if (nargin == 8)
    HaveTrueImages = true;
    Mtrue = abs(Mtrue); % absolute value
    %Mtrue = Mtrue - min(Mtrue(:)); % min value = 0
    %Mtrue = Mtrue / max(Mtrue(:)); % max value = 1
    NRMSEfcn = @(M) 100 * norm(M(:) - Mtrue(:),2) / norm(Mtrue(:),2);
    SNRfcn = @(M) 20 * log10(norm(Mtrue(:),2) / norm(M(:) - Mtrue(:),2));
else
    HaveTrueImages = false;
end

%--------------------------------------------------------------------------
% Parse inputs
%--------------------------------------------------------------------------
[nx,ny,nt,~] = size(d);
nd = nx * ny;

M = reshape(E' * d,[nd nt]);
Lpre = M;
L = zeros(nd,nt);
S = zeros(nd,nt);
%--------------------------------------------------------------------------

% L + S iterations
its = 0;
cost = [];
deltaM = [];
time = [];
NRMSE = [];
SNR = [];
fprintf('\n*** SLRN MRI ***\n')
while(1)
    % Record iteration time
    tic;
    
    % Update iteration count
	its = its + 1;
    
    % Save last iterate
	Mlast = M;
    
    %----------------------------------------------------------------------
	% Low-rank update
    %----------------------------------------------------------------------
    % Empirical OptShrink shrinkage
	[UL SL VL] = svd(M - S,0);
    sL = diag(SL);
    if (r == 0)
        % Use RankDetect to determine rank
        tau = @(q) min(1,2 / q^(1/4));
        [~,Signals] = RankDetectSVD(sL,nd,nt,tau);
    else
        % Use fixed rank from user
        Signals = 1:r;
    end
    sL = ComputeWeightsSVD(sL,nd,nt,Signals);
    %sL = sL(1:r); % PCA
	L = UL(:,Signals) * diag(sL) * VL(:,Signals)';
    %----------------------------------------------------------------------
    
	% Sparse update
	%S = reshape(T' * soft(T * reshape(M - L,[nx ny nt]),lambdaS),[nd nt]);
	S = reshape(T' * soft(T * reshape(M - Lpre,[nx ny nt]),lambdaS),[nd nt]);
    
	% Data consistency update
	resid = E * reshape(L + S,[nx ny nt]) - d;
	M = L + S - reshape(E' * resid,[nd nt]);
    
    % L_{k-1} for the next iteration
    Lpre = L;
    
    % Print algorithm status
    cost(its,1) = nan; %#ok  cost function
    deltaM(its,1) = norm(M(:) - Mlast(:),2) / norm(Mlast(:),2); %#ok relative change
    time(its,1) = toc; %#ok
    fprintf('Iter: %d, deltaM: %f3\n',its,deltaM(its,1)); 
    
    % Compute NRMSE/SNRs
    if HaveTrueImages
        Mhat = abs(L + S);
        %Mhat = Mhat - min(Mhat(:)); % min value = 0
        %Mhat = Mhat / max(Mhat(:)); % max value = 1
        NRMSE(its,1) = NRMSEfcn(Mhat); %#ok
        SNR(its,1) = SNRfcn(Mhat); %#ok
    end
    
    % Stopping criteria 
    if ((its >= MaxIters) || (deltaM(its,1) < tol))
        break;
    end
end

% Return L and S as images over time
L = reshape(L,[nx ny nt]);
S = reshape(S,[nx ny nt]);

end

% Soft-thresholding
function sY = soft(Y,tau)

sY = max(abs(Y) - tau,0) .* sign(Y);
sY(isnan(sY)) = 0;

end
