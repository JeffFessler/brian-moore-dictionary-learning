function [L S cost deltaL time its NRMSE SNR] = SLRN_MRI_FS(E,T,d,r,lambdaS,tol,MaxIters,Mtrue,M0)
% 
% [L S cost deltaL time its] = SLRN_MRI_FS(E,T,d,r,lambdaS,tol,MaxIters,M0);
% [L S cost deltaL time its NRMSE SNR] = SLRN_MRI_FS(E,T,d,r,lambdaS,tol,MaxIters,Mtrue,M0);
% 
% d: FULLY sampled k-t data (nx,ny,nt,nc)
% E: data acquisition operator
% T: sparsifying transform
% r: desired rank of low-rank component (r = 0 employs RankDetect)
% lambdaS: l1-norm weight
% MaxIters: max number of allowed iterations
% tol: stopping tolerance
% [OPTIONAL] Mtrue: true target image
% [OPTIONAL] M0: E' * d
% 
% New sparse + low-rank decomposition algorithm with optimal low-rank
% shrinkage
% 
% Brian Moore (2013)
% 

% Check to see if we should compute NRMSEs and SNRs
if (nargin >= 8)
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

if (nargin == 9)
    M = reshape(M0,[nd nt]);
else
    M = reshape(E' * d,[nd nt]);
end
L = M;
S = zeros(nd,nt);
%--------------------------------------------------------------------------

% L + S iterations
its = 0;
cost = [];
deltaL = [];
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
	Llast = L;
    
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
	L = UL(:,Signals) * diag(sL) * VL(:,Signals)';
    %----------------------------------------------------------------------
    
	% Sparse update
	S = reshape(T' * soft(T * reshape(M - L,[nx ny nt]),lambdaS),[nd nt]);
	%S = reshape(T' * soft(T * reshape(M - Llast,[nx ny nt]),lambdaS),[nd nt]);
    
    % Print algorithm status
    cost(its,1) = nan; %#ok  cost function
    deltaL(its,1) = norm(L(:) - Llast(:),2) / norm(Llast(:),2); %#ok relative change
    time(its,1) = toc; %#ok
    fprintf('Iter: %d, deltaL: %f3\n',its,deltaL(its,1)); 
    
    % Compute NRMSE/SNRs
    if HaveTrueImages
        Mhat = abs(L + S);
        %Mhat = Mhat - min(Mhat(:)); % min value = 0
        %Mhat = Mhat / max(Mhat(:)); % max value = 1
        NRMSE(its,1) = NRMSEfcn(Mhat); %#ok
        SNR(its,1) = SNRfcn(Mhat); %#ok
    end
    
    % Stopping criteria 
    if ((its >= MaxIters) || (deltaL(its,1) < tol))
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
