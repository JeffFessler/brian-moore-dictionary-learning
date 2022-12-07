function [L S cost deltaM time its NRMSE SNR] = lps_ist(E,T,d,lambdaL,lambdaS,tol,MaxIters,Mtrue)
% 
% [L S cost deltaM time its] = lps_ist(E,T,d,lambdaL,lambdaS,tol,MaxIters);
% [L S cost deltaM time its NRMSE SNR] = lps_ist(E,T,d,lambdaL,lambdaS,tol,MaxIters,Mtrue);
% 
% d: undersampled k-t data (nx,ny,nt,nc)
% E: data acquisition operator
% T: sparsifying transform
% lambdaL: nuclear-norm weight
% lambdaS: l1-norm weight
% tol: stopping tolerance
% MaxIters: max number of allowed iterations
% [OPTIONAL] Mtrue: true target image
% 
% L + S reconstruction of undersampled dynamic MRI data using iterative
% soft-thresholding of singular values of L and entries of TS
% 
% Ricardo Otazo (2013)
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

% Initializations
M = E' * d;
[nx ny nt] = size(M);
M = reshape(M,[nx * ny,nt]);
Lpre = M;
L = zeros(nx * ny,nt); 
S = zeros(nx * ny,nt); 
its = 0;

% L + S iterations
cost = [];
deltaM = [];
time = [];
NRMSE = [];
SNR = [];
fprintf('\n*** L + S IST ***\n')
while(1)
    % Record iteration time
    tic;
    
    % Update iteration count
	its = its + 1;
    
    % Save last iterate
	Mlast = M;
    
	% Low-rank update
	[UL SL VL] = svd(M - S,0);
    sL = SoftThresh(diag(SL),lambdaL);
    %sL = SoftThresh(diag(SL),SL(1) * lambdaL); % This is not what your paper said you were doing
	L = UL * diag(sL) * VL';
    
	% Sparse update
	%S = reshape(T' * SoftThresh(T * reshape(M - L,[nx ny nt]),lambdaS),[nx * ny,nt]);
	S = reshape(T' * SoftThresh(T * reshape(M - Lpre,[nx ny nt]),lambdaS),[nx * ny,nt]);
    
	% Data consistency update
	resid = E * reshape(L + S,[nx ny nt]) - d;
	M = L + S - reshape(E' * resid,[nx * ny,nt]);
    
    % L_{k-1} for the next iteration
    Lpre = L;
    
    % Print algorithm status
	temp = T * reshape(S,[nx ny nt]);
    cost(its,1) = norm(resid(:),2)^2 + lambdaL * sum(sL) + lambdaS * norm(temp(:),1); %#ok  cost function
    deltaM(its,1) = norm(M(:) - Mlast(:),2) / norm(Mlast(:),2); %#ok relative change
    time(its,1) = toc; %#ok
    fprintf('Iter: %d, Cost: %f3, deltaM: %f3\n',its,cost(its,1),deltaM(its,1)); 
    
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

% Soft-thresholding function
function y = SoftThresh(x,p)

y = max(abs(x) - p,0) .* sign(x);
y(isnan(y)) = 0;

end
