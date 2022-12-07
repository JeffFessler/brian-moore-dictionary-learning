function [L,S,ite,time,cost,mse,delta,amse] = lps_ist2(param)
% Syntax:   [L,S,ite,time,cost,mse,delta,amse] = lps_ist2(param);

% Check if we're using OptShrink
useOptShrink = isfield(param,'r');

% Initialize L and S
L = param.L0;
S = param.S0;
[nx, ny, nt] = size(L);
L = reshape(L,[nx * ny,nt]);
S = reshape(S,[nx * ny,nt]);

% Initialize M
resk = param.E * reshape(L + S,[nx,ny,nt]) - param.d;
M    = L + S - reshape(param.E' * resk,[nx * ny,nt]);

% Iterations
if useOptShrink
    disp('***** L + S (OptShrink) *****')
else
    disp('***** L + S (SVT) *****')
end
ite   = 0;
time  = nan(param.nite,1);
cost  = nan(param.nite,1);
mse   = nan(param.nite,1);
delta = nan(param.nite,1);
amse  = nan(param.nite,1);
while true
    tic;
	ite = ite + 1;
    
    % Last iterates
    M0    = M;
    Ltild = M - S;
    Stild = M - L;
    
	% L update
    if useOptShrink
        % OptShrink
        L = OptShrink(Ltild,param.r);
    else
        % SVT
        [Ut, St, Vt] = svd(Ltild,0);
        St = diag(SoftThresh(diag(St),param.lambda_L));
        L  = Ut * St * Vt';
    end
    
	% S update
	S = reshape(param.T' * (SoftThresh(param.T * reshape(Stild,[nx,ny,nt]),param.lambda_S)),[nx*ny,nt]);
    
	% M update
	resk = param.E * reshape(L + S,[nx,ny,nt]) - param.d;
	M    = L + S - reshape(param.E' * resk,[nx * ny,nt]);
    
    % Save stats
    time(ite)     = toc;
    if ~useOptShrink
        cost(ite) = 0.5 * norm(resk(:),2)^2 + param.lambda_L * sum(diag(St)) + param.lambda_S * norm(vec(param.T * reshape(S,[nx,ny,nt])),1);
    end
    mse(ite)      = norm(param.Lfull(:) - (L(:) + S(:)));
    delta(ite)    = norm(M(:) - M0(:)) / norm(M0(:));
    amse(ite)     = norm(abs(param.Lfull(:)) - abs(L(:) + S(:)));
    
    % Print stats
	fprintf(' ite: %03d, time: %.2f, cost: %.3f, mse: %.3f, amse: %.3f, delta: %.2e\n',ite,time(ite),cost(ite),mse(ite),amse(ite),delta(ite));
    
    % Check stopping criteria
    if (ite >= param.nite) || (delta(ite) < param.tol)
        break;
    end
end
L = reshape(L,nx,ny,nt);
S = reshape(S,nx,ny,nt);

time  = time(1:ite);
cost  = cost(1:ite);
mse   = mse(1:ite);
delta = delta(1:ite);
amse  = amse(1:ite);

% Soft thresholding
function y = SoftThresh(x,lambda)
y = max(abs(x) - lambda,0) .* sign(x);
y(isnan(y)) = 0;

% Vectorize data
function vx = vec(x)
vx = x(:);
