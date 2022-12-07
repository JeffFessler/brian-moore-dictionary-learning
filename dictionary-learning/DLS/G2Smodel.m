function [A, b, normA, f] = G2Smodel(N,M,type)
%
% Syntax:       [A, b, normA, f] = G2Smodel(N,M);
%               [A, b, normA, f] = G2Smodel(N,M,type);
%               
% Inputs:       N is either a 3 x mn matrix or an m x n x 3 matrix
%               containing the surface normal vectors
%               
%               M is an m x n logical mask defining the support of the 
%               surface being reconstructed
%               
%               [OPTIONAL] type = {'1c','1n','2a','2b','2c'} determines
%               which surface model to apply. The default value is
%               type = '2a', and the options are
%               
%                   type = {'1c',       First order differences models with
%                           '1n'}       without circulant boundaries
%                   
%                   type = {'2a',       Second-order integrability models
%                           '2b',       with Neumann (natural) boundary
%                           '2c'}       conditions
%               
% Outputs:      (A,b,f) define a least-squares model of the form
%               
%               f = \argmin_f 0.5 \|b - Af\|_2^2
%               
%               that produces an estimate f = F(:) of the surface
%               corresponding to the input normal vectors, where
%               size(F) = [m, n]
%               
%               normA is (an upper bound on) the largest singular value
%               of A
%               
% Description:  Generates a least-squares model for the gradients-to-
%               surface (G2S) photometric stereo problem
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         May 12, 2016
%

% Parse inputs
[m, n] = size(M);
if ~exist('type','var') || isempty(type)
    % Default type
    type = '2a';
end

% Reshape normals, if necessary
if ndims(N) ~= 3
    % N is currently 3 x mn
    N = reshape(N',[m, n, 3]);
end

% Compute gradients
P = -N(:,:,1) ./ N(:,:,3);                  % x-gradients      
Q =  N(:,:,2) ./ N(:,:,3);                  % y-gradients
P(isnan(P) | ~M) = 0;
Q(isnan(Q) | ~M) = 0;

% Generate model
switch lower(type)
    case '1c'
        %
        % First order model
        %
        % Circulant boundary conditions
        %
        A     = [kron(Dc(n), speye(m));     % First differences across rows
                 kron(speye(n), Dc(m))];    % First differences down cols
        b     = [P(:); Q(:)];
        normA = sqrt(8);
    case '1n'
        %
        % First order model
        %
        % Non-circulant boundary conditions
        %
        A     = [kron(Dnc(n), speye(m));    % First differences across rows
                 kron(speye(n), Dnc(m))];   % First differences down cols
        b     = [vec(P(:,1:(n - 1))); vec(Q(1:(m - 1),:))];
        normA = sqrt(8);
    case '2a'
        %
        % Second order model
        % 
        % A matrix from pg 438 of simchony1990 with Neumann boundary
        % conditions
        %
        B     = Ba(m);
        d1    = [nan; ones(n - 1,1)];
        A     = kron(speye(n),B) + ...
                kron(spdiags([flipud(d1), d1],[-1, 1],n,n),speye(m));
        b     = integrabilityNeumann(P,Q);
        normA = 8;
    case '2b'
        %
        % Second order model
        % 
        % A matrix from pg 439 of simchony1990 with Neumann boundary
        % conditions
        %
        B     = Bb(m);
        d1    = [nan; 2; ones(n - 2,1)];
        A     = kron(speye(n),B) + ...
                kron(spdiags([flipud(d1), d1],[-1, 1],n,n),speye(m));
        b     = integrabilityNeumann(P,Q);
        normA = 8.21;
    case '2c'
        %
        % Second order model
        % 
        % Equivalent A matrix from
        % from http://ubee.enseeiht.fr/photometricstereo/Matlab/simchony.m
        % with Neumann boundary conditions
        %
        B1    = B1c(m);
        B2    = B2c(m);
        A     = blkdiag(B1,kron(speye(n - 2),B2),B1) + ...
                kron(spdiags(ones(n,2),[-1, 1],n,n),speye(m));
        b     = integrabilityNeumann(P,Q);
        normA = 8;
    otherwise
        % Unsupported type
        error('Type "%s" not supported',type);
end

% Solve model, if requested
if nargout >= 4
    switch lower(type)
    case '1c'
        % Must use lsqr, since \|b - Af\|^2 > 0
        f = lsqr(A,b,[],2000);
    case '1n'
        % Must use lsqr, since \|b - Af\|^2 > 0
        f = lsqr(A,b,[],2000);
    case '2a'
        % Fast solver via 2D discrete sine transform
        f = vec(simchony2a(P,Q));
    case '2b'
        % lsqr doesn't converge... write fast solver for this?
        f = nan(m * n,1);
    case '2c'
        % Fast solver via 2D discrete cosine transform
        f = vec(simchony2c(P,Q));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectorize data
function x = vec(X)
x = X(:);

% Circulant first-differences matrix
function D = Dc(n)
D      = spdiags([-ones(n,1), ones(n,1)],[0, 1],n,n);
D(n,1) = 1;

% Non-circulant first-differences matrix
function D = Dnc(n)
D = spdiags([-ones(n - 1,1), ones(n - 1,1)],[0, 1],n - 1,n);

% B matrix for method 2a
function B = Ba(m)
d0 = repmat(-4,m,1);                        % Main diagonal
d1 = [nan; ones(m - 1,1)];                  % Off diagonals
B  = spdiags([flipud(d1), d0, d1],[-1, 0, 1],m,m);

% B matrix for method 2b
function B = Bb(m)
d0 = repmat(-4,m,1);                        % Main diagonal
d1 = [nan; 2; ones(m - 2,1)];               % Off diagonals
B  = spdiags([flipud(d1), d0, d1],[-1, 0, 1],m,m);

% B1 matrix for model 2c
function B1 = B1c(m)
d0 = [-2; repmat(-3,m - 2,1); -2];          % Main diagonal
d1 = [nan; ones(m - 1,1)];                  % Off diagonals
B1 = spdiags([flipud(d1), d0, d1],[-1, 0, 1],m,m);

% B2 matrix for model 2c
function B2 = B2c(m)
d0 = [-3; repmat(-4,m - 2,1); -3];          % Main diagonal
d1 = [nan; ones(m - 1,1)];                  % Off diagonals
B2 = spdiags([flipud(d1), d0, d1],[-1, 0, 1],m,m);

% Integrability model with Neumann boundary conditions
function b = integrabilityNeumann(P,Q)

% Parse inputs
[m, n] = size(Q);

% Integrability model
dP = 0.5 * (P(:,[2:n, n]) - P(:,[1, 1:(n - 1)]));
dQ = 0.5 * (Q([2:m, m],:) - Q([1, 1:(m - 1)],:));
B  = dP + dQ;

% Neumann boundary condition
B(2:(m - 1),1) = 0.5 * ( P(2:(m - 1),1) + P(2:(m - 1),2)    );
B(2:(m - 1),n) = 0.5 * (-P(2:(m - 1),n) - P(2:(m - 1),n - 1));
B(1,2:(n - 1)) = 0.5 * ( Q(1,2:(n - 1)) + Q(2,2:(n - 1))    );
B(m,2:(n - 1)) = 0.5 * (-Q(m,2:(n - 1)) - Q(m - 1,2:(n - 1)));
B(1,1) = 0.5 * ( P(1,1) + P(1,2)     + Q(1,1) + Q(2,1)    );
B(m,1) = 0.5 * ( P(m,1) + P(m,2)     - Q(m,1) - Q(m - 1,1));
B(1,n) = 0.5 * (-P(1,n) - P(1,n - 1) + Q(1,n) + Q(2,n)    );
B(m,n) = 0.5 * (-P(m,n) - P(m,n - 1) - Q(m,n) - Q(m - 1,n));

% Vectorize
b = B(:);
