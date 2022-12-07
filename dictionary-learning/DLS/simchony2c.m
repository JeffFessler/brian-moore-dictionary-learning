function F = simchony2c(P,Q)
%
% Syntax:       F = simchony2c(P,Q);
%               
% Inputs:       P and Q are m x n matrices containing x and y gradients,
%               respectively
%               
% Outputs:      F is an m x n matrix containing the estimated surface
%               
% Description:  Computes the solution to the least squares problem
%               
%               \min_f 0.5\|b - Af\|_2^2
%               
%               where A is the equivalent matrix from
%               http://ubee.enseeiht.fr/photometricstereo/Matlab/simchony.m
%               with Neumann boundary conditions
%               
% Reference:    Simchony, T., Chellappa, R., and Shao, M., "Direct
%               analytical methods for solving Poisson equations in
%               computer vision problems", IEEE Transactions on Pattern
%               Analysis and Machine Intelligence, 12(5), 435-446, 1990
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         May 12, 2016
%

% Parse inputs
[m, n] = size(P);

% Integrability + Neumann boundary conditions
B = integrabilityNeumann(P,Q);

% Solve least squares problem
[X, Y] = meshgrid(0:(n - 1),0:(m - 1));
Z      = dct2(B) ./ (2 * (cos(pi * X / n) + cos(pi * Y / m)) - 4);
Z(1,1) = 0; % Zero offset
F      = idct2(Z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Integrability model with Neumann boundary conditions
function B = integrabilityNeumann(P,Q)

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
