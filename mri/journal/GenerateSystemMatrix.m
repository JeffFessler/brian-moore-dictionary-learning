function A = GenerateSystemMatrix(ny,nx,nt)
%--------------------------------------------------------------------------
% Syntax:       A = GenerateSystemMatrix(ny,nx,nt);
%               
% Inputs:       (ny,nx,nt,nc) = (#rows,#cols,#frames)
%               
% Outputs:      A is an multicoil encoding operator 
%               
% Description:  This function generates a multicoil encoding operator of
%               the given input dimensions that mimics the sampling
%               pattern and coil senstitivies from the Otazo L + S dataset
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         August 13, 2014
%--------------------------------------------------------------------------

% Load Otazo data
data = load('../cardiac_perf_R8.mat');

% Compute empirical sampling CDF
Fx = mean(cumsum(squeeze(data.kdata(1,:,:,1)) ~= 0,1),2);
n = length(Fx);
F = CDFgen({1:nx,CDFeval(CDFgen({1:n,Fx}),linspace(1,n,nx))});

% Generate new masks
M = false(ny,nx,nt);
randf = @(n) min(max(round(CDFsample(F,n)),1),nx);
for i = 1:nt
    while (nnz(M(1,:,i)) < 16)
        M(:,randf(1),i) = true;
    end
end

%--------------------------------------------------------------------------
% Generate sensitivity maps
%--------------------------------------------------------------------------
% Method #1: based on Otazo's maps
nc = size(data.b1,3);
S = zeros(ny,nx,nc);
for k = 1:nc
    S(:,:,k) = imresize(data.b1(:,:,k),[ny nx]);
end

% Method #2: 1 coil w/ rerfect sensitivity
%S = ones(ny,nx);
%--------------------------------------------------------------------------

% Generate system matrix
A = Emat_xyt(M,S);
