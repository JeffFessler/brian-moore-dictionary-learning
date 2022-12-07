function normA = get_normA(nt,nf)
% Syntax:   normA = get_normA(nt,nf);

% Process based on total # frames
switch nf
    case 480
        % ir_mri_dce_obj1
        ntl = [2 4 6 8 10 12 16 20 24 30 32 40];
        normAl = [22.985 21.143 20.805 20.638 20.543 20.526 20.402 20.378 20.317 20.299 20.344 20.274];
    case 168
        % ir_mri_dce_obj2
        ntl = [2 4 6 8 12 14 21 24 28 42];
        normAl = [20.780 20.516 20.419 20.307 20.263 20.253 20.266 20.235 20.224 20.197];
    case 336
        % ir_mri_dce_obj3
        ntl = [2 4 6 8 12 16 21 24 28 42];
        normAl = [21.586 20.841 20.565 20.516 20.419 20.307 20.292 20.297 20.275 20.266];
    otherwise
    % Pre-compute error
    error('normA(%i,%i) wasn''t precomputed...',nt,nf);
end

% Return ||A||
normA = normAl(nt == ntl);
if isempty(normA)
    % Pre-compute error
    error('normA(%i,%i) wasn''t precomputed...',nt,nf);
end

%{
% Compute normA
AAtfcn = @(x) A' * (A * x);
opts = struct('issym',true,'isreal',false,'maxit',5,'disp',1);
normA = sqrt(eigs(AAtfcn,size(A,2),1,'LM',opts)) %#ok
%}

%{
% Feasible scan times: 5 * t must be integer-valued
% Valid #frames
n = 168;
tmp = n ./ (1:n);
valid_nt = find(tmp == round(tmp))
%}
