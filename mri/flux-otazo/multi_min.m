function [Y I] = multi_min(X,dims)
% Syntax:   [Y I] = multi_min(X,dims);

% Rearrange so that averaging dimensions are last
sz = [size(X) ones(1,max(dims) - ndims(X))];
Nd = length(sz);
indims = setdiff(1:Nd,dims);
Y = permute(X,[indims dims]);

% Iteratively minimize over dimensions
Ndims = length(dims);
inds = cell(Ndims,1);
for i = Ndims:-1:1
    [Y inds{i}] = min(Y,[],ndims(Y));
end

% Construct minimizing index list
sz = size(Y);
Nd = length(sz);
if ((Nd == 2) && (size(Y,2) == 1))
    Nd = 1;
end
NY = numel(Y);
I = cell(sz);
for idx = 1:NY
    Yinds = cell(1,Nd);
    [Yinds{:}] = ind2sub(sz,idx);
    for i = 1:Ndims
        Yinds{Nd + i} = inds{i}(Yinds{:});
    end
    I{Yinds{1:Nd}} = [Yinds{(Nd + 1):end}];
end
