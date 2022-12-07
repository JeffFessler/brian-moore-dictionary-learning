function M = cell2mov(C,gap,val)
% Syntax:   M = cell2mov(C,gap);
%           M = cell2mov(C,gap,val);
%
% NOTE: cell2mov(C,0) == cell2mat(C)
%

% Parse C
[m n] = size(C);
type = class(C{1});
if strcmpi(type,'logical')
    % Convert booleans to doubles
    type = 'double';
end
[~,~,dim3] = size(C{1});
width = sum(cellfun(@(c)size(c,2),C(1,:))) + (n - 1) * gap;

% Parse val
if ~exist('val','var') || isempty(val)
    % Default = max value
    maxVals = cellfun(@(c)max(c(:)),C);
    val = max(maxVals(:));
end

% Create movie
HORZ = val * ones(gap,width,dim3,type);
M = row2mov(C(1,:));
for i = 2:m
    M = cat(1,M,HORZ,row2mov(C(i,:)));
end

% Concatenate row
function row = row2mov(R)
    row = R{1};
    VERT = val * ones(size(R{1},1),gap,dim3,type);
    for j = 2:n
        row = cat(2,row,VERT,R{j});
    end
end

end
