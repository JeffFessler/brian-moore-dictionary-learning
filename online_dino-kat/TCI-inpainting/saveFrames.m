function saveFrames(X,basename,dim)
% Syntax:   saveFrames(X,basename);
% Syntax:   saveFrames(X,basename,dim);

% Parse basename
[path, name, ext] = fileparts(basename);
if ~isempty(path) && ~exist(path,'dir')
    mkdir(path);
end
RESIZE = exist('dim','var') && ~isempty(dim);

% Save frames
nt = size(X,3);
for k = 1:nt
    namek = sprintf('%s/%s-%d%s',path,name,k,ext);
    Xk = X(:,:,k);
    if RESIZE
        Xk = imresize(Xk,dim);
    end
    imwrite(Xk,namek,ext(2:end));
end
