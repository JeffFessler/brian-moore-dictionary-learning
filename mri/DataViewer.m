function DataViewer(M,fps,res,titlestr)
% Syntax:   DataViewer(M);
%           DataViewer(M,fps);
%           DataViewer(M,fps,res);
%           DataViewer(M,fps,res,titlestr);

% Knobs
FontSize = 12;

% Parse inputs
if (nargin < 2)
    fps = 10;
end
if (nargin < 3)
    res = 500;
end
if (nargin < 4)
    titlestr = 'Frame';
end

% Show video
figure;
if ~iscell(M)
    % One data set
    minM = min(M(:));
    maxM = max(M(:));
    nt = size(M,3);
    for i = 1:nt
        imshow(imresize(M(:,:,i),[res NaN],'nearest'),[minM maxM]);
        title(sprintf('%s %i/%i',titlestr,i,nt),'FontSize',FontSize);
        pause(1 / fps);
    end
    close(gcf);
else
    % Compare two data sets
    M1 = M{1};
    M2 = M{2};
    minM1 = min(M1(:));
    maxM1 = max(M1(:));
    minM2 = min(M2(:));
    maxM2 = max(M2(:));
    nt1 = size(M1,3);
    nt2 = size(M2,3);
    for i = 1:max(nt1,nt2)
        if (i <= nt1)
            subplot(1,2,1);
            imshow(imresize(M1(:,:,i),[res NaN],'nearest'),[minM1 maxM1]);
            title(sprintf('%s %i/%i',titlestr{1},i,nt1),'FontSize',FontSize);
        end
        
        if (i <= nt2)
            subplot(1,2,2);
            imshow(imresize(M2(:,:,i),[res NaN],'nearest'),[minM2 maxM2]);
            title(sprintf('%s %i/%i',titlestr{2},i,nt2),'FontSize',FontSize);
        end
        
        pause(1 / fps);
    end
    close(gcf);
end
