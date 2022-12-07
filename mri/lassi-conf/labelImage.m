function labelImage(I,xlabels,ylabels,gap,fontSize)
% Syntax:   labelImage(I,xlabels,ylabels);
%           labelImage(I,xlabels,ylabels,gap,fontSize);
%           labelImage(I,xlabels,ylabels,gap);

% Parse inputs
if ~exist('gap','var') || isempty(gap)
    % Default gap
    gap = 2; 
end
if ~exist('fontSize','var') || isempty(fontSize)
    % Default font size
    fontSize = 12;
end

% Label function
txt = @(x,y,str,ha,va) text(x,y,str,'HorizontalAlignment',ha, ...
                                    'VerticalAlignment',va, ...
                                    'FontName','Helvetica', ...
                                    'FontSize',fontSize, ...
                                    'FontWeight','normal', ...
                                    'Rotation',0);

% Display image
figure;
imshow(I);
hold on;
lim = axis();

% x labels
Nx  = numel(xlabels);
dx  = (lim(2) - lim(1)) / Nx;
for i = 1:Nx
    txt((i - 0.5) * dx,-gap,xlabels{i},'center','bottom');
end

% y labels
Ny  = numel(ylabels);
dy  = (lim(4) - lim(3)) / Ny;
for j = 1:Ny
    txt(-gap,(j - 0.5) * dy,ylabels{j},'right','middle');
end
