function GenerateMovie_labels(X,outpath,fps,lims)
%--------------------------------------------------------------------------
% Syntax:       GenerateMovie_labels(X,outpath,fps);
%               GenerateMovie_labels(X,outpath,fps,lims);
%               
% Inputs:       X is an (nx x ny x nt) cube containing nt video frames
%               
%               outpath is the desired output filename/path (without
%               extension)
%               
%               fps is the desired frames per second for the movie
%               
%               [OPTIONAL] lims = [minlim maxlim] contains the desired
%               limits for image display. The default value is
%               lims = [min(X(:)) max(X(:))]
%               
% Description:  This function generates an .avi video of the input data X
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         January 6, 2015
%--------------------------------------------------------------------------

% Knobs
fontSize = 14;
ds = 30;
%format = 'Motion JPEG AVI';
format = 'Uncompressed AVI';
%format = 'Archival';
%format = 'Motion JPEG 2000';
%format = 'MPEG-4';

% Parse inputs
if ((nargin < 4) || isempty(lims))
    lims = [min(X(:)) max(X(:))];
end
[ny nx nt] = size(X);
dim = 2 * ds + [nx ny];

% Create video object - NOTE: could also use avifile()
vidObj = VideoWriter([outpath '.avi'],format);
vidObj.FrameRate = fps;

% Open .avi file
vidObj.open();

% Set up nice figure
scrsz = get(0,'ScreenSize');
xyc = 0.5 * scrsz(3:4);
f = figure('Name','', ...
           'MenuBar','None', ...
           'NumberTitle','off', ...
           'Position',[(xyc - 0.5 * dim) dim], ...
           'Renderer','zbuffer', ...
           'Color','w', ...
           'Resize','off');
ax = gca;
set(ax,'Units','pixels');
set(ax,'Position',[ds ds nx ny]);
axis(ax,[0 nx 0 ny]);
axis(ax,'off');

% Text label function
txt = @(x,y,str,theta) text(x,y,str,'HorizontalAlignment','center', ...
                                    'VerticalAlignment','middle', ...
                                    'FontSize',fontSize, ...
                                    'FontWeight','normal', ...
                                    'rotation',theta);
                                    %'FontName','Helvetica', ...

% Record movie
for k = 1:nt
    % Display kth frame
    imshow(flipud(X(:,:,k)),lims,'Border','tight','Parent',ax);
    set(ax,'YDir','normal');
    
    % Add column labels
    txt(1 * nx / 6,ny + ds / 2,'L + S',0);
    txt(3 * nx / 6,ny + ds / 2,'L',0);
    txt(5 * nx / 6,ny + ds / 2,'S',0);
    
    % Add row labels
    txt(-ds / 2,1 * ny / 4,'OptShrink',90);
    txt(-ds / 2,3 * ny / 4,'SVT',90);
    
    % Write frame to file
    drawnow;
    frame = getframe(f);
    vidObj.writeVideo(frame);
end
close(f);

% Close .avi file
vidObj.close();
