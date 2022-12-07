function GenerateMovie1(X,outpath,fps,watermark,res,lims)
%--------------------------------------------------------------------------
% Syntax:       GenerateMovie1(X,outpath,fps,watermark);
%               GenerateMovie1(X,outpath,fps,watermark,res);
%               GenerateMovie1(X,outpath,fps,watermark,res,lims);
%               
% Inputs:       X is an (ny x nx x nt) cube containing nt video frames
%               
%               outpath is the desired output filename/path (without
%               extension)
%               
%               fps is the desired frames per second for the movie
%               
%               watermark is a string (or cell array of strings) to display
%               at the bottom of each frame
%               
%               [OPTIONAL] res is the desired vertical resolution, in
%               pixels. The default value (or when res = []) is ny
%               
%               [OPTIONAL] lims = [minlim maxlim] contains the desired
%               limits for image display. The default value (or when lims 
%               = []) is lims = [min(X(:)) max(X(:))]
%               
% Description:  This function generates an .avi video of the input data X
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         August 25, 2014
%--------------------------------------------------------------------------

% Watermark knobs
wpad = 2;
color = 'w';
fontSize = 8;
fontName = 'Helvetica';
fontWeight = 'normal';

% Video formats
%format = 'Motion JPEG AVI';
format = 'Uncompressed AVI';
%format = 'Archival';
%format = 'Motion JPEG 2000';
%format = 'MPEG-4';

% Parse inputs
if ischar(watermark)
    watermark = {watermark};
end
if ~exist('res','var') || isempty(res)
    % Default resolution
    res = size(X,1);
end
[ny nx] = size(imresize(X(:,:,1),[res nan],'nearest'));
nt = size(X,3);
if ~exist('lims','var') || isempty(lims)
    % Default display limits
    lims = [min(X(:)) max(X(:))];
end

% Create video object (could also use avifile() function)
vidObj = VideoWriter(outpath,format);
vidObj.FrameRate = fps;

% Open the video file
vidObj.open();

% Set up a nice figure for recording
scrsz = get(0,'ScreenSize');
xyc = 0.5 * scrsz(3:4);
dim = [nx ny];
f = figure('Name','', ...
           'MenuBar','None', ...
           'NumberTitle','off', ...
           'Position',[(xyc - 0.5 * dim) dim], ...
           'Renderer','zbuffer', ...
           'Resize','off');
ax = axes('Position',[0 0 1 1]);
axis(ax,'off');

% Record the movie
for k = 1:nt
    % Display the kth frame
    hold off;
    Xk = imresize(X(:,:,k),[res nan],'nearest');
    imshow(Xk,lims,'Border','tight','Parent',ax);
    hold on;
    
    % Display watermark, if necessary
    if ~isempty(watermark)
        text(nx - wpad - 1,ny - wpad,watermark, ...
            'Color',color, ...
            'FontSize',fontSize, ...
            'FontWeight',fontWeight, ...
            'FontName',fontName, ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','bottom');
    end
    
    % Write frame to file
    drawnow; pause(0.01);
    frame = getframe(f);
    vidObj.writeVideo(frame);
end
close(f);

% Close the .avi file
vidObj.close();
