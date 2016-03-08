function fh = nicefigure(varargin)

% Creates nice figure (or makes an existing figure niiiiice)
%
% fh = nicefigure(title)
% fh = nicefigure(fh, title)
% fh = nicefigure(fh)

%   parse the input
if nargin == 0 
    title = 'Figure';
elseif nargin == 1 
    if ishandle(varargin{1})
        %existing window was passed
        fh = varargin{1};
        title = 'Figure';
    else
        % desired figure name was passed
        title = varargin{1};
    end
else
    % 2 arguments passed
    fh = varargin{1};
    title = varargin{2};
end

%   Make the nicefigure
if ~exist('fh','var')
    fh = figure('Color','white','Name',title,'NumberTitle','off');
else
    set(fh,'Color','white','Name',title,'NumberTitle','off');
end

