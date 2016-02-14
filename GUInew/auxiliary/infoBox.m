function varargout = infoBox(varargin)
% GUI to display information to inform the user with a picture
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Input Argument:   Cell array with one cell containing the path to 
%                   the picture that should be displayed on the 
%                   Info Box
% ------------------------------------------------------------------
% Authors:      Niklas Kochdumper
% Last Change:  13 Feb 2015
% ------------------------------------------------------------------

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @infoBox_OpeningFcn, ...
                   'gui_OutputFcn',  @infoBox_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


function infoBox_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for infoBox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

test = nargin;

if nargin ~= 4
   errordlg('Wrong number of input arguments','Error Dialog','modal'); 
end

path = cell2mat(varargin{1,1});

h = imread(path);

pos = get(hObject,'Position');
pos(1,3) = size(h,2);
pos(1,4) = size(h,1);
set(hObject,'Position',pos);

posAxes = pos;
posAxes(1,1) = 0;
posAxes(1,2) = 0;
set(handles.axes,'Position',posAxes);

posPb = get(handles.pb,'Position');
posPb(1,1) = round((pos(1,3) - posPb(1,3))/2);
posPb(1,2) = 25;
set(handles.pb,'Position',posPb);

h = image(h);

set(handles.axes,'XTick',[]);
set(handles.axes,'YTick',[]);


function varargout = infoBox_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function pb_Callback(hObject, eventdata, handles)
    delete(handles.infoBox)
