function varargout = truncateModel(varargin)
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                    -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Input Arguments:
%       *Required Input Arguments:*
%       -varargin{1}:       cell-array of size 1x1
%           {1,1}:          model that should be truncated
%             
% ------------------------------------------------------------------
% Authors:      Niklas Kochdumper
% Last Change:  23 June 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @truncateModel_OpeningFcn, ...
                   'gui_OutputFcn',  @truncateModel_OutputFcn, ...
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


%--------------------------------------------------------------------------
%                      OPENING AND CLOSING FUNCTION
%--------------------------------------------------------------------------

function truncateModel_OpeningFcn(hObject, eventdata, handles, varargin)

    % Choose default command line output for truncateModel
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    if nargin ~= 4
       errordlg('Wrong number of input arguments','Error Dialog','modal'); 
    end

    system = varargin{1}{1,1};

    if ~isa(system,'ss') && ~isa(system,'sss')
       errordlg('First input argument is an object of wrong type. Only ss-objects and sss-objects are supported','Error Dialog','modal');
    end
    
    %Set the selectable options for the inputs
    
    inputs = {};
    
    for i = 1:size(system.B,2)
       inputs{end+1,1} = i; 
    end

    set(handles.lb_input,'String',inputs);
    
    %Set the selectable options for the outputs
    
    outputs = {};
    
    for i = 1:size(system.C,1)
       outputs{end+1,1} = i; 
    end

    set(handles.lb_output,'String',outputs);
    
    %Set the default name for the truncated model
    
    setDefaultName(handles,system.Name);
    
    %Define a new field in the handles structure where the truncated model
    %is stored 
    
    handles.truncatedModel = [];
    handles.system = system;
    guidata(hObject,handles);
    
function varargout = truncateModel_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;

    
%--------------------------------------------------------------------------
%                               CALLBACKS
%--------------------------------------------------------------------------

function et_Callback(hObject, eventdata, handles)
    % check variable name
    isvalidvarname(hObject)

function pb_Callback(hObject, eventdata, handles)
% Truncate the model and save the truncated model to the workspace

    %Check if a valid name for the truncated system was selected by the
    %user
    
    if get(handles.et,'UserData')==1
        errordlg('Please correct name for reduced system first','Error Dialog','modal')
        uiwait
        return
    end
    
    %Truncate the model
    
    indicesInput = get(handles.lb_input,'Value');
    indicesOutput = get(handles.lb_output,'Value');
    
    sysTruncated = handles.system(indicesOutput,indicesInput);
    
    %Save the truncated model to the workspace
    
    name = get(handles.et,'String');
    assignin('base',name,sysTruncated);
    msgbox('The truncation of the model was successfull!','Information','modal');
    
    %Set a new default name for the system
    
    setDefaultName(handles,handles.system.Name);
    
    
    
%--------------------------------------------------------------------------
%                          AUXILIARY FUNCTIONS
%--------------------------------------------------------------------------
    
function setDefaultName(handles,systemName)
%This functino sets the default name for the truncated Model

    defaultName = strcat('sys_',systemName,'_truncated');
    newName = defaultName;
    counter = 1;
    
    while existInBaseWs(newName)
        
        newName = strcat(defaultName,num2str(counter));
        counter = counter + 1;      
    end
        
    defaultName = newName;
    set(handles.et,'String',defaultName);
        
        