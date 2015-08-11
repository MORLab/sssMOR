function varargout = entermatrixname(varargin)
% GUI zur Eingabe eines Matrixnames für load ascii data

% letzte Änderung 08.10.2010

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @entermatrixname_OpeningFcn, ...
                   'gui_OutputFcn',  @entermatrixname_OutputFcn, ...
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


% --- Executes just before entermatrixname is made visible.
function entermatrixname_OpeningFcn(hObject, temp1, handles, varargin)
handles.output = '?';
set(handles.ed_matname,'String',varargin{1})
set(handles.ed_matname,'UserData',0)
% Update handles structure
guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = entermatrixname_OutputFcn(hObject, eventdata, handles) 
if ~isempty(handles)
    varargout{1} = handles.output;
    delete(hObject);
else
    varargout{1} =[];
end

function ed_matname_Callback(hObject, eventdata, handles)
matname=get(hObject,'String');
if isvarname(matname)
    set(hObject,'UserData',0)
else
    errordlg('You have to enter a valid variable name','Error Dialog','modal')
    set(hObject,'UserData',1)
end

function ed_matname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pb_ok_Callback(hObject, eventdata, handles)
if get(handles.ed_matname,'UserData')==0
    k=[];
    p=get(handles.ed_matname,'String');
    %überprüfen ob name bereits im workspace vorhanden
    if exist_in_base_ws(p)
      s=sprintf('%s already exists in base workspace. Do you want to overwrite it?',p);
      k=stqd('String',s,'Title','Question Dialog'); %Aufruf von Standard Question Dialogbox
    end
    if ~isempty(k) && strcmp(k,'No')
        return % Workspace Variable soll nicht überschrieben werden
    end
    handles.output=p; % Variablennamen an die Aufrufende Funktion übergeben
    guidata(hObject, handles); % guidata aktualisieren (damit output gespeichert wird)
    uiresume(handles.figure1) % figure schließen
else
    errordlg('You have to enter a valid variable name','Error Dialog','modal')
end


