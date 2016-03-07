function varargout = mimoKrylov(varargin)
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                    -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Input Arguments:
%       *Required Input Arguments:*
%       -varargin{1}:       cell-array of size 1x3
%           {1,1}:          struct containing the following fields
%             -.algorithm:  selected algorithm in the main-GUI
%                           ['irka'/'rk']
%             -.side:       selection for matched moments in the main-GUI
%                           ['input' / 'output' / 'twosided']
%             -.system:     selected system in the main-GUI
%                           (instance of class sss)
%           {1,2}:          cell array containing the entries of the
%                           input-table of the main-GUI
%           {1,3}:          cell array containing the entries of the
%                           output-table of the main-GUI
% Output Arguments:
%       -varargout{1}:      struct containing the following fields
%           -.inputData:    cell-array containing the new shifts and
%                           directions (input-space)
%           -.outputData:   cell-array containing the new shifts and
%                           directions (output-space)
%           -.hermite:      hermite-interpolation selected or not
%                           [0 / 1]
% ------------------------------------------------------------------
% Authors:      Heiko Panzer, Sylvia Cremer, Maria Cruz Varona, 
%               Alessandro Castagnotto
% Last Change:  07 March 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------


%Initialisation code

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mimoKrylov_OpeningFcn, ...
                   'gui_OutputFcn',  @mimoKrylov_OutputFcn, ...
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


function mimoKrylov_OpeningFcn(hObject, eventdata, handles, varargin)

    %Read out the data passed to the Sub-GUI

    parameter = varargin{1}{1,1};
    tableIn = varargin{1}{1,2};
    tableOut = varargin{1}{1,3};

    %Save the parameters in the GUI-handles

    handles.parameter = parameter;
    guidata(hObject,handles);

    %List all possible matrices that might contain tangential directions
    %(Do it here and not in a create function because the handles are needed to
    %evalutae if a matric is appropriate or not)

    set(handles.pu_importDirections,'String',listDirectionsInWorkspace(handles));

    %Set the GUI-parameters dependent on the passed data

    set(handles.uitable_input,'Data',tableIn);
    set(handles.uitable_output,'Data',tableOut);

    if strcmp(parameter.algorithm,'irka')          %IRKA

        set(handles.rb_hermite,'Value',1);
        set(handles.rb_hermite,'Enable','off');

        set(handles.uitable_input,'Visible','on');
        set(handles.uitable_output,'Visible','off');

    else                                            %RK

        if strcmp(parameter.side,'input')

            set(handles.rb_hermite,'Value',0);
            set(handles.rb_hermite,'Enable','off');

            set(handles.rb_importDirections_input,'Value',1);
            set(handles.rb_importDirections_output,'Enable','off');

            set(handles.uitable_input,'Visible','on');
            set(handles.uitable_output,'Visible','off');

            set(handles.uitable_input,'ColumnWidth',[199,199,199,0]);

        elseif strcmp(parameter.side,'output')

            set(handles.rb_hermite,'Value',0);
            set(handles.rb_hermite,'Enable','off');

            set(handles.rb_importDirections_output,'Value',1);
            set(handles.rb_importDirections_input,'Enable','off');

            set(handles.uitable_input,'Visible','off');
            set(handles.uitable_output,'Visible','on');       

        else        %twosided

            set(handles.rb_hermite,'Value',1);      

            set(handles.uitable_input,'Visible','on');
            set(handles.uitable_output,'Visible','off');

        end
    end
    
    %Update the displays for reduced Order and matched moments
    
    countMatchedMoments(handles);
    
    %Prevent running the output function
    
    uiwait;

function varargout = mimoKrylov_OutputFcn(hObject, eventdata, handles)

    try

        data.inputData = [];
        data.outputData = [];

        try
           data.inputData = handles.dataInput; 
        end

        try
           data.outputData = handles.dataOutput; 
        end

        data.hermite = get(handles.rb_hermite,'Value');

        varargout{1} = data;

        delete(handles.mimoKrylov);
    
    catch ex
        varargout{1} = [];
    end


%Callbacks Tables

function uitable_output_CellSelectionCallback(hObject, eventdata, handles)
%Reads out the selected indices and stores it in a the handles structure

    handles.SelectedIndex = eventdata.Indices;
    guidata(hObject,handles);

function uitable_input_CellSelectionCallback(hObject, eventdata, handles)
%Reads out the selected indices and stores it in a the handles structure

    handles.SelectedIndex = eventdata.Indices;
    guidata(hObject,handles);

function uitable_input_CellEditCallback(hObject, eventdata, handles)

    %Update the displays for reduced Order and matched moments

    countMatchedMoments(handles);

function uitable_output_CellEditCallback(hObject, eventdata, handles)
    
    %Update the displays for reduced Order and matched moments

    countMatchedMoments(handles);

    
%Callbacks Buttons Input/Output

function pb_input_Callback(hObject, eventdata, handles)

    %Set the right table visible

    set(handles.uitable_input,'Visible','on');
    set(handles.uitable_output,'Visible','off');

    %Change the color of the buttons to signal witch table is selected

    set(handles.pb_input,'BackgroundColor',[1;0.843;0]);
    set(handles.pb_output,'BackgroundColor',[0.94;0.94;0.94]);
    
    %Enable the correct radioButton for the direction input
    
    set(handles.rb_importDirections_input,'Value',1);
    set(handles.rb_importDirections_input,'Enable','on');
    set(handles.rb_importDirections_output,'Enable','off');

function pb_output_Callback(hObject, eventdata, handles)

    %Set the right table visible

    set(handles.uitable_input,'Visible','off');
    set(handles.uitable_output,'Visible','on');

    %Change the color of the buttons to signal witch table is selected

    set(handles.pb_output,'BackgroundColor',[1;0.843;0]);
    set(handles.pb_input,'BackgroundColor',[0.94;0.94;0.94]);
    
    %Enable the correct radioButton for the direction input
    
    set(handles.rb_importDirections_output,'Value',1);
    set(handles.rb_importDirections_output,'Enable','on');
    set(handles.rb_importDirections_input,'Enable','off');
    
    
%Callbacks Panel Paramters

function rb_hermite_Callback(hObject, eventdata, handles)

    if get(handles.rb_hermite,'Value') == 1

        set(handles.pb_input,'Visible','off');
        set(handles.pb_output,'Visible','off');

        if get(handles.rb_blockKrylov,'Value') ~= 1

            set(handles.uitable_input,'ColumnWidth',{150 150 150 150});

        end

        set(handles.uitable_input,'Visible','on');
        set(handles.uitable_output,'Visible','off');

        set(handles.rb_importDirections_input,'Enable','on');
        set(handles.rb_importDirections_output,'Enable','on');

    else

        set(handles.pb_input,'Visible','on');
        set(handles.pb_output,'Visible','on');

        if get(handles.rb_blockKrylov,'Value') ~= 1

            set(handles.uitable_input,'ColumnWidth',{199 199 199 0});
        end

        if strcmp(get(handles.uitable_input,'Visible'),'on')
            set(handles.rb_importDirections_input,'Value',1);
            set(handles.rb_importDirections_input,'Enable','on');
            set(handles.rb_importDirections_output,'Enable','off');

            set(handles.pb_input,'BackgroundColor',[1;0.843;0]);
            set(handles.pb_output,'BackgroundColor',[0.94;0.94;0.94]);
        else
            set(handles.pb_output,'BackgroundColor',[1;0.843;0]);
            set(handles.pb_input,'BackgroundColor',[0.94;0.94;0.94]);

            set(handles.rb_importDirections_output,'Value',1);
            set(handles.rb_importDirections_output,'Enable','on');
            set(handles.rb_importDirections_input,'Enable','off');
        end
    end

    %Update the displays for reduced Order and matched moments

    countMatchedMoments(handles);

function bg_blockKrylov_SelectionChangedFcn(hObject, eventdata, handles)

    if get(handles.rb_blockKrylov,'Value') == 1
       
        set(handles.uitable_input,'ColumnWidth',{300,300,0,0});
        set(handles.uitable_output,'ColumnWidth',{300,300,0});
        
        set(handles.panel_importDirections,'Visible','off');
        
    else
        
        if get(handles.rb_hermite,'Value') == 1
            set(handles.uitable_input,'ColumnWidth',{150,150,150,150});
        else
            set(handles.uitable_input,'ColumnWidth',{199,199,199,0});
        end
        set(handles.uitable_output,'ColumnWidth',{199,199,199});
        
        set(handles.panel_importDirections,'Visible','on');
    end

    
%Callbacks Panel EditTable

function pb_addRow_Callback(hObject, eventdata, handles)
% add empty row to table

    if strcmp(get(handles.uitable_input,'Visible'),'on')
        x=get(handles.uitable_input,'Data');
        x{size(x,1)+1,1}=[];
        set(handles.uitable_input,'Data',x)
    else
        x=get(handles.uitable_output,'Data');
        x{size(x,1)+1,1}=[];
        set(handles.uitable_output,'Data',x)
    end

function pb_deleteRow_Callback(hObject, eventdata, handles)

    %Read out the selected indices from the handles

    currentCell = handles.SelectedIndex;
    currentCell = currentCell(:,1);
    indices = unique(currentCell);  

    %Delete selected rows

    if strcmp(get(handles.uitable_input,'Visible'),'on')
        x=get(handles.uitable_input,'Data');

        for i = 1:size(indices,1)
           x(indices(i,1)-i+1,:) = []; 
        end

        set(handles.uitable_input,'Data',x)
    else
        x=get(handles.uitable_output,'Data');

        for i = 1:size(indices,1)
           x(indices(i,1)-i+1,:) = []; 
        end

        set(handles.uitable_output,'Data',x)
    end
    
    %Update the displays for reduced Order and matched moments
    
    countMatchedMoments(handles);

function pb_reset_Callback(hObject, eventdata, handles)

    if strcmp(get(handles.uitable_input,'Visible'),'on')
        set(handles.uitable_input,'Data',cell(1,4))
    else
        set(handles.uitable_output,'Data',cell(1,3))
    end
    
    %Update the displays for reduced Order and matched moments
    
    countMatchedMoments(handles);

function pb_complexConjugate_Callback(hObject, eventdata, handles)

    %Read out the selected indices from the handles

    currentCell = handles.SelectedIndex;
    currentCell = currentCell(:,1);
    indices = unique(currentCell); 

    %Add the complex conjugate of the selected cells as new rows

    if strcmp(get(handles.uitable_input,'Visible'),'on')
        x=get(handles.uitable_input,'Data');
        xConj = cell(size(indices,1),4);

        for i = 1:size(indices,1)      
           xConj{i,1} = conj(x{indices(i,1),1});
           xConj{i,2} = x{indices(i,1),2};
           xConj{i,3} = x{indices(i,1),3};
           xConj{i,4} = x{indices(i,1),4};
        end

        x = [x;xConj];
        set(handles.uitable_input,'Data',x);
    else
        x=get(handles.uitable_output,'Data');
        xConj = cell(size(indices,1),3);

        for i = 1:size(indices,1)      
           xConj{i,1} = conj(x{indices(i,1),1});
           xConj{i,2} = x{indices(i,1),2};
           xConj{i,3} = x{indices(i,1),3};
        end

        x = [x;xConj];
        set(handles.uitable_output,'Data',x);
    end
    
    %Update the displays for reduced Order and matched moments
    
    countMatchedMoments(handles);
    
    
%Callbacks Panel ImportShifts

function pb_importShifts_Callback(hObject, eventdata, handles)
%If a vector from workspace was selected, add it to table

    %Read out the vector from workspace

    contents = cellstr(get(handles.pu_importShifts,'String'));
    try
        x=evalin('base',contents{get(handles.pu_importShifts,'Value')});
    catch ex
       errordlg('No valid matrix selected');
       return;
    end
    
    writeShiftsToTable(handles,x);
    
    %Update the displays for reduced Order and matched moments
    
    countMatchedMoments(handles);

function et_importShifts_Callback(hObject, eventdata, handles)

    %Execute the typed in command in the workspace

    text = get(handles.et_importShifts,'String');

    try
        m = evalin('base',text);
    catch er
        errordlg(er.message);
        return;
    end

    %Read out the expension points from the vector

    if (size(m,1) == 1 || size(m,2) == 1) || isCorrectExpensionPointMatrix(m) == 1
        writeShiftsToTable(handles,m);
    else
        errordlg('The matrix created with the command line has not the correct format.');
        return;
    end
    
    %Update the displays for reduced Order and matched moments
    
    countMatchedMoments(handles);

function pu_importShifts_CreateFcn(hObject, eventdata, handles)
%List of vectors in workspace that might be expansion points

    set(hObject,'String',listShiftsInWorkspace)


%Callbacks Panel ImportDirectinos

function et_importDirections_Callback(hObject, eventdata, handles)
%Execute the typed in command in the workspace

    text = get(handles.et_importDirections,'String');

    try
        m = evalin('base',text);
    catch er
        errordlg(er.message);
        return;
    end

    %Read out the expension points from the vector

    if size(m,1) == handles.parameter.system.p || size(m,2) == handles.parameter.system.p || ...
            size(m,1) == handles.parameter.system.m || size(m,2) == handles.parameter.system.m
        
        writeDirectionsToTable(handles,m);
    else
        errordlg('The matrix created with the command line has not the correct format.');
        return;
    end
    
    %Update the displays for reduced Order and matched moments
    
    countMatchedMoments(handles);

function pb_importDirections_Callback(hObject, eventdata, handles)
%If a matrix from workspace was selected, add it to table

    %Read out the matrix from workspace

    contents = cellstr(get(handles.pu_importDirections,'String'));
    try
        x=evalin('base',contents{get(handles.pu_importDirections,'Value')});
    catch ex
       errordlg('No valid matrix selected');
       return;
    end
    
    writeDirectionsToTable(handles,x);
    
    %Update the displays for reduced Order and matched moments
    
    countMatchedMoments(handles);


%Callback Refresh

function pb_refresh_Callback(hObject, eventdata, handles)
%Update the lists for possible expension-point-vectors and possible
%tangential direction Matrices

    set(handles.pu_importShifts,'String',listShiftsInWorkspace)
    set(handles.pu_importDirections,'String',listDirectionsInWorkspace(handles))

    
%Callback Ok

function pb_ok_Callback(hObject, eventdata, handles)
%Check if the table-entries are correct and close the GUI if that is
%correct

    if ~dataCorrect(handles)
       return;        
    else
       uiresume;
       return; 
    end
  
    
%Auxiliary functions

function x = listShiftsInWorkspace
% list all vectors from workspace that might be expansion points

    s=evalin('base', 'whos');
    
    % preallocate memory
    x=cell(length(s),1);
    
    for i=1:length(s)
        
        if strcmp(s(i).class,'double') && length(s(i).size)==2 && ...
                (s(i).size(1)==2 || s(i).size(2)==2)    %Matrix

            %check wheather one dimension is integer
            vec = evalin('base',s(i).name);

            % save name

            if isCorrectExpensionPointMatrix(vec) == 1
                x{i}=s(i).name;
            end

        elseif strcmp(s(i).class,'double') && length(s(i).size)==2 && ...
                (s(i).size(1)==1 || s(i).size(2)==1)    %Vector

            x{i}=s(i).name;

        end
    end
    % remove empty (non-system) entries
    x(cellfun(@isempty,x)) = [];
    if isempty(x)
       x = {''}; 
    end

function x = listDirectionsInWorkspace(handles)
% list all matrices from workspace that might contain possible tangential
% directions

    s=evalin('base', 'whos');
    
    % preallocate memory
    x=cell(length(s),1);
    
    for i=1:length(s)
        
        if strcmp(s(i).class,'double') && length(s(i).size)==2 && ...
           (s(i).size(1) == handles.parameter.system.m || ...
           s(i).size(2) == handles.parameter.system.m || ...
           s(i).size(1) == handles.parameter.system.p || ...
           s(i).size(2) == handles.parameter.system.p)

            x{i}=s(i).name;

        end
    end
    % remove empty (non-system) entries
    x(cellfun(@isempty,x)) = [];
    if isempty(x)
       x = {''}; 
    end    
    
    
function writeShiftsToTable(handles,m)
%Takes a matrix m with expension points and creates a cell array that can
%be used as Data for an uitable out of the values of m  

    %Create a matrix in which the expension-points are stored in the
    %correct format

    if size(m,1) == 1 || size(m,2) == 1     %only expension points
        
        if size(m,1) == 1
            m = m';
        end
        
        values = unique(m);
        
        table = zeros(length(values),2);
        
        for i = 1:size(table,1)
           table(i,1) = values(i);
           xTemp = find(m == values(i));
           table(i,2) = length(xTemp);
        end
        
        m = table;
        
    elseif size(m,2)~=2
        m=m';
    end
    
    if size(m,2)==2 && size(m,1)==2        %expension points with #moments
        
       if ~(mod(m(1,2),1)==0 && mod(m(2,2),1)==0)
           
          if mod(m(1,1),1)==0 && mod(m(2,1),1)==0
              xTemp = m;
              m(1,1) = xTemp(1,2);
              m(2,1) = xTemp(2,2);
              m(1,2) = xTemp(1,1);
              m(2,2) = xTemp(2,1);
          elseif mod(m(1,1),1)==0 && mod(m(1,2),1)==0
              xTemp = m;
              m(1,1) = xTemp(2,1);
              m(2,1) = xTemp(2,2);
              m(1,2) = xTemp(1,1);
              m(2,2) = xTemp(1,2);
          elseif mod(m(2,1),1)==0 && mod(m(2,2),1)==0
              xTemp = m;
              m(1,2) = xTemp(2,1);
              m(2,1) = xTemp(1,2);
          end
       end        
    end
    
    %Create Cell-Array
    
    tableData=num2cell(m);
    
    %Write the data to the table
 
    if strcmp(get(handles.uitable_input,'Visible'),'on')
        
        momData = get(handles.uitable_input,'Data');
        
        if size(tableData,1) > size(momData,1)
           momData = [momData;cell(size(tableData,1)-size(momData,1),4)];
           momData(:,1:2) = tableData(:,1:2);
        elseif size(tableData,1) < size(momData,1)
           momData = [tableData, momData(1:size(tableData,1),3:4)];
        else
           momData = [tableData, momData(:,3:4)];
        end
        
        set(handles.uitable_input,'Data',momData);
    else
        momData = get(handles.uitable_output,'Data');
        
        if size(tableData,1) > size(momData,1)
           momData = [momData;cell(size(tableData,1)-size(momData,1),4)];
           momData(:,1:2) = tableData(:,1:2);
        elseif size(tableData,1) < size(momData,1)
           momData = [tableData, momData(1:size(tableData,1),3:4)];
        else
           momData = [tableData, momData(:,3:4)];
        end
        
        set(handles.uitable_output,'Data',momData);
    end
    
function writeDirectionsToTable(handles,m)
%Takes a matrix m with tangentail directions and creates a cell array that 
%can be used as Data for an uitable out of the values of m

    if get(handles.rb_importDirections_input,'Value') == 1      %Input
       
        %Write the tangential directions from the matrix m to a cell array
        %witch has the correct format
        
        if size(m,1) == handles.parameter.system.m
            
            data = cell(size(m,2),1);
            for i = 1:size(m,2)
               data{i,1} = vec2string(m(:,i)); 
            end
            
        elseif size(m,2) == handles.parameter.system.m
            
            data = cell(size(m,1),1);
            for i = 1:size(m,1)
               data{i,1} = vec2string(m(i,:)); 
            end
            
        else
            error('Wrong matrix dimension for tangential directions for input subspace!');
        end  
        
        %Add the tangential directions to the table
        
        momData = get(handles.uitable_input,'Data');
        
        if size(momData,1) > size(data,1)
            
            momData(:,3) = cell(size(momData,1),1);
            momData(1:size(data,1),3) = data(:,1);
           
        elseif size(momData,1) < size(data,1)
            
            momData = [momData;cell(size(data,1)-size(momData,1),4)];
            momData(:,3) = data(:,1);
            
        else
            
            momData(:,3) = data(:,1);
             
        end
        
        set(handles.uitable_input,'Data',momData);
        
    else                                                        %Output
        
        %Write the tangential directions from the matrix m to a cell array
        %witch has the correct format       
        
        if size(m,1) == handles.parameter.system.p
            
            data = cell(size(m,2),1);
            for i = 1:size(m,2)
               data{i,1} = vec2string(m(:,i)); 
            end
            
        elseif size(m,2) == handles.parameter.system.p
            
            data = cell(size(m,1),1);
            for i = 1:size(m,1)
               data{i,1} = vec2string(m(i,:)); 
            end
            
        else
            error('Wrong matrix dimension for tangential directions for input subspace!');
        end 
        
        %Add the tangential directions to the table
        
        if strcmp(get(handles.uitable_input,'Visible'),'on')    %Input-Tanble
            
            momData = get(handles.uitable_input,'Data');
            
            if size(momData,1) > size(data,1)
            
                momData(:,4) = cell(size(momData,1),1);
                momData(1:size(data,1),4) = data(:,1);

            elseif size(momData,1) < size(data,1)

                momData = [momData;cell(size(data,1)-size(momData,1),4)];
                momData(:,4) = data(:,1);

            else

                momData(:,4) = data(:,1);

            end

            set(handles.uitable_input,'Data',momData);

        else                                                %Output-Table
            
            momData = get(handles.uitable_output,'Data');
        
            if size(momData,1) > size(data,1)

                momData(:,3) = cell(size(momData,1),1);
                momData(1:size(data,1),3) = data(:,1);

            elseif size(momData,1) < size(data,1)

                momData = [momData;cell(size(data,1)-size(momData,1),3)];
                momData(:,3) = data(:,1);

            else

                momData(:,3) = data(:,1);

            end

            set(handles.uitable_output,'Data',momData);
            
        end
        
        
    end
    
    
function [] = countMatchedMoments(handles)
%Updates the number of matched moemts if for MOR with Krylov

    try
        %Read out the data from the tables for input and output shifts
        
        x=get(handles.uitable_input,'Data');
        data1=cell2mat(x(:,1:2));
        x = get(handles.uitable_output,'Data');
        data2=cell2mat(x(:,1:2));
        
        if ~isempty(data1)                  %Input
           data1 = sum(data1(:,2)); 
        else
           data1 = 0;
        end

        if ~isempty(data2)                  %Output
           data2 = sum(data2(:,2)); 
        else
           data2 = 0;
        end
        
        %Updata the displays of reduced order and matched moments
        
        redOrder = 0;

        if strcmp(handles.parameter.side,'twosided')     %Two sided
            if strcmp(handles.parameter.algorithm,'rk') && ...
                    get(handles.rb_hermite,'Value') == 0 %Input and Output specifiable

                set(handles.st_matchedMom,'String',data2+data1);
            else
                set(handles.st_matchedMom,'String',2*data1);
            end
            
            redOrder = data1;
            
        elseif strcmp(handles.parameter.side,'input')    %Input
            set(handles.st_matchedMom,'String',data1);
            redOrder = data1;
        else                                             %Output
            set(handles.st_matchedMom,'String',data2);
            redOrder = data2;
        end
        
        set(handles.st_redOrder,'String',redOrder);

        %Check weather the reduced order is bigger than the original order


       if handles.parameter.system.n < redOrder
          msgbox('Reduced order is bigger than the original order of the system. Please correct that before reducing the system.','Warning','Warn');
          uiwait;
       end
    end
    
    
function t=isCorrectExpensionPointMatrix(m)
%Checks if a matric has the korrekt format to contain possible expension
%points for the Krylov-methods

    t = 0;

    if size(m,1)==2 || size(m,2)==2    
        %check wheather one dimension is integer
        isInt1 = 1;
        isInt2 = 1;
        if size(m,1) == 2 && size(m,2) == 2
            isInt1 = 0;
            isInt2 = 0;
            if (mod(m(1,1),1)==0 && mod(m(2,1),1)==0) || ...
               (mod(m(1,2),1)==0 && mod(m(2,2),1)==0)
           
                isInt1 = 1;
            end
            if (mod(m(1,1),1)==0 && mod(m(1,2),1)==0) || ...
               (mod(m(2,1),1)==0 && mod(m(2,2),1)==0)
           
                isInt2 = 1;
            end
     
        elseif size(m,1) == 2
            for j = 1:size(m,2)
               if mod(m(1,j),1)~=0
                  isInt1 = 0; 
               end
               if mod(m(2,j),1)~=0
                  isInt2 = 0; 
               end
            end
            
        else
            for j = 1:size(m,1)
               if mod(m(j,1),1)~=0
                  isInt1 = 0; 
               end
               if mod(m(j,2),1)~=0
                  isInt2 = 0; 
               end
            end
        end
               
        if isInt1 == 1 || isInt2 == 1
            t = 1;
        end
    end
    
function [string] = vec2string(vect)
%Gives a string-representation of a vector

    string = '[';
    for i = 1:length(vect)-1
        string = strcat(strcat(string,num2str(vect(i))),';');
    end
    string = strcat(strcat(string,num2str(vect(length(vect)))),']');
 
    
function t = dataCorrect(handles)
%Checks if all the entrys in the table are correct

    t = 1;

    try
    
        if strcmp(handles.parameter.algorithm,'rk')     %Rational Krylov

            if strcmp(handles.parameter.side,'twosided')    %Twosided

                if get(handles.rb_hermite,'Value') == 1         %Hermite

                    if get(handles.rb_blockKrylov,'Value') == 0  %Tangential

                        data = get(handles.uitable_input,'Data');
                        tableContentCorrect(data,handles.parameter.system.m,handles.parameter.system.p); 
                        handles.dataInput = data;

                    else                                        %Block Krylov

                        data = get(handles.uitable_input,'Data');
                        tableContentCorrect(data(:,1:2),0,0);
                        handles.dataInput = data(:,1:2);
                 
                    end

                else                                %No Hermite

                    if get(handles.rb_blockKrylov,'Value') == 0  %Tangential

                        data = get(handles.uitable_input,'Data');
                        tableContentCorrect(data(:,1:3),handles.parameter.system.m,0); 
                        handles.dataInput = data(:,1:3);
                        
                        data = get(handles.uitable_output,'Data');
                        tableContentCorrect(data,handles.parameter.system.p,0); 
                        handles.dataOutput = data;
                        
                    else                                        %Block Krylov

                        data = get(handles.uitable_input,'Data');
                        tableContentCorrect(data(:,1:2),0,0);
                        handles.dataInput = data(:,1:2);
                        
                        data = get(handles.uitable_output,'Data');
                        tableContentCorrect(data(:,1:2),0,0);                         
                        handles.dataOutput = data(:,1:2);    
                        
                    end
                end

            elseif strcmp(handles.parameter.side,'input')       %Input

                if get(handles.rb_blockKrylov,'Value') == 0  %Tangential

                    data = get(handles.uitable_input,'Data');
                    tableContentCorrect(data(:,1:3),handles.parameter.system.m,0);
                    handles.dataInput = data(:,1:3);
                    
                else                                %Block Krylov

                    data = get(handles.uitable_input,'Data');
                    tableContentCorrect(data(:,1:2),0,0);
                    handles.dataInput = data(:,1:2);
                    
                end

            else                                            %Output
                
                if get(handles.rb_blockKrylov,'Value') == 0  %Tangential

                    data = get(handles.uitable_output,'Data');
                    tableContentCorrect(data(:,1:3),handles.parameter.system.p,0);
                    handles.dataOutput = data(:,1:3);
                    
                else                                %Block Krylov

                    data = get(handles.uitable_output,'Data');
                    tableContentCorrect(data(:,1:2),0,0);
                    handles.dataOutput = data(:,1:2);
                    
                end               
            end

        else                                            %IRKA

            if get(handles.rb_blockKrylov,'Value') == 0  %Tangential

                data = get(handles.uitable_input,'Data');
                tableContentCorrect(data,handles.parameter.system.m,handles.parameter.system.p); 
                handles.dataInput = data;
                
            else                                        %Block Krylov

                data = get(handles.uitable_input,'Data');
                tableContentCorrect(data(:,1:2),0,0); 
                handles.dataInput = data(:,1:2);
                
            end          
        end
        
        %Store data in handles
        
        guidata(handles.uitable_input,handles);

    catch ex
        errordlg(ex.message);
        t = 0;
        return;
    end   
    
function tableContentCorrect(c,m3,m4)
%Check if the content of the cell-array is valid to be used for MOR

    if size(c,2) == 2
        
        if arrayEmpty(c)
            error('At least one entry of the table is empty. Please correct this first!');
        end

        if ~matrixIsNumeric(c(:,1))
            error('All shifts have to be numeric values. Please correct this first!');
        end

        if ~matrixIsInteger(c(:,2))
            error('The number of moments has to be a integer-value for each shift. Please correct this first!');
        end
        
    elseif size(c,2) == 3
        
        if arrayEmpty(c)
            error('At least one entry of the table is empty. Please correct this first!');
        end

        if ~matrixIsNumeric(c(:,1))
            error('All shifts have to be numeric values. Please correct this first!');
        end

        if ~matrixIsInteger(c(:,2))
            error('The number of moments has to be a integer-value for each shift. Please correct this first!');
        end
        
        if ~correctTangentialDirection(c(:,3),m3)
            error('Not all specified tangential directions are compartiple with the current system. Please correct this first!');
        end
        
    elseif size(c,2) == 4
        
        if arrayEmpty(c)
            error('At least one entry of the table is empty. Please correct this first!');
        end

        if ~matrixIsNumeric(c(:,1))
            error('All shifts have to be numeric values. Please correct this first!');
        end

        if ~matrixIsInteger(c(:,2))
            error('The number of moments has to be a integer-value for each shift. Please correct this first!');
        end
        
        if ~correctTangentialDirection(c(:,3),m3)
            error('Not all specified tangential directions are compartiple with the current system. Please correct this first!');
        end
        
        if ~correctTangentialDirection(c(:,4),m4)
            error('Not all specified tangential directions are compartiple with the current system. Please correct this first!');
        end
        
    end
        
             
function t = arrayEmpty(c)
%Check if a cell array has an empty cell
    
    t = 0;

    for i = 1:size(c,1)
       for j = 1:size(c,2)
          if isempty(c{i,j})
             t = 1;
             break;
          end
       end
    end
        
function t = matrixIsInteger(c)
%Check if all entrys of an cell-array are integer values
   
    t = 1;

    for i = 1:size(c,1)
       for j = 1:size(c,2)
          if mod(c{i,j},1) ~= 0
             t = 0;
             break;
          end
       end
    end
        
function t = matrixIsNumeric(c)
%Check if all entrys of an cell-array are numeric values

    t = 1;
    
    for i = 1:size(c,1)
       for j = 1:size(c,2)
          if isnumeric(c{i,j}) == 0
             t = 0;
             break;
          end
       end
    end
    
function t = correctTangentialDirection(c,dimension)
%Check if all entrys of an cell-array can be evaluated to an vector of the
%size dimension

    t = 1;
    
    for i = 1:size(c,1)
       for j = 1:size(c,2)
          vec = evalin('base',c{i,j});
          
          if size(vec,1) ~= dimension
             t = 0;
             break;
          end
       end
    end
    
