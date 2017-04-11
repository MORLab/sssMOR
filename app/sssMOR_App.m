function varargout = sssMOR_App(varargin)
% SSSMOR_APP - sssMOR toolbox App 
%
% Syntax:
%       SSSMOR_APP
%
% Description:
%       The *sssMOR App* is a Graphical User Interface for Model Order Reduction,
%       which uses the functions from the *sss Toolbox* and from the *sssMOR Toolbox*.  
%
%       The main menu of the user interface contains the five items *About*,
%       *Loading and Setting up Models*, *Model Order Reduction*, *Postprocessing
%       and Visualization* and *System Analysis*. The sequence of these five items
%       follows the usual workflow for Model Order Reduction, which starts with
%       loading the desired model, continues with the reduction and ends with an
%       evaluation of the results.
%
% See Also:
%       sss_gettingStarted, sssMOR_gettingStarted
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Sylvia Cremer, Niklas Kochdumper,
%               Maria Cruz Varona, Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 May 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


%Structure:
%
%   INITIALIZATION
%
%   FIGURE(TOP-LEVEL)
%   
%       MAIN-MENUE
%
%           ABOUT
%
%           POSTPROCESSING AND VISUALISATION
%
%           LOADING AND SETTING UP MODELS
%
%           MODEL ORDER REDUCTION
%
%           SYSTEM ANALYSIS
%
%       FOOTER
%
%   UTILITY FUNCTIONS

%--------------------------------------------------------------------------
%                           INITIALIZATION
%--------------------------------------------------------------------------


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sssMOR_App_OpeningFcn, ...
                   'gui_OutputFcn',  @sssMOR_App_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
% try
    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
% catch ex
%     errordlg(ex.message,'Error Dialog','modal')
%     throw(ex)   %***
% end
        
        
function sssMOR_App_OpeningFcn(hObject, eventdata, handles, varargin)  %#ok<*INUSL>

    %Make latex code possible for all static text-fields whose tag names start
    %with "latex"  
    
    lbls = findall(hObject);
    
    for i=1:length(lbls)
        
          l = lbls(i);
          t = get(l,'tag');
          
          if strfind(t,'latex_') == 1
              
             % Get current text, position and tag
             t = 'nothing'; 
             set(l,'units','characters');
             s = get(l,'string');
             p = get(l,'position');
             parent = get(l,'parent');
             
             % Remove the UICONTROL
             delete(l);
             
             % Replace it with an axis with an text object on it
             handles.(t) = axes('parent',parent,'units','characters','position',[p(1) p(2) p(3) p(4)],'visible','off');
             handles.(t) = text(0,0.3,s,'interpreter','latex');            
          end
    end
    
    %Add the subfolder from the GUI-folder-structure to the search-paths
    
    guiPath = which('sssMOR_App.m');
    seperator = filesep;
    pathArray = strsplit(guiPath,seperator);
    path = '';
    
    for i = 1:(size(pathArray,2)-1)
       
        if i > 1
           path = strcat(path,seperator,pathArray{1,i});
        else
           path = strcat(path,pathArray{1,i}); 
        end       
    end
    
    %Add path with all subfolders
    
    addpath(genpath(path));
    
    %Save the Width of all Tables in pixels in the handles-structure and
    %correctly set the width of the columns of all tables
    %(Nescesarray because the column width has to be in 'pixels', but because
    % the hole GUI is in 'characters' this leads to visualisation errors 
    % on different operating systems)
    
    p = get(handles.figure1,'Position');
    widthChar = p(1,3);
    set(handles.figure1,'Units','pixels');
    p = get(handles.figure1,'Position');
    widthPix = p(1,3);
    charToPix = widthPix/widthChar;
    
    w = charToPix*74;
    set(handles.uitable_mor_krylov_MimoExps,'ColumnWidth',{round(w/4),round(w/4),round(w/4),round(w/4)});
    handles.widthTableMimoKrylovInput = w;
    
    w = charToPix*74;
    set(handles.uitable_mor_krylov_MimoExps_output,'ColumnWidth',{round(w/3),round(w/3),round(w/3)});
    handles.widthTableMimoKrylovOutput = w;
    
    w = charToPix*40;
    set(handles.uitable_mor_krylov,'ColumnWidth',{round(w/2),round(w/2)});
    
    w = charToPix*40;
    set(handles.uitable_mor_krylov_output,'ColumnWidth',{round(w/2),round(w/2)});
    
    
    %Add the pictures for Header and footer
    
    A = imread('Heading.png');
    %Scale image to the correct size
    set(handles.logo_tum,'Units','pixels');
    p = get(handles.logo_tum,'Position');
    set(handles.logo_tum,'Units','characters');
    axes(handles.logo_tum);
    %Show image and set properties
    h=image(imresize(A,[round(p(1,4)),round(p(1,3))]));
    set(h,'ButtonDownFcn',@(hObject,eventdata)sssMOR_App('logo_tum_ButtonDownFcn',handles.logo_tum,eventdata,guidata(hObject)));
    set(handles.logo_tum,'XTick',[]);
    set(handles.logo_tum,'YTick',[]);
    
    A=imread('Footer.png');
    %Scale image to the correct size
    set(handles.logos_footer,'Units','pixels');
    p = get(handles.logos_footer,'Position');
    set(handles.logos_footer,'Units','characters');
    axes(handles.logos_footer);
    %Show image and set properties
    h=image(imresize(A,[round(p(1,4)),round(p(1,3))]));
    set(h,'ButtonDownFcn',@(hObject,eventdata)sssMOR_App('logos_footer_ButtonDownFcn',handles.logos_footer,eventdata,guidata(hObject)));
    set(handles.logos_footer,'XTick',[]);
    set(handles.logos_footer,'YTick',[]);

    
    %Set default-values for the variables saved in handles
    
    handles.splash = splash('splash.png');
    handles.output = hObject;
    handles.virtgr_an_red_buttons=[handles.pb_an_sys1_stability,handles.pb_an_sys1_h2,handles.pb_an_sys1_hinf,...
                                   handles.pb_an_sys1_decaytime,handles.pb_an_sys1_dissipativity...
                                   handles.pb_an_sys2_stability,handles.pb_an_sys2_h2,handles.pb_an_sys2_hinf,...
                                   handles.pb_an_sys2_decaytime,handles.pb_an_sys2_dissipativity,...
                                   handles.pb_an_compare_h2, handles.pb_an_compare_hinf,...
                                   handles.pb_an_sys1_calcall,handles.pb_an_sys2_calcall];
    handles.allbuttons=[handles.pb_save,handles.pb_load,handles.pb_create,handles.pb_readascii];
    handles.zoom=[];
    handles.plotData = [];
    handles.chosenSystems = 0;
    handles.storedHsv = {};
    
    %Set factors to convert from character-units to pixel units
    
    set(handles.figure1,'Units','characters');
    posChar = get(handles.figure1,'Position');
    set(handles.figure1,'Units','pixels');
    posPix = get(handles.figure1,'Position');
    
    handles.PixToCharWidth = posChar(1,3)/posPix(1,3);
    handles.PixToCharHeight = posChar(1,4)/posPix(1,4);
    
    %Set the default-folder for opening data-files
    
    handles.letzterpfad='*';
    
    path = mfilename('fullpath');
    pathDirectories = strsplit(path,seperator);
    newPath = '';    
    
    for i = 1:length(pathDirectories)-2
       newPath = strcat(newPath,strcat(pathDirectories{1,i},seperator));
    end
    
    newPath = strcat(newPath,'src',seperator,'sss',seperator,'benchmarks',seperator);
    
    if exist(newPath,'dir') == 7
        %Path exists
        handles.letzterpfad = newPath;
        
    else
        %Path don't exists: Try to find the path in the matlab search path.
        %This is important if the GUI is installed as an App
        
        [newPath,name,ending] = fileparts(which('CDplayer.mat'));
        newPath = strcat(newPath,seperator);
        
        if exist(newPath,'dir') == 7
           handles.letzterpfad = newPath; 
        end        
    end
    
    %Set the footer visible or not dependent on the size of the screen
    %(to convert from pixel units to chararcter units the values are
    %multiplied with 0.198 for horizontal dimensions and with 0.0747 for
    %vertical dimensions)
    
    screensize = get( 0, 'Screensize' );
    
    set(handles.figure1,'Units','pixels');
    position = get(handles.figure1,'Position');
    position(1,1) = max(round((screensize(1,3)-position(1,3))/2),0);
    position(1,2) = max(round((screensize(1,4)-position(1,4))/2),0);
    set(handles.figure1,'Position',position);
    set(handles.figure1,'Units','characters');
    
    if screensize(1,4) < position(1,4)+100    %Footer ausblenden
        
        position = get(handles.logos_footer,'Position');
        diff = position(1,4);
        position(1,2) = position(1,2) - diff;
        set(handles.logos_footer,'Position',position);
        
        position = get(handles.figure1,'Position');
        position(1,4) = position(1,4) - diff;
        set(handles.figure1,'Position',position);
        
        position = get(handles.pb_loading,'Position');
        position(1,2) = position(1,2) - diff;
        set(handles.pb_loading,'Position',position);
        
        position = get(handles.pb_postprocessing,'Position');
        position(1,2) = position(1,2) - diff;
        set(handles.pb_postprocessing,'Position',position);
        
        position = get(handles.pb_analysis,'Position');
        position(1,2) = position(1,2) - diff;
        set(handles.pb_analysis,'Position',position);
        
        position = get(handles.pb_about,'Position');
        position(1,2) = position(1,2) - diff;
        set(handles.pb_about,'Position',position);
        
        position = get(handles.pb_mor,'Position');
        position(1,2) = position(1,2) - diff;
        set(handles.pb_mor,'Position',position);
        
        position = get(handles.panel_load,'Position');
        position(1,2) = position(1,2) - diff;
        set(handles.panel_load,'Position',position);
        
        position = get(handles.panel_post,'Position');
        position(1,2) = position(1,2) - diff;
        set(handles.panel_post,'Position',position);
        
        position = get(handles.panel_about,'Position');
        position(1,2) = position(1,2) - diff;
        set(handles.panel_about,'Position',position);
        
        position = get(handles.panel_mor,'Position');
        position(1,2) = position(1,2) - diff;
        set(handles.panel_mor,'Position',position);
        
        position = get(handles.panel_analysis,'Position');
        position(1,2) = position(1,2) - diff;
        set(handles.panel_analysis,'Position',position);
        
        position = get(handles.logo_tum,'Position');
        position(1,2) = position(1,2) - diff;
        set(handles.logo_tum,'Position',position);
        
        %Set the header visible or not dependent on the size of the screen
        
        set(handles.figure1,'Units','pixels');
        position = get(handles.figure1,'Position');
        set(handles.figure1,'Units','characters');
        
        if screensize(1,4) < position(1,4)+50    %Header ausblenden
            
            position = get(handles.logo_tum,'Position');
            diff = position(1,4);
            
            position = get(handles.pb_loading,'Position');
            diff = diff - position(1,4);
            
            position = get(handles.figure1,'Position');
            position(1,4) = position(1,4) - diff;
            set(handles.figure1,'Position',position);
                  
        end       
    end
    

    % Update handles structure
    guidata(hObject, handles);
    t=timer('TimerFcn',@(t,y)timer_Callback(handles),'ExecutionMode','singleShot','StartDelay',2);
    start(t)


    % obj = fieldnames(handles);
    % for i=1:length(obj)
    %     h=obj(i); h=h{1};
    %     h = eval(['handles.' h]);
    %     if ishandle(h)
    %         get(h, 'Type')
    %         if strcmp(get(h, 'Type'), 'figure')
    %             set(h, 'Color', [240,240,240]/255);
    %         elseif any(strcmp(get(h, 'Type'), 'uicontrol')) || any(strcmp(get(h, 'Type'), 'uipanel')) || any(strcmp(get(h, 'Type'), 'text'))
    %             set(h, 'BackgroundColor', [240,240,240]/255);
    %         end
    %     end
    % end

    %set(handles.figure1,'Color',[240,240,240]/255)

function timer_Callback(handles)
    % removes the splashscreen after given time
    splash(handles.splash,'off')
    
    
    % Display an error message if footer and header are already
    % hidden, but the app-window is still to small to fit on the 
    % screen
    screensize = get( 0, 'Screensize' );
    
    set(handles.figure1,'Units','pixels');
    position = get(handles.figure1,'Position');
    set(handles.figure1,'Units','characters');

    if screensize(1,4) < position(1,4)+20    
        errordlg('The screen is to small to correctly display the app. This will eventually restrict the functionality!');
    end    


function varargout = sssMOR_App_OutputFcn(hObject, eventdata, handles)
    varargout{1} = handles.output;



%--------------------------------------------------------------------------
%                           FIGURE(TOP-LEVEL)
%--------------------------------------------------------------------------

function figure1_CloseRequestFcn(hObject, eventdata, handles)
    % delete temp variables in workspace
    if existInBaseWs('GUI_dont_delete')
        evalin('base','clear GUI_dont_delete')
    end
    % delete GUI
    delete(hObject);

function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
    %Show a hand as the mouse-symbol if the mouse is over a logo with a weblink
    %behind
    %(to convert from pixel units to chararcter units for Windows OS the 
    %values are multiplied with 0.2 for horizontal dimensions and with 
    %0.0769 for vertical dimensions)
    
    p = get(hObject,'CurrentPoint');
    
    p(1,1) = p(1,1)/0.2;
    p(1,2) = p(1,2)/0.0769;
    
    logosVisible = 1;
    set(handles.figure1,'Units','characters');
    position = get(handles.figure1,'Position');
    if position(1,4) < 58       %Footer not shown
       p(1,2) = p(1,2) + 40; 
       if position(1,4) < 51    %Header not shown
           logosVisible = 0;
       end
    end
    
    if logosVisible
        if p(1,1)>925 && p(1,1)<1040 && p(1,2)>700 && p(1,2)<750
            %TUM-Logo (Header)
            set(gcf,'Pointer','hand'); 
        elseif p(1,1)> 5 && p(1,1)<75 && p(1,2)>700 && p(1,2)<750 
            %Lehrstuhl-Logo (Header)
            set(gcf,'Pointer','hand');
        elseif p(1,1)>5 && p(1,1)<80 && p(1,2)>0 && p(1,2)<40
            %sssMOR-Logo (Footer)
            set(gcf,'Pointer','hand');
        elseif p(1,1)>990 && p(1,1)<1040 && p(1,2)>0 && p(1,2)<40
            %TUM-Logo (Footer)
            set(gcf,'Pointer','hand');
        elseif p(1,1)>305 && p(1,1)<420 && p(1,2)>0 && p(1,2)<40
            %MorLab-Logo (Footer)
            set(gcf,'Pointer','hand');
        elseif p(1,1)>650 && p(1,1)<765 && p(1,2)>0 && p(1,2)<40
            %Lehrstuhl-Logo (Footer)
            set(gcf,'Pointer','hand');
        else
            %Not over any logo
            set(gcf,'Pointer','arrow');
        end   
    end




%--------------------------------------------------------------------------
%                             MAIN-MENUE
%--------------------------------------------------------------------------

%Callbacks of the panel

function logo_tum_ButtonDownFcn(hObject, eventdata, handles)
    % link to web page
    %(to convert from pixel units to chararcter units for Windows OS the 
    %values are multiplied with 0.2 for horizontal dimensions and with 
    %0.0769 for vertical dimensions)
    p=get(hObject,'CurrentPoint');
    p(1,1) = p(1,1)*handles.PixToCharWidth/0.2;
    p(1,2) = p(1,2)*handles.PixToCharHeight/0.0769;
    
    if p(1,1)>925 && p(1,1)<1040 && p(1,2)>5 && p(1,2)<80
        web www.tum.de
    elseif p(1,1)> 5 && p(1,1)<75 && p(1,2)>5 && p(1,2)<80  
        web www.rt.mw.tum.de
    end


%Callbacks of the buttons

function pb_loading_Callback(hObject, eventdata, handles)  %#ok<*DEFNU>
    % make load/save tab visible
    set(handles.panel_load,'Visible','on')
    set(handles.panel_mor,'Visible','off')
    set(handles.panel_post,'Visible','off')
    set(handles.panel_analysis,'Visible','off')
    set(handles.panel_about,'Visible','off')
    % refresh list of systems
    pb_refreshsys_Callback(hObject, eventdata, handles)

function pb_mor_Callback(hObject, eventdata, handles)
    % make mor tab visible
    set(handles.panel_load,'Visible','off')
    set(handles.panel_mor,'Visible','on')
    set(handles.panel_post,'Visible','off')
    set(handles.panel_analysis,'Visible','off')
    set(handles.panel_about,'Visible','off')
    % refresh list of systems
    pb_refreshsys_Callback(hObject, eventdata, handles)

function pb_postprocessing_Callback(hObject, eventdata, handles)
    % make post processing tab visible
    set(handles.panel_load,'Visible','off')
    set(handles.panel_mor,'Visible','off')
    set(handles.panel_post,'Visible','on')
    set(handles.panel_analysis,'Visible','off')
    set(handles.panel_about,'Visible','off')
    % refresh list of open figures and systems
    pb_refreshsys_Callback(hObject, eventdata, handles)
    pb_PaV_refeshObjects_Callback(handles.pb_PaV_refeshObjects, eventdata, handles)

function pb_analysis_Callback(hObject, eventdata, handles)
    % make analysis tab visible
    set(handles.panel_load,'Visible','off')
    set(handles.panel_mor,'Visible','off')
    set(handles.panel_post,'Visible','off')
    set(handles.panel_analysis,'Visible','on')
    set(handles.panel_about,'Visible','off')
    % refresh list of systems
    pb_refreshsys_Callback(hObject, eventdata, handles)

function pb_about_Callback(hObject, eventdata, handles)
    set(handles.panel_load,'Visible','off')
    set(handles.panel_mor,'Visible','off')
    set(handles.panel_post,'Visible','off')
    set(handles.panel_analysis,'Visible','off')
    set(handles.panel_about,'Visible','on')  



%--------------------------------------------------------------------------
%                                 ABOUT
%--------------------------------------------------------------------------

function logo_about_CreateFcn(hObject, eventdata, handles)
    % load About logo
    addRelativePaths();
    A=imread('MOR.jpg');
    image(A);
    set(hObject,'YDir','reverse');

function text_about_weblink_ButtonDownFcn(hObject, eventdata, handles)
    % Open a weblink if the mouse is klicked over the text element
    web https://www.rt.mw.tum.de/?sssMOR

    
    
%--------------------------------------------------------------------------
%                    POSTPROCESSING AND VISUALISATION
%--------------------------------------------------------------------------

function pb_PaV_remove_Callback(hObject, eventdata, handles)
%Remove the selected system from the list

    list = get(handles.lb_PaV_selectedSystems,'String');
    index = get(handles.lb_PaV_selectedSystems,'Value');
    
    if isempty(list)
       set(handles.lb_PaV_selectedSystems,'String','');
       set(handles.lb_PaV_selectedSystems,'Value',1);
       return;
    end
    
    sysName = list{index};
    
    list{index,1} = [];
    list = list(~cellfun(@isempty, list));
    
    if index > length(list)
       set(handles.lb_PaV_selectedSystems,'Value',length(list));
    end
    
    set(handles.lb_PaV_selectedSystems,'String',list);
    
    %Remove the stored data for this system
    
    data = handles.plotData;
    
    for i = 1:length(data)
       if strcmp(data{i,1}.name,sysName)
          data{i,1} = [];
          break;
       end
    end
    
    data = data(~cellfun(@isempty, data));
    handles.plotData = data;
    guidata(hObject,handles);
    
    lb_PaV_selectedSystems_Callback(handles.lb_PaV_selectedSystems, eventdata, handles)

function pb_PaV_move_Callback(hObject, eventdata, handles)
%Move the selected system from the list on the left side to the right side

    %Check if the list is empty and display warning if the case
    
    list = get(handles.lb_PaV_systemsWs,'String');
    
    if isempty(list{1,1})
       errordlg('No system selected');
       uiwait;
       return;
    end
    
    indices = get(handles.lb_PaV_systemsWs,'Value');
    
    for j = 1:size(indices,2)
    
        %Read out the selected system

        index = indices(1,j);

        selectedSystem = list{index,1};

        %Check if the system exists already in the list on the right side

        listSelected = get(handles.lb_PaV_selectedSystems,'String');
        errorOccured = 0;
        
        if ~isempty(listSelected)
           for i = 1:length(listSelected)
              if strcmp(listSelected{i,1},selectedSystem)
                errordlg(strcat('The system "',selectedSystem,'" already exists in the list.'));
                uiwait;
                errorOccured = 1;
                break;
              end
           end
        end
        
        if errorOccured == 1
           continue; 
        end

        %Get system from workspace
        
        try
            sys = getSysFromWs(selectedSystem);
        catch ex
            errordlg(ex.message,'Error Dialog','modal')
            uiwait
            return
        end

        %Add the system to the list on right side

        if isempty(listSelected)
            set(handles.lb_PaV_selectedSystems,'String',{selectedSystem});
            set(handles.lb_PaV_selectedSystems,'Value',1);
            if get(handles.rb_PaV_plotStyle_manual,'Value')
                set(handles.panel_PaV_plotStyle,'Visible','on');
            else
                set(handles.panel_PaV_plotStyle,'Visible','off');
            end
            set(handles.cb_PaV_SaveData,'Visible','on');
            set(handles.et_PaV_saveData,'Visible','on');
            set(handles.bg_PaV_plotStyle,'Visible','on');
            set(handles.bg_PaV_Resolution,'Visible','on');
        else
            if length(listSelected) >= 10
                errordlg('The maximal number of selected systems is 10. Please remove a system before you add another one.','Error Dialog','modal');
                uiwait
                return
            end
            listSelected{end+1,1} = selectedSystem;
            set(handles.lb_PaV_selectedSystems,'String',listSelected);
        end

        %Set the default-values for all plot-options for this system

        data.name = selectedSystem;
        data.isSystem = 1;

        p = size(sys.B,2);
        m = size(sys.C,1);
        
        if p > 1|| m > 1
           data.in = 'all';
           data.out = 'all';
           data.sizeInputs = m;
           data.sizeOutputs = p;
           %data.legendText = suggestDefaultLegendText(handles,selectedSystem,1);
        else
           data.in = [];
           data.out = [];
           data.sizeInputs = [];
           data.sizeOutputs = [];
           %data.legendText = suggestDefaultLegendText(handles,selectedSystem,[]);
        end

        data.plotStyle = 'auto';
        data.resolution = 'auto';

        data.save = 1;
        data.variableName = suggestPlotDataName(handles,selectedSystem);
        data.saveTimeseries = 0;
        data.variableNameTimeseries = suggestPlotDataNameTimeseries(handles,selectedSystem);
        data.color = mod(handles.chosenSystems,7)+1;
        data.lineStyle = 1;
        data.markerType = 1;

        data.distribution = 'linear';
        data.min = '';
        data.max = '';
        data.steps = 5000;

        handles.plotData{end+1,1} = data;

        handles.chosenSystems = handles.chosenSystems + 1;

        guidata(hObject,handles);

        %Set the default legend Text

        data.legendText = suggestDefaultLegendText(handles,selectedSystem,[]);
        handles.plotData{end,1} = data;
        guidata(hObject,handles);
    
    end
    
    
    lb_PaV_selectedSystems_Callback(handles.lb_PaV_selectedSystems,eventdata,handles);
    
function pb_PaV_moveObjects_Callback(hObject, eventdata, handles)
%Move the selected object from the list on the left side to the right side

    %Check if the list is empty and display warning if the case
    
    list = get(handles.lb_PaV_objectsWs,'String');
    
    if isempty(list{1,1})
       errordlg('No system selected');
       uiwait;
       return;
    end
    
    indices = get(handles.lb_PaV_objectsWs,'Value');
    
    for j = 1:size(indices,2)
    
        %Read out the selected system

        index = indices(1,j);
    
        %Read out the selected system

        selectedObject = list{index,1};

        %Check if the system exists already in the list on the right side

        listSelected = get(handles.lb_PaV_selectedSystems,'String');
        errorOccured = 0;
        
        if ~isempty(listSelected)
           for i = 1:length(listSelected)
              if strcmp(listSelected{i,1},selectedObject)
                errordlg(strcat('The object "',selectedObject, '" already exists in the list.'));
                uiwait;
                errorOccured = 1;
                break;
              end
           end
        end
        
        if errorOccured
           continue; 
        end

        %Get system from workspace

        try
            object = evalin('base',selectedObject);
        catch ex
            errordlg(ex.message,'Error Dialog','modal')
            uiwait
            return
        end

        %Add the system to the list on right side

        if isempty(listSelected)
            set(handles.lb_PaV_selectedSystems,'String',{selectedObject});
            set(handles.lb_PaV_selectedSystems,'Value',1);
            if get(handles.rb_PaV_plotStyle_manual,'Value')
                set(handles.panel_PaV_plotStyle,'Visible','on');
            else
                set(handles.panel_PaV_plotStyle,'Visible','off');
            end
            set(handles.cb_PaV_SaveData,'Visible','on');
            set(handles.et_PaV_saveData,'Visible','on');
            set(handles.bg_PaV_plotStyle,'Visible','on');
            set(handles.bg_PaV_Resolution,'Visible','on');
        else
            if length(listSelected) >= 10
                errordlg('The maximal number of selected systems is 10. Please remove a system before you add another one.','Error Dialog','modal');
                uiwait
                return
            end
            listSelected{end+1,1} = selectedObject;
            set(handles.lb_PaV_selectedSystems,'String',listSelected);
        end

        %Set the default-values for all plot-options for this system

        data.name = selectedObject;
        data.isSystem = 0;

        if size(object,1) > 1|| size(object,2) > 1
           data.in = 'all';
           data.out = 'all';
           data.sizeInputs = size(object,2);
           data.sizeOutputs = size(object,1);
           data.legendText = suggestDefaultLegendText(handles,selectedObject,1);
        else
           data.in = [];
           data.out = [];
           data.sizeInputs = [];
           data.sizeOutputs = [];
           data.legendText = suggestDefaultLegendText(handles,selectedObject,[]);
        end

        data.plotStyle = 'auto';
        data.resolution = [];

        data.save = [];
        data.variableName = [];
        data.saveTimeseries = [];
        data.variableNameTimeseries = [];
        data.color = mod(handles.chosenSystems,7)+1;
        data.lineStyle = 1;
        data.markerType = 1;

        data.distribution = [];
        data.min = [];
        data.max = [];
        data.steps = [];

        handles.plotData{end+1,1} = data;

        handles.chosenSystems = handles.chosenSystems + 1;

        guidata(hObject,handles);
    
    end
    
    lb_PaV_selectedSystems_Callback(handles.lb_PaV_selectedSystems,eventdata,handles);
    
function pb_PaV_refeshObjects_Callback(hObject, eventdata, handles)
%Refresh the list of data-objects (frd, tf, etc.) existing in workspace

    %Get all objects existing in workspace
    
    switch get(handles.plot_type,'Value')
        
        case 1      %Bode 
            class = 'frd';
        case 2      %Magitude
            class = 'frd';
        case 3      %Singular Values
            class = 'frd';
        case 4      %Impulse
            class = 'tf';
        case 5      %Step
            class = 'tf';
        case 6      %Pole-Zero-Map
            class = 'zpk';
    end
    
    l = listClassesInWorkspace(class);
    
    %Check if the previously selected object still exists in workspace and
    %select it if the case

    l_alt = get(handles.lb_PaV_objectsWs,'String');

    if ~isempty(l_alt)

        sOld=l_alt{get(handles.lb_PaV_objectsWs,'Value')};
        indexOld=find(strcmp(sOld,l));

        if ~isempty(indexOld)
            set(handles.lb_PaV_objectsWs,'Value',indexOld);
        else
            set(handles.lb_PaV_objectsWs,'Value',1);
        end

    else
      set(handles.lb_PaV_objectsWs,'Value',1);
    end
    
    %If there are no instances of the given class in workspace, set the 
    %list of selectable options to the empty string, else set it to the 
    %list of systems

    if ~isempty(l) && size(l,1) >= 1
        
        set(handles.lb_PaV_objectsWs,'String',l)
        
    else
        set(handles.lb_PaV_objectsWs,'Value',1);         
        set(handles.lb_PaV_objectsWs,'String', [{''}; l]);

    end
    
    
function bg_PaV_plotStyle_SelectionChangedFcn(hObject, eventdata, handles)
%Set the panel with the options for the plot style visible or invisible
%depending on the specified value

    list = get(handles.lb_PaV_selectedSystems,'String');

    if get(handles.rb_PaV_plotStyle_manual,'Value') && ~isempty(list)
        set(handles.panel_PaV_plotStyle,'Visible','on');
    else
        set(handles.panel_PaV_plotStyle,'Visible','off');
    end
    
    savePlotData(handles);

    
function pb_PaV_down_Callback(hObject, eventdata, handles)
%Change the order of the selected Systems

    list = get(handles.lb_PaV_selectedSystems,'String');
    index = get(handles.lb_PaV_selectedSystems,'Value');
    newList = list;
    
    if ~isempty(list) && index < length(list)
       
        newList{index,1} = list{index+1,1};
        newList{index+1,1} = list{index,1};
        
        set(handles.lb_PaV_selectedSystems,'String',newList);
        set(handles.lb_PaV_selectedSystems,'Value',index+1);
    end

function pb_PaV_up_Callback(hObject, eventdata, handles)
%Change the order of the selected Systems

    list = get(handles.lb_PaV_selectedSystems,'String');
    index = get(handles.lb_PaV_selectedSystems,'Value');
    newList = list;
    
    if ~isempty(list) && index > 1
       
        newList{index,1} = list{index-1,1};
        newList{index-1,1} = list{index,1};
        
        set(handles.lb_PaV_selectedSystems,'String',newList);
        set(handles.lb_PaV_selectedSystems,'Value',index-1);
    end    

function pb_PaV_showObject_Callback(hObject, eventdata, handles)
%Show the selected object in the variable editor

    %Get the selected object
    
    list = get(handles.lb_PaV_selectedSystems,'String');
    selectedObject = list{get(handles.lb_PaV_selectedSystems,'Value')};
    
    %Show object in the variable explorer
    
    try
       openvar(selectedObject); 
    catch ex 
       errordlg(strcat('The object "',selectedObject,'" does not exist in workspace!'));
       uiwait;
       return;
    end

    
function pu_PaV_color_Callback(hObject, eventdata, handles)

    savePlotData(handles);
    
function pu_PaV_lineStyle_Callback(hObject, eventdata, handles)

    savePlotData(handles);
    
function pu_PaV_plotStyle_marker_Callback(hObject, eventdata, handles)

    savePlotData(handles);

function ed_legend_Callback(hObject, eventdata, handles)
   
    savePlotData(handles);

function cb_PaV_SaveData_Callback(hObject, eventdata, handles)

    if get(hObject,'Value') == 1
       set(handles.et_PaV_saveData,'Enable','on'); 
    else
       set(handles.et_PaV_saveData,'Enable','off'); 
    end

    savePlotData(handles);

function et_PaV_saveData_Callback(hObject, eventdata, handles)

    savePlotData(handles);

function cb_PaV_SaveTimeseries_Callback(hObject, eventdata, handles)

    if get(hObject,'Value') == 1
       set(handles.et_PaV_saveTimeseries,'Enable','on'); 
    else
       set(handles.et_PaV_saveTimeseries,'Enable','off'); 
    end

    savePlotData(handles);

function et_PaV_saveTimeseries_Callback(hObject, eventdata, handles)
    
    savePlotData(handles);

function pu_in_Callback(hObject, eventdata, handles)

    handles = savePlotData(handles);
    set(handles.ed_legend,'String',suggestDefaultLegendText(handles,'',[]));
    savePlotData(handles);

function pu_out_Callback(hObject, eventdata, handles)

    handles = savePlotData(handles);
    set(handles.ed_legend,'String',suggestDefaultLegendText(handles,'',[]));
    savePlotData(handles);

    
function lb_PaV_selectedSystems_Callback(hObject, eventdata, handles)
%Selection changed. Write the stored values for the now selected system to
%the controls

    list = get(handles.lb_PaV_selectedSystems,'String');
    
    if ~isempty(list)
    
        selectedSystem = list{get(handles.lb_PaV_selectedSystems,'Value'),1};

        %Get the stored data for the system

        for i = 1:length(handles.plotData)
           if strcmp(selectedSystem,handles.plotData{i,1}.name)
              data = handles.plotData{i,1};
              break; 
           end
        end

        if ~isempty(selectedSystem)

            %Siso vs. Mimo
            
            if isempty(data.in) && isempty(data.out)    %Siso
               set(handles.panel_intoout,'Visible','off');
            else                                        %Mimo

               if get(handles.pu_in,'Value') > data.sizeInputs
                    set(handles.pu_in,'Value', 1)
               end
               if get(handles.pu_out,'Value') > data.sizeOutputs
                    set(handles.pu_out,'Value', 1)
               end
               set(handles.panel_intoout,'Visible','on')
               in={num2str((1:data.sizeInputs)'),'all'};
               out={num2str((1:data.sizeOutputs)'),'all'};
               set(handles.pu_in,'String',in)
               set(handles.pu_out,'String',out)

               set(handles.panel_intoout,'Visible','on');
               if ~strcmp(data.in,'all')
                  set(handles.pu_in,'Value',str2num(data.in));
               else
                  set(handles.pu_in,'Value',length(get(handles.pu_in,'String')));
               end
               if ~strcmp(data.out,'all')
                  set(handles.pu_out,'Value',str2num(data.out));
               else
                  set(handles.pu_out,'Value',length(get(handles.pu_out,'String')));
               end
            end
            
            %Plot style
            
            if strcmp(data.plotStyle,'manual')
               set(handles.panel_PaV_plotStyle,'Visible','on'); 
               set(handles.rb_PaV_plotStyle_auto,'Value',0);
               set(handles.rb_PaV_plotStyle_manual,'Value',1);
            else
               set(handles.panel_PaV_plotStyle,'Visible','off');
               set(handles.rb_PaV_plotStyle_auto,'Value',1);
               set(handles.rb_PaV_plotStyle_manual,'Value',0);
            end
            
            set(handles.pu_PaV_color,'Value',data.color);
            set(handles.pu_PaV_lineStyle,'Value',data.lineStyle);
            set(handles.pu_PaV_plotStyle_marker,'Value',data.markerType);
            set(handles.ed_legend,'String',data.legendText);
            
            %System vs Object
            
            if data.isSystem == 1;            
                set(handles.cb_PaV_SaveData,'Value',data.save);
                set(handles.et_PaV_saveData,'String',data.variableName);
                
                if data.save
                   set(handles.et_PaV_saveData,'Enable','on'); 
                else
                   set(handles.et_PaV_saveData,'Enable','off'); 
                end
                
                set(handles.cb_PaV_SaveTimeseries,'Value',data.saveTimeseries);
                set(handles.et_PaV_saveTimeseries,'String',data.variableNameTimeseries);
            
                if data.saveTimeseries
                   set(handles.et_PaV_saveTimeseries,'Enable','on'); 
                else
                   set(handles.et_PaV_saveTimeseries,'Enable','off'); 
                end
                
                set(handles.cb_PaV_SaveData,'Visible','on');
                set(handles.et_PaV_saveData,'Visible','on');
                set(handles.pb_PaV_infoSaveData,'Visible','on');
                
                if get(handles.plot_type,'Value') == 4 || ...
                   get(handles.plot_type,'Value') == 5     %Step or Impulse              
                    set(handles.cb_PaV_SaveTimeseries,'Visible','on');
                    set(handles.et_PaV_saveTimeseries,'Visible','on');
                    set(handles.pb_PaV_infoSaveTimeseries,'Visible','on');
                else
                    set(handles.cb_PaV_SaveTimeseries,'Visible','off');
                    set(handles.et_PaV_saveTimeseries,'Visible','off');
                    set(handles.pb_PaV_infoSaveTimeseries,'Visible','off');
                end
                
                set(handles.bg_PaV_Resolution,'Visible','on');
                set(handles.pb_PaV_showObject,'Visible','off');
                
                if strcmp(data.distribution,'linear')
                  set(handles.rb_distlin,'Value',1);
                  set(handles.rb_distlog,'Value',0);
                else
                  set(handles.rb_distlin,'Value',0); 
                  set(handles.rb_distlog,'Value',1);
                end

                set(handles.ed_min,'String',data.min);
                set(handles.ed_max,'String',data.max);
                set(handles.sl_steps,'Value',data.steps);
                set(handles.et_curstep,'String',num2str(data.steps));
                
                if strcmp(data.resolution,'manual')
                   set(handles.panel_manual,'Visible','on');                   
                   set(handles.rb_manual,'Value',1);
                   set(handles.rb_auto,'Value',0);
                   
                else
                   set(handles.panel_manual,'Visible','off');
                   set(handles.rb_manual,'Value',0);
                   set(handles.rb_auto,'Value',1);
                end
            else
                set(handles.cb_PaV_SaveData,'Visible','off');
                set(handles.et_PaV_saveData,'Visible','off');
                set(handles.pb_PaV_infoSaveData,'Visible','off');
                
                set(handles.cb_PaV_SaveTimeseries,'Visible','off');
                set(handles.et_PaV_saveTimeseries,'Visible','off');
                set(handles.pb_PaV_infoSaveTimeseries,'Visible','off');
                
                set(handles.bg_PaV_Resolution,'Visible','off');
                set(handles.panel_manual,'Visible','off');
                
                set(handles.pb_PaV_showObject,'Visible','on');
            end

        end    
    else
        
        set(handles.panel_intoout,'Visible','off');
        set(handles.panel_PaV_plotStyle,'Visible','off');
        set(handles.cb_PaV_SaveData,'Visible','off');
        set(handles.et_PaV_saveData,'Visible','off'); 
        set(handles.cb_PaV_SaveTimeseries,'Visible','off');
        set(handles.et_PaV_saveTimeseries,'Visible','off');
        set(handles.pb_PaV_infoSaveTimeseries,'Visible','off');
        set(handles.bg_PaV_plotStyle,'Visible','off');
        set(handles.bg_PaV_Resolution,'Visible','off');        
        set(handles.pb_PaV_showObject,'Visible','off');
        set(handles.pb_PaV_infoSaveData,'Visible','off');
    end
  
function lb_PaV_systemsWs_Callback(hObject, eventdata, handles)
%Update the system information

    x = get(hObject,'String');
    index = get(hObject,'Value');
    y = x{index};
    
    if ~isempty(y) && length(index) == 1
        
        %Load system from workspace

        sys = evalin('base', y);
        if ~isa(sys, 'ss') && ~isa(sys, 'sss')
            errordlg('Variable is not a valid state space model.','Error Dialog','modal')
            set(handles.sysinfo,'String','Invalid model')
            uiwait
            return
        end
        
        %Display system information

        displaySystemInformation(handles.sysinfo,sys);
    else
        set(handles.sysinfo,'String','');
    end


function plot_type_Callback(hObject, eventdata, handles)

if get(hObject,'Value')==6
    % pzmap
    set(handles.rb_auto,'Value',1)
    set(handles.rb_auto,'Enable','inactive')
    set(handles.panel_manual,'Visible','off')
    set(handles.bg_distribution,'Visible','off')
    set(handles.rb_manual,'Enable','off')
else
    set(handles.rb_auto,'Enable','on')
    set(handles.rb_manual,'Enable','on')
end

%Adapt the values for the selectable objects (frd, fir, etc.)

list = get(handles.lb_PaV_selectedSystems,'String');

switch get(handles.plot_type,'Value')
    
    case 1      %Bode
        list = removeObjectsFromList(list,'tf');
        list = removeObjectsFromList(list,'zpk');
        set(handles.panel_PaV_objectsWs,'Title','Frd-objects');
    case 2      %Magnitude
        list = removeObjectsFromList(list,'tf');
        list = removeObjectsFromList(list,'zpk');
        set(handles.panel_PaV_objectsWs,'Title','Frd-objects');
    case 3      %Singular Values
        list = removeObjectsFromList(list,'tf');
        list = removeObjectsFromList(list,'zpk');
        set(handles.panel_PaV_objectsWs,'Title','Frd-objects');
    case 4      %Impulse
        list = removeObjectsFromList(list,'frd');
        list = removeObjectsFromList(list,'zpk');
        set(handles.panel_PaV_objectsWs,'Title','Tf-objects');
    case 5      %Step
        list = removeObjectsFromList(list,'frd');
        list = removeObjectsFromList(list,'zpk');
        set(handles.panel_PaV_objectsWs,'Title','Tf-objects');
    case 6      %Pole-Zero Map
        list = removeObjectsFromList(list,'frd');
        list = removeObjectsFromList(list,'tf');    
        set(handles.panel_PaV_objectsWs,'Title','Zpk-objects');

end

l_alt = get(handles.lb_PaV_selectedSystems,'String');

if ~isempty(l_alt)

    sOld=l_alt{get(handles.lb_PaV_selectedSystems,'Value')};
    indexOld=find(strcmp(sOld,list));

    if ~isempty(indexOld)
        set(handles.lb_PaV_selectedSystems,'Value',indexOld);
    else
        set(handles.lb_PaV_selectedSystems,'Value',1);
    end

else
  set(handles.lb_PaV_selectedSystems,'Value',1);
end

if ~isempty(list) && size(list,1) >= 1

    set(handles.lb_PaV_selectedSystems,'String',list)

else
    set(handles.lb_PaV_selectedSystems,'Value',1);         
    set(handles.lb_PaV_selectedSystems,'String', {});
end

%Set option to save timeseries data visible if "step" is selected

if get(handles.plot_type,'Value') == 4 || ...
   get(handles.plot_type,'Value') == 5                    %Step or Impulse
    set(handles.cb_PaV_SaveTimeseries,'Visible','on');
    set(handles.et_PaV_saveTimeseries,'Visible','on');
    set(handles.pb_PaV_infoSaveTimeseries,'Visible','on');
else
    set(handles.cb_PaV_SaveTimeseries,'Visible','off');
    set(handles.et_PaV_saveTimeseries,'Visible','off');
    set(handles.pb_PaV_infoSaveTimeseries,'Visible','off');
end


%Change the variable name for the stored plot-data

for i = 1:length(handles.plotData)
   handles.plotData{i,1}.variableName = suggestPlotDataName(handles, ...
       handles.plotData{i,1}.name);
   handles.plotData{i,1}.variableNameTimeseries = suggestPlotDataNameTimeseries(handles, ...
       handles.plotData{i,1}.name);
end

guidata(hObject,handles);
pb_PaV_refeshObjects_Callback(handles.pb_PaV_refeshObjects, eventdata, handles)
lb_PaV_selectedSystems_Callback(handles.lb_PaV_selectedSystems,eventdata,handles);

    
function rb_auto_Callback(hObject, eventdata, handles)
% automatic choice of time/frequency
set(handles.panel_manual,'Visible','off')
set(handles.bg_distribution,'Visible','off')

savePlotData(handles);

function rb_manual_Callback(hObject, eventdata, handles)
% manual choice of time/frequency
set(handles.panel_manual,'Visible','on')
set(handles.bg_distribution,'Visible','on')

savePlotData(handles);

function ed_min_Callback(hObject, eventdata, handles)
% check user input
testMin(hObject,handles) % if ok, save 0 to userdata, otherwise 1
if get(handles.rb_manual,'Value')==0 % manual selected but inactive
    set(handles.rb_manual,'Value',1)
    testMax(handles.ed_max,handles)
end

savePlotData(handles);

function ed_max_Callback(hObject, eventdata, handles)
% see ed_min_Callback
testMax(hObject,handles)
if get(handles.rb_manual,'Value')==0
    set(handles.rb_manual,'Value',1)
    testMin(handles.ed_min,handles)
end

savePlotData(handles);

function sl_steps_Callback(hObject, eventdata, handles)
% round to 100
h=get(hObject,'Value');
k=mod(h,100);
if k>50
    h=(100-k)+h;
else
    h=h-k;
    if h<10
        h=10;
    end
end
set(hObject,'Value',h)
% show value in textfield
set(handles.et_curstep,'String',h)
% check min and max
if get(handles.rb_manual,'Value')==0
    set(handles.rb_manual,'Value',1)
    testMax(handles.ed_max,handles)
    testMin(handles.ed_min,handles)
end

savePlotData(handles);

function et_curstep_Callback(hObject, eventdata, handles)
%Check if the value lies in the allowed interval and update the slider

    testWidth(hObject);
    
    dataIncorrect = get(hObject,'UserData');
    
    value = str2num(get(hObject,'String'));
    
    if dataIncorrect == 0
       if value < 10
           errordlg('Value has to be greater or equal than 10','Error Dialog','modal')
           uiwait
           set(hObject,'String',num2str(get(handles.sl_steps,'Value')));
           return
       elseif value > 10000
           errordlg('Value has to be smaller or equal than 10000','Error Dialog','modal')
           uiwait
           set(hObject,'String',num2str(get(handles.sl_steps,'Value')));
           return
       end
       
       set(handles.sl_steps,'Value',str2num(get(hObject,'String'))); 
    else
       set(hObject,'String',num2str(get(handles.sl_steps,'Value')));
    end
    
    savePlotData(handles);

function bg_distribution_SelectionChangedFcn(hObject, eventdata, handles)

savePlotData(handles);


function pb_PaV_infoSaveData_Callback(hObject, eventdata, handles)
%Show a information-box with information about the selectable options

    infoBox({'InfoSaveData.png'});
    uiwait;

function pb_PaV_infoSaveTimeseries_Callback(hObject, eventdata, handles)
%Show a information-box with information about the selectable options

    infoBox({'InfoSaveTimeseries.png'});
    uiwait;
    
    
function pb_plot_Callback(hObject, eventdata, handles)
%Plot the graph with the choosen systems

    set(hObject,'String','Busy');
    set(hObject,'Enable','off'); 
    set(handles.figure1,'Pointer','watch')
    drawnow();
    
    try
    
        %Read out all choosen systems from workspace

        list = get(handles.lb_PaV_selectedSystems,'String');

        if isempty(list)
           error('No system selected!'); 
        end
        
        systemList = {};
    
        for i = 1:length(list)
            
            sysname = list{i,1};
            
            try
                sys = evalin('base', sysname);            
            catch ex
               pb_PaV_refeshObjects_Callback(handles.pb_PaV_refeshObjects, eventdata, handles);
               pb_refreshsys_Callback(handles.pb_refreshsys, eventdata, handles)
               error(strcat('The system "',sysname,'" could not be loaded from workspace')); 
            end
            
            systemList{i,1} = sysname;
            systemList{i,2} = sys;
            
            for j = 1:length(handles.plotData)
               if strcmp(handles.plotData{j,1}.name,sysname)
                  systemList{i,3} = handles.plotData{j,1};
                  break; 
               end
            end
        end
        
        %Chosen input- and output- channels for MIMO-systems
        
        for i = 1:size(systemList,1);
           if (isa(systemList{i,2},'sss') || isa(systemList{i,2},'ss')) && ...
               (size(systemList{i,2}.B,2) > 1 || size(systemList{i,2}.C,1) > 1) %Mimo
              in = systemList{i,3}.in;
              out = systemList{i,3}.out; 
              if ~strcmp(in,'all') && ~strcmp(out,'all')
                 systemList{i,2} = systemList{i,2}(str2num(in),str2num(out)); 
              elseif ~strcmp(in,'all')
                 systemList{i,2} = systemList{i,2}(1:size(systemList{i,2}.B,2),str2num(out));
              elseif ~strcmp(out,'all')
                 systemList{i,2} = systemList{i,2}(1:size(systemList{i,2}.C,1),str2num(out));
              end              
           end
        end
        
        %Frequency/Time-vector
        
        for i = 1:size(systemList,1)
            if (isa(systemList{i,2},'sss') || isa(systemList{i,2},'ss')) && ... 
                    strcmp(systemList{i,3}.resolution,'manual')

                data = systemList{i,3};
                
                if isempty(data.min)
                    error('Please correct minimal frequency/time first');
                elseif isempty(data.max)
                    error('Please correct maximal frequency/time first');
                end
                minimum=str2num(data.min);
                maximum=str2num(data.max);
                steps=data.steps;
                if strcmp(data.distribution,'linear')
                    frequency=linspace(minimum,maximum,steps);
                else
                    if minimum <= 0
                        error('Minimum frequency/time must be greater than zero if you choose logarithmical distribution!')
                    end
                    frequency=logspace(log10(minimum),log10(maximum),steps);
                end
                
                systemList{i,6} = frequency;
            else
                systemList{i,6} = [];
            end
        end
        
        %Plot-style
           
        for i = 1:size(systemList,1)

            if strcmp(systemList{i,3}.plotStyle,'manual')
            
                switch systemList{i,3}.lineStyle;

                    case 1
                        lineStyle = '-';
                    case 2
                        lineStyle = '--';
                    case 3
                        lineStyle = ':';
                    case 4
                        lineStyle = '.-';                    
                end

                switch systemList{i,3}.color

                    case 1
                        color = 'b';
                    case 2
                        color = 'r';
                    case 3
                        color = 'g';
                    case 4
                        color = 'm';
                    case 5
                        color = 'c';
                    case 6
                        color = 'k';
                    case 7
                        color = 'y';                   
                end
                
                switch systemList{i,3}.markerType
                  
                    case 1
                        markerType = '';
                    case 2
                        markerType = '+';
                    case 3
                        markerType = 'o';
                    case 4
                        markerType = '*';
                    case 5
                        markerType = '.';
                    case 6
                        markerType = 'x';
                    case 7
                        markerType = 's';
                    case 8
                        markerType = 'd';
                    case 9
                        markerType = '^';
                    case 10
                        markerType = 'v';
                    case 11
                        markerType = '>';
                    case 12
                        markerType = '<';
                    case 13
                        markerType = 'p';
                    case 14
                        markerType = 'h';
                end
                
                if strcmp(markerType,'')
                   systemList{i,4} = strcat(lineStyle,color);
                else
                   systemList{i,4} = strcat(lineStyle,markerType,color);
                end
            else
                systemList{i,4} = [];
            end
        end       
        
        % Plot style manual or auto
        
        plotStyleManual = 0;
        for i = 1:size(systemList,1)
            if strcmp(systemList{i,3}.plotStyle,'manual')
                plotStyleManual = 1;
                break;
            end
        end
        
        %Graph-Type
        
        figure;
        figureHandles = gcf;
        set(figureHandles,'NumberTitle','off');
        
        switch get(handles.plot_type,'Value')
                
            case 1      %Bode
                
                set(figureHandles,'Name','Bode Diagram');
                
                % convert to frd-objects
                for i = 1:size(systemList,1)
                   if isa(systemList{i,2},'sss')
                       if strcmp(systemList{i,3}.resolution,'manual') 
                            systemList{i,5} = freqresp(systemList{i,2},systemList{i,6},struct('frd',1));
                       else
                            systemList{i,5} = freqresp(systemList{i,2},struct('frd',1));
                       end
                   elseif isa(systemList{i,2},'ss')
                       if strcmp(systemList{i,3}.resolution,'manual') 
                           frequencyResponse = freqresp(systemList{i,2},systemList{i,6});
                           systemList{i,5} = frd(frequencyResponse,systemList{i,6},systemList{i,2});
                       else
                           [frequencyResponse,omega] = freqresp(systemList{i,2});
                           systemList{i,5} = frd(frequencyResponse,omega,systemList{i,2});
                       end
                   else
                       systemList{i,5} = systemList{i,2}; 
                   end
                end
                
                % join arguments together
                if plotStyleManual == 1                    
                    arguments = cell(2*size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{2*(i-1)+1}=systemList{i,5};
                        arguments{2*(i-1)+2}=systemList{i,4};
                    end
                    arguments(cellfun(@isempty,arguments)) = [];                 
                else                    
                    arguments = cell(size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{i}=systemList{i,5};
                    end                   
                end
                
                % create bode-plot
                bode(arguments{:});
            
            case 2      %Magnitude
                
                set(figureHandles,'Name','Bode Diagram');
                
                % convert to frd-objects
                for i = 1:size(systemList,1)
                   if isa(systemList{i,2},'sss')
                       if strcmp(systemList{i,3}.resolution,'manual') 
                            systemList{i,5} = freqresp(systemList{i,2},systemList{i,6},struct('frd',1));
                       else
                            systemList{i,5} = freqresp(systemList{i,2},struct('frd',1));
                       end
                   elseif isa(systemList{i,2},'ss')
                       if strcmp(systemList{i,3}.resolution,'manual') 
                           frequencyResponse = freqresp(systemList{i,2},systemList{i,6});
                           systemList{i,5} = frd(frequencyResponse,systemList{i,6},systemList{i,2});
                       else
                           [frequencyResponse,omega] = freqresp(systemList{i,2});
                           systemList{i,5} = frd(frequencyResponse,omega,systemList{i,2});
                       end
                   else
                       systemList{i,5} = systemList{i,2}; 
                   end
                end
                
                % join arguments together
                if plotStyleManual == 1                    
                    arguments = cell(2*size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{2*(i-1)+1}=systemList{i,5};
                        arguments{2*(i-1)+2}=systemList{i,4};
                    end
                    arguments(cellfun(@isempty,arguments)) = [];                 
                else                    
                    arguments = cell(size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{i}=systemList{i,5};
                    end                   
                end
                
                % create bodemag-plot
                bodemag(arguments{:});
                
            case 3      %Singular values
                
                set(figureHandles,'Name','Singular Values');
                
                % convert to frd-objects
                for i = 1:size(systemList,1)
                   if isa(systemList{i,2},'sss')
                       if strcmp(systemList{i,3}.resolution,'manual') 
                            systemList{i,5} = freqresp(systemList{i,2},systemList{i,6},struct('frd',1));
                       else
                            systemList{i,5} = freqresp(systemList{i,2},struct('frd',1));
                       end
                   elseif isa(systemList{i,2},'ss')
                       if strcmp(systemList{i,3}.resolution,'manual') 
                           frequencyResponse = freqresp(systemList{i,2},systemList{i,6});
                           systemList{i,5} = frd(frequencyResponse,systemList{i,6},systemList{i,2});
                       else
                           [frequencyResponse,omega] = freqresp(systemList{i,2});
                           systemList{i,5} = frd(frequencyResponse,omega,systemList{i,2});
                       end
                   else                 %Frd-object
                       systemList{i,5} = systemList{i,2}; 
                   end
                end
                
                % join arguments together
                if plotStyleManual == 1                    
                    arguments = cell(2*size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{2*(i-1)+1}=systemList{i,5};
                        arguments{2*(i-1)+2}=systemList{i,4};
                    end
                    arguments(cellfun(@isempty,arguments)) = [];                 
                else                    
                    arguments = cell(size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{i}=systemList{i,5};
                    end                   
                end
                
                % create sigma-plot
                sigma(arguments{:});
                                
            case 4      %Impulse Response
                
                set(figureHandles,'Name','Impulse Response');
                
                % convert to tf-objects
                
                Tmax = 0;       %Plot from 0 to Tmax 
                Ttemp = 0;
                
                for i = 1:size(systemList,1)
                    if isa(systemList{i,2},'sss') || isa(systemList{i,2},'ssRed')
                       if strcmp(systemList{i,3}.resolution,'manual') 
                            [systemList{i,5},systemList{i,7}.h,systemList{i,7}.t] ...
                                = impulse(systemList{i,2},systemList{i,6},struct('tf',1));
                            Ttemp = systemList{i,6}(end);
                       else
                            [systemList{i,5},systemList{i,7}.h,systemList{i,7}.t] ...
                                = impulse(systemList{i,2},struct('tf',1));
                            Ttemp = max(cellfun(@length,systemList{i,5}.num(:)))*systemList{i,5}.Ts;
                       end
                   else
                       systemList{i,5} = systemList{i,2};
                       Ttemp = max(cellfun(@length,systemList{i,5}.num(:)))*systemList{i,5}.Ts;
                    end
                    if Ttemp > Tmax
                        Tmax = Ttemp;
                    end
                end
                
                % join arguments together
                if plotStyleManual == 1                    
                    arguments = cell(2*size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{2*(i-1)+1}=systemList{i,5};
                        arguments{2*(i-1)+2}=systemList{i,4};
                    end
                    arguments(cellfun(@isempty,arguments)) = [];                 
                else                    
                    arguments = cell(size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{i}=systemList{i,5};
                    end                   
                end
                
                % create impulse-plot
                impulse(arguments{:});
                
            case 5      %Step Response
                
                set(figureHandles,'Name','Step Response');
                
                % convert to tf-objects
                
                Tmax = 0;       %Plot from 0 to Tmax 
                Ttemp = 0;
                
                for i = 1:size(systemList,1)
                    if isa(systemList{i,2},'sss') || isa(systemList{i,2},'ssRed')
                       if strcmp(systemList{i,3}.resolution,'manual') 
                            [systemList{i,5},systemList{i,7}.h,systemList{i,7}.t] ...
                                = step(systemList{i,2},systemList{i,6},struct('tf',1));
                            Ttemp = systemList{i,6}(end);
                       else
                            [systemList{i,5},systemList{i,7}.h,systemList{i,7}.t] ...
                                = step(systemList{i,2},struct('tf',1));
                            Ttemp = max(cellfun(@length,systemList{i,5}.num(:)))*systemList{i,5}.Ts;
                       end
                   else
                       systemList{i,5} = systemList{i,2}; 
                       Ttemp = max(cellfun(@length,systemList{i,5}.num(:)))*systemList{i,5}.Ts;
                   end
                   if Ttemp > Tmax
                      Tmax = Ttemp; 
                   end
                end
                
                % join arguments together
                if plotStyleManual == 1                    
                    arguments = cell(2*size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{2*(i-1)+1}=systemList{i,5};
                        arguments{2*(i-1)+2}=systemList{i,4};
                    end
                    arguments(cellfun(@isempty,arguments)) = [];                 
                else                    
                    arguments = cell(size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{i}=systemList{i,5};
                    end                   
                end
                
                % create step-plot
                step(arguments{:});
                
            case 6      %Pole-Zero-Map
                
                set(figureHandles,'Name','Pole-Zero Map');
                
                % convert to zpk-objects
                for i = 1:size(systemList,1)
                   if isa(systemList{i,2},'sss')
                       systemList{i,5} = zpk(systemList{i,2},struct('zpk',1));
                   elseif isa(systemList{i,2},'ssRed')
                       systemList{i,5} = zpk(systemList{i,2});
                   else
                       systemList{i,5} = systemList{i,2}; 
                   end
                end
                
                % join arguments together
                if plotStyleManual == 1                    
                    arguments = cell(2*size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{2*(i-1)+1}=systemList{i,5};
                        arguments{2*(i-1)+2}=systemList{i,4};
                    end
                    arguments(cellfun(@isempty,arguments)) = [];                 
                else                    
                    arguments = cell(size(systemList,1),1);
                    for i = 1:size(systemList,1)
                        arguments{i}=systemList{i,5};
                    end                   
                end
                
                % create pzmap-plot
                pzmap(arguments{:});
                                    
        end
        
        %Legend
        
        l = legend('show');
        legendText = get(l,'String');
        
        for i = 1:size(systemList,1)
           if strcmp(systemList{i,3}.plotStyle,'manual')
               legendText{i} = systemList{i,3}.legendText;
           end            
        end
        
        set(l,'String',legendText);
        
        %Save Data
        
        for i = 1:size(systemList,1)
            if systemList{i,3}.save
                assignin('base',systemList{i,3}.variableName,systemList{i,5})
            end
        end

        %Save Timeseries
        
        if get(handles.plot_type,'Value') == 4 || ...
           get(handles.plot_type,'Value') == 5           %Step or Impulse
           for i = 1:size(systemList,1)
               if systemList{i,3}.saveTimeseries
                    assignin('base',systemList{i,3}.variableNameTimeseries, ...
                             systemList{i,7});
            
               end
           end
        end
        
        set(hObject,'String','Plot')
        set(hObject,'Enable','on')
        set(handles.figure1,'Pointer','arrow')
        
        %Update the list of objects (frd, fir, etc.) in workspace, 
        %because there may be new objects after plotting (if "Save data" is
        %selected
        
        pb_PaV_refeshObjects_Callback(handles.pb_PaV_refeshObjects, eventdata, handles)
        
    catch ex
       % check which function produced the error and display the line of
       % the error in the error-dialog that is popping up
       functionList = {'bode','impulse','step','freqresp','pzmap','bodemag'};
       for i = 1:length(functionList)
            [found,str] = extractEntryFromErrorStack(ex.stack,functionList{i});
            if found            
                errordlg({str;ex.message},'Error Dialog','modal');
                break;
            end
       end
       % show the normal error dialog if the error was not produced from a
       % toolbox function
       if ~found
            errordlg(ex.message);
       end
       uiwait
       set(hObject,'String','Plot')
       set(hObject,'Enable','on')
       set(handles.figure1,'Pointer','arrow')
       return
    end
    
        

%--------------------------------------------------------------------------
%                     LOADING AND SETTING UP MODELS
%--------------------------------------------------------------------------

%Callbacks of the buttons

function pb_load_Callback(hObject, eventdata, handles)

    % disable to avoid double call of uigetfile
    set(handles.allbuttons,'Enable','off')
    filename=sprintf('%s.mat',handles.letzterpfad);
    
    loadingSuccessfull = 1;
    
    %Check which option is selected
    
    if get(handles.rb_loadOptions_matrices,'Value') == 1  % Matrices
        
        %Open dialog-box
        [filename,path]=uigetfile(filename);
        set(handles.allbuttons,'Enable','on')
        
        %Check if files are correct
        if filename==0
            return
        end
        handles.letzterpfad=path;
        guidata(hObject, handles);
        if isempty(strfind(filename,'.mat'))  %only mat-files can be loaded
            errordlg('Only .mat files allowed.','Error Dialog','modal')
            uiwait
            return
        end 
       
        %Load matrices to workspace
        
        evalin('base',sprintf('load(''%s%s'');',path,filename));
        
    else                                                % System
        
        %Open dialog-box
        [fileList,path]=uigetfile(filename,'MultiSelect','on');
        set(handles.allbuttons,'Enable','on')
        
        handles.letzterpfad=path;
        guidata(hObject, handles);        
        
        %Create system with the function loadSss
        
        if iscell(fileList)             % Multi-Select
            loadingSuccessfull = 0;
            for i = 1:length(fileList)
                filename = fileList{i};
                if isempty(strfind(filename,'.mat'))  %only mat-files can be loaded
                    errordlg(strcat(filename,': Only .mat files allowed.'), ...
                                    'Error Dialog','modal')
                    uiwait
                else
                    success = loadSystemWithLoadSss(filename,path); 
                    if success
                        loadingSuccessfull = 1;
                    end
                end
            end
        else                            % single System
            loadingSuccessfull = loadSystemWithLoadSss(fileList,path);  
        end
        
        
    end
    
    %Refresh the display of the variables in workspace
    
    set(handles.lb_systems,'Value',[])
    set(handles.lb_systems,'String',systemsInWorkspace)
    set(handles.lb_matrixes,'Value',[])
    set(handles.lb_matrixes,'String',matricesInWorkspace)
    pb_refreshsys_Callback(hObject, eventdata, handles);
    
    if loadingSuccessfull == 1
        
        if get(handles.rb_loadOptions_matrices,'Value') == 1
            
            msgbox('Loading of the matrices was successful','Information','modal')
            uiwait            
        else            
            msgbox('Loading of the system was successful','Information','modal')
            uiwait           
        end
    end
    
function pb_readascii_Callback(hObject, eventdata, handles)
    % import matrix market format

    % disable button to avoid doubleclick
    set(handles.allbuttons,'Enable','off')
    path=sprintf('%s*.*',handles.letzterpfad);
    [filename,path]=uigetfile(path);
    set(handles.allbuttons,'Enable','on')
    if filename==0
        return
    end
    handles.letzterpfad=path;
    guidata(hObject, handles);
    k=strfind(filename,'.');
    % suggest filenames (without extension)
    if length(filename)>k
        suggestion=filename(1:(k-1));
        suggestion=suggestVarname(suggestion,[]);
    else
        suggestion=[];
    end
    matname=inputdlg({'Please enter a name for the selected matrix:'},'Input',...
                    1,{suggestion});
    if isempty(matname)
        errordlg('The choosen name is not a valid variable name.','Error Dialog','modal')
        uiwait
        return        
    else
        if ~isvarname(matname{1,1})
            errordlg('The choosen name is not a valid variable name.','Error Dialog','modal')
            uiwait
            return
        elseif existInBaseWs(matname{1,1})==1
            errordlg('There already exists a varible with this name in the workspace.','Error Dialog','modal')
            uiwait
            return
        end
    end
    % import matrix
    try
        A=readMatrixMarket(strcat(path,filename));
    catch ex
        errordlg(ex.message,'Error Dialog','modal') 
        uiwait
        return
    end
    if isempty(A)
        errordlg('Could not find matrix dimensions.','Error Dialog','modal')
        uiwait
        return
    end
    assignin('base',matname{1,1},A)
    msgbox('Matrix was loaded to workspace.','Information','modal')
    uiwait
    pb_refreshlb_Callback(hObject, eventdata, handles)

function pb_create_Callback(hObject, eventdata, handles)
    % disable to avoid double click
    set(handles.allbuttons,'Enable','off')
    composeModel('load',handles)
    set(handles.allbuttons,'Enable','on')

function pb_infoLoadOptions_Callback(hObject, eventdata, handles)

    %Show a information-box with information about the selectable options

    infoBox({'InfoLoadOptions.png'});
    uiwait;

    
%Callbacks of the list-boxes
    
function lb_systems_Callback(hObject, eventdata, handles)
    set(hObject,'String',systemsInWorkspace)
    set(handles.lb_matrixes,'Value',[]);

function lb_systems_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String',systemsInWorkspace)

function lb_matrixes_Callback(hObject, eventdata, handles)
    set(hObject,'String',matricesInWorkspace)
    set(handles.lb_systems,'Value',[]);

function lb_matrixes_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String',matricesInWorkspace)
    
 
%Callbacks of the buttons in "Save matrices and systems in workspace"

function pb_save_Callback(hObject, eventdata, handles)
% save variables to file

% check whether system or matrix is selected
if isempty(get(handles.lb_systems,'Value')) && isempty(get(handles.lb_matrixes,'Value'))
    errordlg('Please mark variables you want to save in listbox','Error Dialog','modal')
    uiwait
else
    % disable to avoid doubleclick
    set(handles.allbuttons,'Enable','off')
    % (uiputfile sometimes takes a long time to execute)
    [filename,path,filter] = uiputfile('*.mat');
    set(handles.allbuttons,'Enable','on')
    if filename==0
        return
    end
    if filter~=1
        errordlg('Only *.mat files allowed','Error Dialog','modal')
        uiwait
        return
    end
    string=get(handles.lb_systems,'String'); % all in listbox
    value=get(handles.lb_systems,'Value');% all that are selected
    variables=cell(length(value),1);
    for i=1:length(value)
        % write all selected ones to string
        variables{i}=string{value(i)};
    end
    string=get(handles.lb_matrixes,'String');
    value=get(handles.lb_matrixes,'Value');
    for i=1:length(value)
        % add selected matrices
        variables=[variables, string(value(i))]; %#ok<AGROW>
    end
    
    try
        filename=sprintf('save(''%s%s'',temp{:});',path,filename);
        % save variable list in base workspace
        assignin('base','temp', variables)    %*** ??
        % store selected variables
        evalin('base',filename)
        % delete temp variables from base workspace
        evalin('base', 'clear temp')  %*** ??
        msgbox('Saving was successful.','Information','modal')
        uiwait
    catch ex
       errordlg(strcat('Error occured during saving: ',ex.message))
       uiwait
    end
end

function pb_delete_Callback(hObject, eventdata, handles)
% delete variables

string=get(handles.lb_systems,'String'); % all variables in listbox
value=get(handles.lb_systems,'Value'); % all selected variables
variables=cell(length(value),1);
for i=1:length(value)
    % write all selected ones to string
    variables{i}=string{value(i)};
end
string=get(handles.lb_matrixes,'String');
value=get(handles.lb_matrixes,'Value');
for i=1:length(value)
    % add selected matrices
    variables=[variables, string(value(i))]; %#ok<AGROW>
end
k=stqd('String','Are you sure you want to delete all selected systems and matrices from your workspace?','Title','Question Dialog');
if ~isempty(k) && strcmp(k,'No')
    return
end

assignin('base','temp',variables)
evalin('base','clear(temp{:},''temp'')')
pb_refreshsys_Callback(0, 0, handles)
pb_refreshlb_Callback(0, 0, handles)

function pb_refreshlb_Callback(hObject, eventdata, handles)
set(handles.lb_systems,'Value',[])
set(handles.lb_systems,'String',systemsInWorkspace)
set(handles.lb_matrixes,'Value',[])
set(handles.lb_matrixes,'String',matricesInWorkspace)


%Callbacks of the elements in "Extrakt matrices from system"

function ed_extract_A_Callback(hObject, eventdata, handles)
% check variable names
set(hObject,'UserData',0)
if ~isempty(get(hObject,'String'))
    isvalidvarname(hObject)
end
    
function ed_extract_B_Callback(hObject, eventdata, handles)
% check variable names
set(hObject,'UserData',0)
if ~isempty(get(hObject,'String'))
    isvalidvarname(hObject)
end

function ed_extract_C_Callback(hObject, eventdata, handles)
% check variable names
set(hObject,'UserData',0)
if ~isempty(get(hObject,'String'))
    isvalidvarname(hObject)
end

function ed_extract_D_Callback(hObject, eventdata, handles)
% check variable names
set(hObject,'UserData',0)
if ~isempty(get(hObject,'String'))
    isvalidvarname(hObject)
end

function ed_extract_E_Callback(hObject, eventdata, handles)
% check variable names
set(hObject,'UserData',0)
if ~isempty(get(hObject,'String'))
    isvalidvarname(hObject)
end

function pb_load_sug_Callback(hObject, eventdata, handles)
% suggest variable names
x = get(handles.lb_systems,'String');
try
    % return only first selected system
    y=x{get(handles.lb_systems,'Value')};
catch ex
    % no system has been selected
    if strfind(ex.identifier,'MATLAB:TooManyOutputsDueToMissingBraces')
        errordlg('Choose ONE system in the listbox above','Error Dialog','modal')
        uiwait
        return
    else
        throw(ex)
    end
end
% suggest variable names that do not exist in workspace
suggestVarname(sprintf('%s_A',y),handles.ed_extract_A);
suggestVarname(sprintf('%s_B',y),handles.ed_extract_B);
suggestVarname(sprintf('%s_C',y),handles.ed_extract_C);
suggestVarname(sprintf('%s_D',y),handles.ed_extract_D);
suggestVarname(sprintf('%s_E',y),handles.ed_extract_E);

function pb_load_extract_Callback(hObject, eventdata, handles)
x = get(handles.lb_systems,'String');
try
    % return only first selected system
    y=x{get(handles.lb_systems,'Value')};
catch ex
    if strfind(ex.identifier,'MATLAB:TooManyOutputsDueToMissingBraces')
        errordlg('Choose ONE system in the listbox above.','Error Dialog','modal')
        uiwait
        return
    else
        throw(ex)
    end
end
if isempty(y) || length(get(handles.lb_systems,'Value'))>1 
    % more than one system has been selected
    errordlg('Choose ONE system in the listbox above','Error Dialog','modal')
    uiwait
    return
end
% UserData contains 1, if invalid variable name was entered
if get(handles.ed_extract_A,'UserData')==1||...
        get(handles.ed_extract_B,'UserData')==1||...
        get(handles.ed_extract_C,'UserData')==1 ||...
        get(handles.ed_extract_D,'UserData')==1 ||...
        get(handles.ed_extract_E,'UserData')==1
    errordlg('Correct names for matrices first!','Error Dialog','modal')
    uiwait
    return
end
%Check if the a name for all systems is selected
if isempty(get(handles.ed_extract_A,'String'))||...
   isempty(get(handles.ed_extract_B,'UserData'))||...
   isempty(get(handles.ed_extract_C,'UserData')) ||...
   isempty(get(handles.ed_extract_D,'UserData')) ||...
   isempty(get(handles.ed_extract_E,'UserData'))
    errordlg('Correct names for matrices first!','Error Dialog','modal')
    uiwait
    return
end

if get(handles.ed_extract_A,'UserData')==1||...
        get(handles.ed_extract_B,'UserData')==1||...
        get(handles.ed_extract_C,'UserData')==1 ||...
        get(handles.ed_extract_D,'UserData')==1 ||...
        get(handles.ed_extract_E,'UserData')==1
    errordlg('Correct names for matrices first!','Error Dialog','modal')
    uiwait
    return
end

sys=evalin('base',y);

z=[handles.ed_extract_A,handles.ed_extract_B,handles.ed_extract_C,handles.ed_extract_D,handles.ed_extract_E];
w='A';
% if a field is left empty, matrix is not loaded to workspace
for i=z
    if isvarname(get(i,'String'))
        assignin('base',get(i,'String'),evalin('base',sprintf('%s.%s',y,w)))
    end
    w=char(w+1);
end
% refresh systems and matrices lb
pb_refreshlb_Callback(hObject, eventdata, handles)
msgbox('Extraction was successful.','Information','modal')
uiwait


%Callbacks for "Truncate Model"

function pb_LaSuM_truncate_Callback(hObject, eventdata, handles)
%Open the Sub-GUI for Model-Truncation

    list = get(handles.pu_LaSuM_truncate,'String');
    sysName = list{get(handles.pu_LaSuM_truncate,'Value'),1};

    try
        sys = evalin('base',sysName);
    catch
        errordlg('No valid system selected! Please correct this first!','Error Dialog','modal');
        return;
    end
    
    truncateModel({sys});
    uiwait;
    
    % Refresh the displayes of systems in the workspace
    
    pb_refreshsys_Callback(hObject, eventdata, handles);
    set(handles.lb_systems,'Value',[])
    set(handles.lb_systems,'String',systemsInWorkspace)


%--------------------------------------------------------------------------
%                         MODEL ORDER REDUCTION
%--------------------------------------------------------------------------

%Callbacks of the controls on the top

function pu_mor_systems_Callback(hObject, eventdata, handles)
    % original system has been selected

    x = get(hObject,'String');
    y = x{get(hObject,'Value')};
    if isempty(y)
        set(handles.pb_mor_reduce,'Enable','off')
        set(handles.panel_mor_hsv,'Visible','off')
        set(handles.st_mor_sysinfo,'String','Please choose a system.')
        set(handles.ed_mor_q,'Enable','off')
        set(handles.sl_mor_q,'Enable','off')
        return
    end

    try
        sys = evalin('base', y);
    catch ex %#ok<NASGU>
        set(handles.pb_mor_reduce,'Enable','off')
        set(handles.panel_mor_hsv,'Visible','off')
        set(handles.ed_mor_q,'Enable','off')
        set(handles.sl_mor_q,'Enable','off')
        set(handles.st_mor_sysinfo,'String','Invalid model')
        errordlg('Variable is not a valid state space model.','Error Dialog','modal')
        uiwait
        return
    end

    % set max of slider to system order
    if get(handles.sl_mor_q,'Value') > sys.n || ... 
       str2double(get(handles.ed_mor_q,'String')) > sys.n
   
       warndlg('Order of reduced system must not be greater than order of original system.','Warning','modal')
       uiwait
       
       set(handles.sl_mor_q,'Value',1)
       set(handles.ed_mor_q,'String','1')
    end
    
    displaySystemInformation(handles.st_mor_sysinfo,sys);
    
    if sys.n == 1
        set(handles.sl_mor_q,'Enable','off')
    else
        set(handles.sl_mor_q,'Max',size(sys.A,1));
        if get(handles.pu_mor_method,'Value') ~= 3  %not Krylov         
            set(handles.sl_mor_q,'Enable','on')
        end
    end
    set(handles.pb_mor_reduce,'Enable','on');
    
    %If Krylov is selected, set the default values for the expension points and
    %show the correct panel dependent on if the system is siso or mimo

    if get(handles.pu_mor_method,'Value')==3


         if sys.m > 1 || sys.p > 1    %MIMO

           set(handles.panel_mor_krylov_MimoExps,'Visible','on');
           %set(handles.panel_mor_krylov_exps,'Visible','off');

           layoutMimoKrylov(handles);

        else                         %SISO

           set(handles.panel_mor_krylov_MimoExps,'Visible','off');
           %set(handles.panel_mor_krylov_exps,'Visible','on');

           layoutSisoKrylov(handles);
        end

        setKrylovDefaultValues(handles);
    end
    
    %If TBR is selected, check wether the HSV for this system have been
    %calculated already
    
    if get(handles.pu_mor_method,'Value')==1
        
        x = get(handles.pu_mor_systems,'String');
        y = x{get(handles.pu_mor_systems,'Value')};

        if ~isempty(y)
            try
               sys = evalin('base', y);

               if getHankelSingularValues(sys,y,handles)
                   updateTBR(handles.pb_mor_hsv,eventdata,handles);
               elseif isa(sys,'sss') && ~isempty(sys.HankelSingularValues) && ...
                      ~isempty(sys.TBal) && ~isempty(sys.TBalInv)
                   handles = saveHankelSingularValues(sys,y,sys.HankelSingularValues,sys.TBal,sys.TBalInv,handles);
                   updateTBR(handles.pb_mor_hsv,eventdata,handles);
               else                  
                   set(handles.panel_mor_hsv,'Visible','off');
               end 
            catch ex
                set(handles.panel_mor_hsv,'Visible','off');
            end
        else
           set(handles.panel_mor_hsv,'Visible','off');
        end        
    end
    
    % suggest names for reduced system and projection matrices
    
    suggestNamesMOR(y,handles);    

function pb_refreshsys_Callback(hObject, eventdata, handles)

    %Refresh list of open figures and systems
    
    %Get all systems from workspace
    
    l=systemsInWorkspace();
    
    %Check wether the previous selected systems exist in workspace
    
    
    %Postprocessing and Visualisation
    
    l_alt = get(handles.lb_PaV_systemsWs,'String');

    if ~isempty(l_alt)

        sOld=l_alt{get(handles.lb_PaV_systemsWs,'Value')};
        indexOld=find(strcmp(sOld,l));

        if ~isempty(indexOld)
            set(handles.lb_PaV_systemsWs,'Value',indexOld);
        else
            set(handles.lb_PaV_systemsWs,'Value',1);
        end

    else
      set(handles.lb_PaV_systemsWs,'Value',1);
    end
    
    %System Analysis

    l_alt = get(handles.pu_an_sys1,'String');

    if ~isempty(l_alt)

        sOld=l_alt{get(handles.pu_an_sys1,'Value')};
        indexOld=find(strcmp(sOld,l));

        if ~isempty(indexOld)
            set(handles.pu_an_sys1,'Value',indexOld);
        else
            set(handles.pu_an_sys1,'Value',1);
        end

    else
      set(handles.pu_an_sys1,'Value',1);
    end 
    
    
    %Model Order Reduction

    l_alt = get(handles.pu_mor_systems,'String');

    if ~isempty(l_alt)

        sOld=l_alt{get(handles.pu_mor_systems,'Value')};
        indexOld=find(strcmp(sOld,l));

        if ~isempty(indexOld)
            set(handles.pu_mor_systems,'Value',indexOld);
        else
            set(handles.pu_mor_systems,'Value',1);
        end

    else
      set(handles.pu_mor_systems,'Value',1);
    end
    
    
    %Loading and Setting up Models

    l_alt = get(handles.pu_LaSuM_truncate,'String');
    l_mimo = {};
    
    for i = 1:length(l)      %Show only MIMO-models
        try
            sysTemp = evalin('base',l{i,1});
            if size(sysTemp.B,2) > 1 || size(sysTemp.C,1) > 1
               l_mimo{end+1,1} = l{i,1}; 
            end
        end          
    end
    
    if ~isempty(l_alt)

        sOld=l_alt{get(handles.pu_LaSuM_truncate,'Value')};
        indexOld=find(strcmp(sOld,l_mimo));

        if ~isempty(indexOld)
            set(handles.pu_LaSuM_truncate,'Value',indexOld);
        else
            set(handles.pu_LaSuM_truncate,'Value',1);
        end

    else
      set(handles.pu_LaSuM_truncate,'Value',1);
    end
    
    if ~isempty(l_mimo) && size(l_mimo,1) >= 1
        set(handles.pu_LaSuM_truncate,'String',l_mimo);
        set(handles.pb_LaSuM_truncate,'Enable','on');
    else
        set(handles.pu_LaSuM_truncate,'String', [{''}; l_mimo]);
        set(handles.pb_LaSuM_truncate,'Enable','off');
    end
    
    
    %If there are no systems in workspace, set the list of selectable
    %options to the empty string, else set it to the list of systems

    if ~isempty(l) && size(l,1) >= 1
        
        set(handles.lb_PaV_systemsWs,'String',l);
        set(handles.pu_mor_systems,'String', l);
        set(handles.pu_an_sys1,'String',l);
        
        set(handles.virtgr_an_red_buttons,'Enable','on');
        
        if get(handles.pu_mor_method,'Value') ~= 3    %not Krylov
            set(handles.ed_mor_q,'Enable','on')
            set(handles.sl_mor_q,'Enable','on')
        end        
        
        pu_mor_systems_Callback(handles.pu_mor_systems, eventdata, handles);
        lb_PaV_systemsWs_Callback(handles.lb_PaV_systemsWs,eventdata,handles);
        pu_an_sys1_Callback(handles.pu_an_sys1,-1,handles);
        
    else
        set(handles.pu_mor_systems,'Value',1);  
        
        set(handles.lb_PaV_systemsWs,'String', [{''}; l]);
        set(handles.pu_mor_systems,'String', [{''}; l]);
        set(handles.pu_an_sys1,'String',[{''}; l]);

        set(handles.pb_mor_reduce,'Enable','off');
        set(handles.panel_mor_hsv,'Visible','off');
        set(handles.virtgr_an_red_buttons,'Enable','off');
        
        pu_mor_systems_Callback(handles.pu_mor_systems, eventdata, handles);
        lb_PaV_systemsWs_Callback(handles.lb_PaV_systemsWs,eventdata,handles);
        pu_an_sys1_Callback(handles.pu_an_sys1,eventdata,handles);
        pu_an_sys2_Callback(handles.pu_an_sys2,eventdata,handles);
    end
    
    %list vectors in workspace that might be s0
    
    set(handles.pu_mor_krylov_s0,'String',listS0InWorkspace);

function pu_mor_method_Callback(hObject, eventdata, handles)
% selection of reduction method
if get(hObject,'Value')==1
    set(handles.panel_mor_tbr,'Visible','on')
    set(handles.panel_mor_modal,'Visible','off')
    set(handles.panel_mor_krylov,'Visible','off')
    set(handles.panel_mor_hsv,'Visible','off')
elseif get(hObject,'Value')==2
    set(handles.panel_mor_tbr,'Visible','off')
    set(handles.panel_mor_modal,'Visible','on')
    set(handles.panel_mor_krylov,'Visible','off')
elseif get(hObject,'Value')==3
    set(handles.panel_mor_tbr,'Visible','off')
    set(handles.panel_mor_modal,'Visible','off')
    set(handles.panel_mor_krylov,'Visible','on')
end

%If TBR is selected, show the option for saving the hankel singular values

if get(hObject,'Value')==1
    set(handles.uipanel_mor_saveHsv,'Visible','on');
else
    set(handles.uipanel_mor_saveHsv,'Visible','off');
end

%If TBR is selected, check wether the HSV for this system have been
%calculated already
    
if get(handles.pu_mor_method,'Value')==1

    x = get(handles.pu_mor_systems,'String');
    y = x{get(handles.pu_mor_systems,'Value')};

    if ~isempty(y)
        try
           sys = evalin('base', y);

           if getHankelSingularValues(sys,y,handles)
               updateTBR(handles.pb_mor_hsv,eventdata,handles);
           elseif isa(sys,'sss') && ~isempty(sys.HankelSingularValues) && ...
                  ~isempty(sys.TBal) && ~isempty(sys.TBalInv)
               handles = saveHankelSingularValues(sys,y,sys.HankelSingularValues,sys.TBal,sys.TBalInv,handles);
               updateTBR(handles.pb_mor_hsv,eventdata,handles);   
           else
               set(handles.panel_mor_hsv,'Visible','off');
           end 
        catch ex
            set(handles.panel_mor_hsv,'Visible','off');
        end
    else
       set(handles.panel_mor_hsv,'Visible','off');
    end        
end

%If Krylov is selected, the reduced order can't be specified by the slider 
%anymore

if get(hObject,'Value')==3
    set(handles.ed_mor_q,'Enable','off')
    set(handles.ed_mor_q,'String','')
    set(handles.sl_mor_q,'Enable','off')
else
    set(handles.ed_mor_q,'Enable','on')
    set(handles.ed_mor_q,'String',num2str(get(handles.sl_mor_q,'Value')))
    set(handles.sl_mor_q,'Enable','on')
end

%If Krylov is selected, set the default values for the expension points and
%show the right panel dependent on if the system is siso or mimo

if get(hObject,'Value')==3

    x = get(handles.pu_mor_systems,'String');
    y = x{get(handles.pu_mor_systems,'Value')};
    
    if ~isempty(y)
       sys = evalin('base', y);
       
       if sys.m > 1 || sys.p > 1    %MIMO
           
           set(handles.panel_mor_krylov_MimoExps,'Visible','on');
           %set(handles.panel_mor_krylov_exps,'Visible','off');
           
           layoutMimoKrylov(handles);
           
       else                         %SISo
           
           set(handles.panel_mor_krylov_MimoExps,'Visible','off');
           %set(handles.panel_mor_krylov_exps,'Visible','on');
           
           layoutSisoKrylov(handles);
       end
    else
        
        set(handles.panel_mor_krylov_MimoExps,'Visible','off');
        %set(handles.panel_mor_krylov_exps,'Visible','on');
        
        layoutSisoKrylov(handles);
    end
    
    setKrylovDefaultValues(handles);
end

%If Krylov is selected, show the riht

% if get(hObject,'Value')==3
%    
%     value = get(handles.uitable_mor_krylov,'Data');
%    
%     if isempty(cell2mat(value(1,1))) && isempty(cell2mat(value(1,2)))
%         
%         x = get(handles.pu_mor_systems,'String');
%         y = x{get(handles.pu_mor_systems,'Value')};
%         
%         if ~isempty(y)
%            sys = evalin('base', y);
%            if size(sys.A,1) > 6
%                value{1,1} = 0;
%                value{1,2} = 2;
%                value{2,1} = 1+i;
%                value{2,2} = 2;
%                value{3,1} = 1-i;
%                value{3,2} = 2;           
%            else
%                value{1,1} = 0;
%                value{1,2} = 1;
%            end
%            set(handles.uitable_mor_krylov,'Data',value);
%            set(handles.uitable_mor_krylov_output,'Data',value);
%            countMatchedMoments(handles);
%         end        
%     end    
% end

%Set the "Save Shifts" option visible if krylov and irka is selected

if get(hObject,'Value')==3 && get(handles.pu_mor_krylov_algorithm,'Value')==1
    set(handles.uipanel_mor_optimalShifts,'Visible','on');
else
    set(handles.uipanel_mor_optimalShifts,'Visible','off');
end

%Change the suggested Names for the reducted system

x = get(handles.pu_mor_systems,'String');
y = x{get(handles.pu_mor_systems,'Value')};

suggestNamesMOR(y,handles);

function sl_mor_q_Callback(hObject, eventdata, handles)
    % order selection slider has been changed
    q=get(hObject,'Value'); 
    % round to integer number
    q=round(q);
    set(hObject,'Value',q)
    % write to textfield
    set(handles.ed_mor_q,'String',q);
    updateTBR(hObject, eventdata, handles)

function ed_mor_q_Callback(hObject, eventdata, handles)
% order selection textfield has been changed
sys_x=get(handles.pu_mor_systems,'String');
sysname=sys_x{get(handles.pu_mor_systems,'Value')};
if isempty(sysname)
    errordlg('Please choose a system first!','Error Dialog','modal')
    uiwait
    set(hObject,'String',get(handles.sl_mor_q,'Value'))
    return
end
x=str2num(get(hObject,'String')); %#ok<ST2NM>
x=round(x);
if  isempty(x) || length(x)>1 || x<1 || imag(x)~=0
    errordlg('Has to be an integer','Error Dialog','modal')
    uiwait
    set(hObject,'String',get(handles.sl_mor_q,'Value'))
    return
elseif get(handles.sl_mor_q,'Max')<x
    errordlg('Must be smaller than the size of the original system','Error Dialog','modal')
    uiwait
    set(hObject,'String',get(handles.sl_mor_q,'Value'))
    return
end
set(hObject,'String',x)
set(handles.sl_mor_q,'Value',x)
updateTBR(hObject, eventdata, handles)

%Callbacks of the controls on the bottom

function pb_mor_reduce_Callback(hObject, eventdata, handles)

    %Check if ths spectified names for the system and for the projection
    %matrices are correct

    if (get(handles.ed_mor_w,'UserData')==1 && get(handles.cb_mor_savew,'Value')==1)...
     ||(get(handles.ed_mor_v,'UserData')==1 && get(handles.cb_mor_savev,'Value')==1)
        errordlg('Please correct names for projection matrices first','Error Dialog','modal')
        uiwait
        return
    end
    if get(handles.cb_mor_saveShifts,'Value')==1 && get(handles.pu_mor_krylov_algorithm,'Value')==1 && ...
            get(handles.pu_mor_method,'Value')==3 && get(handles.cb_mor_saveShifts,'Value')==1 && ...
            get(handles.ed_mor_saveShifts,'UserData')==1
        errordlg('Please correct names for the optimal shifts first','Error Dialog','modal')
        uiwait
        return
    end    
    if get(handles.cb_mor_saveHsv,'Value')==1 && get(handles.pu_mor_method,'Value')==1 && ...
            get(handles.cb_mor_saveHsv,'Value')==1 && get(handles.ed_mor_saveHsv,'UserData')==1
        errordlg('Please correct names for the Hankel Singular Values','Error Dialog','modal')
        uiwait
        return
    end
    q=str2double(get(handles.ed_mor_q,'String'));
    if (isnan(q) || isempty(q)) && get(handles.pu_mor_method,'Value')~=3
            errordlg('Please choose order of reduced system first.','Error Dialog','modal')
            uiwait
            return
    end
    if get(handles.ed_mor_sysred,'UserData')==1
        errordlg('Please correct name for reduced system first','Error Dialog','modal')
        uiwait
        return
    elseif existInBaseWs(get(handles.ed_mor_sysred,'String'))==1
        s=sprintf('%s already exists in base workspac. Do you want to overwrite it?',get(handles.ed_mor_sysred,'String'));
        % show dialog box
        k=stqd('String',s,'Title','Question Dialog');
        if ~isempty(k) && strcmp(k,'No')
            return
        end
    end
        
    %Get the system that should be reduced
    
    set(handles.figure1,'Pointer','watch')
    set(hObject,'Enable','off')
    drawnow % this makes the watch appear
    sysname = get(handles.pu_mor_systems,'String');
    sysname = sysname{get(handles.pu_mor_systems,'Value')};
    try
       sys=evalin('base',sysname);
    catch ex
       errordlg('Selected system could not be found in the workspace.')
       uiwait
       pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles)
       return
    end

    % Reduce
    
    switch get(handles.pu_mor_method,'Value')
        
    case 1 %TBR
        
        [hsvStored,hsv,R,L] = getHankelSingularValues(sys,sysname,handles);
        
        if ~hsvStored || isempty(hsv) || isempty(R) || isempty(L)
            set(handles.figure1,'Pointer','arrow')
            set(hObject,'Enable','on')
            errordlg('Please calculate Hankel Singular Values first.','Error Dialog','modal')
            uiwait
            return
        end
        
        %Reduce
        
        lastwarn('');
        
        try
        
            if get(handles.rb_mor_tbrtruncate,'Value')==1
                % truncation

                Opts.type = 'tbr';

                [sysr, V, W] = tbr(sys,q,R,L,Opts);

            else
                % match DC gain

                Opts.type  = 'matchDcGain';

                [sysr, V, W] = tbr(sys,q,R,L,Opts);

            end
            
            error('TBR:WarningOccured',lastwarn);
        
        catch ex
            set(handles.figure1,'Pointer','arrow')
            set(hObject,'Enable','on')
            if strcmp(ex.identifier,'TBR:WarningOccured')
                if ~isempty(ex.message)
                    msgbox(ex.message,'Warning','warn');
                    uiwait
                end
            else
                % display the line of the error if "tbr" produced the error
                [found,str] = extractEntryFromErrorStack(ex.stack,'tbr');
                if found            
                    errordlg({str;ex.message},'Error Dialog','modal');
                else
                    errordlg(ex.message,'Error Dialog','modal')
                end
                uiwait
                return
            end
        end
        
    case 2 % modal truncation
        
        %Check if the value specified for neighborhood is correct if the
        %option "in neighborhood of" is specified
        
        if get(handles.rb_mor_modal_neig,'Value')==1
            if get(handles.ed_mor_modal_neig,'UserData')==1
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg('Please correct s0 first','Error Dialog','modal')
                uiwait
                return
            end
        end
        
        %Check if the A-matrix of the system is symmetric if SA or LA is
        %selected
        
        if get(handles.rb_mor_modal_LA,'Value')==1 || ...
                get(handles.rb_mor_modal_SA,'Value')==1
          
            if issymmetric(sys.A) == 0 || isreal(sys.A) == 0
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg('Options "LA" and "SA" are only valid for systems with real symmetric A!','Error Dialog','modal')
                uiwait
                return
            end            
        end
        
        %Options for eigenvalues
        
        if get(handles.rb_mor_modal_neig,'Value')==1
            opts.type = str2double(get(handles.ed_mor_modal_neig,'String'));
        elseif get(handles.rb_mor_modal_large,'Value')==1
            opts.type ='LM';
        elseif get(handles.rb_mor_modal_small,'Value')==1
            opts.type ='SM';
        elseif get(handles.rb_mor_modal_LA,'Value')==1
            opts.type ='LA';
        elseif get(handles.rb_mor_modal_SA,'Value')==1
            opts.type ='SA';
        elseif get(handles.rb_mor_modal_LR,'Value')==1
            opts.type ='LR';
        elseif get(handles.rb_mor_modal_SR,'Value')==1
            opts.type ='SR';
        end
        
        %Options for orthogonalization
        
        if get(handles.rb_mor_modal_orth_y,'Value')==1
            opts.orth = 'qr';
        else
            opts.orth = '0';
        end
        
        %Options for real reduced system
        
        if get(handles.rb_mor_modal_real_y,'Value')==1
            opts.real = 'real';
        else
            opts.real = '0'; 
        end
        
        %Reduce system

        lastwarn('');
        
        try
            [sysr, V, W] = modalMor(sys, q, opts);
            %sysr.mor_info=struct('time', clock, 'method', 'Modal Truncation', 'orgsys', sysname);
            error('MODAL:WarningOccured',lastwarn);
        catch ex
            set(handles.figure1,'Pointer','arrow')
            set(hObject,'Enable','on')
            if strcmp(ex.identifier,'MODAL:WarningOccured')
                if ~isempty(ex.message)
                    msgbox(ex.message,'Warning','warn');
                    uiwait
                end
            else
                % display the line of the error if "modalMor" produced the 
                % error
                [found,str] = extractEntryFromErrorStack(ex.stack,'modalMor');
                if found            
                    errordlg({str;ex.message},'Error Dialog','modal');
                else
                    errordlg(ex.message,'Error Dialog','modal')
                end
                uiwait
                return
            end
        end
        
        
    case 3 % Krylov
        
        
        %Get expansion points
        try
            if sys.m > 1 || sys.p > 1       %MIMO
                [s_inp,s_out,Rt,Lt] = getMimoExpensionPoints(handles);
            else                            %SISO
                [s_inp,s_out] = getExpansionPoints(handles);
            end
        catch ex
           errordlg(ex.message)
           uiwait
           set(handles.figure1,'Pointer','arrow')
           set(hObject,'Enable','on')
           return
        end
        
        if isempty(s_inp) && isempty(s_out)
            errordlg('The specified shifts and moments have the wrong format. Please correct this first!')
            uiwait
            set(handles.figure1,'Pointer','arrow')
            set(hObject,'Enable','on')
            return
        end
        
        %Reduce with the specified Algorithm
        
        switch get(handles.pu_mor_krylov_algorithm,'Value')
            
        case 2 % explicit moment matching
                        
            lastwarn('');

            try
                
                if (sys.m > 1 || sys.p > 1) && handles.MimoParam.block == 0  %MIMO
                    if get(handles.rb_mor_krylov_input,'Value')==1
                        [sysr, V, W] = rk(sys, s_inp,Rt);
                    else
                        [sysr, V, W] = rk(sys, s_inp, s_out, Rt, Lt);
                    end
                else                        %SISO or block-krylov
                    if get(handles.rb_mor_krylov_input,'Value')==1
                        [sysr, V, W] = rk(sys, s_inp);
                    else
                        [sysr, V, W] = rk(sys, s_inp, s_out);
                    end
                end
                
                error('RK:WarningOccured',lastwarn);
                
            catch ex
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                if strcmp(ex.identifier,'RK:WarningOccured')
                    if ~isempty(ex.message)
                        msgbox(ex.message,'Warning','warn');
                        uiwait
                    end
                else
                    % display the line of the error if "rk" produced the 
                    % error
                    [found,str] = extractEntryFromErrorStack(ex.stack,'rk');
                    if found            
                        errordlg({str;ex.message},'Error Dialog','modal');
                    else
                        errordlg(ex.message,'Error Dialog','modal')
                    end
                    uiwait
                    return
                end
            end

        case 3 %ICOP
            
            if sys.isMimo
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg('MIMO is not supported yet, coming soon','Error Dialog','modal')
                uiwait
                return
            end
            if get(handles.rb_mor_krylov_twosided,'Value')==1
                projection='two';
            elseif get(handles.rb_mor_krylov_twosided,'Value')==1
                projection='out';
            else
                projection='in';
            end
            maxiter=str2double(get(handles.ed_mor_krylov_max,'String'));
            epsilon=str2double(get(handles.ed_mor_krylov_epsilon,'String'));
            if isempty(maxiter) || maxiter==0 || isempty(epsilon) %|| epsilon==0 ??
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg('Please choose the maximum number of iterations and epsilon first','Error dialog','modal')
                uiwait
                return
            end
            
            %Reduce
            
            lastwarn('');
            
            try
                [sysr,V,W,alpha_opt] = RK_ICOP(sys, s0,'maxiter',maxiter,'epsilon',epsilon,'projection',projection );
                %sysr.mor_info=struct('time', clock, 'method', 'Krylov, ICOP', 'orgsys', sysname);
                %sysr.mor_info.alpha_opt=alpha_opt;
                error('ICOP:WarningOccured',lastwarn);
            catch ex
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                
                if strcmp(ex.identifier,'ICOP:WarningOccured')
                    if ~isempty(ex.message)
                        msgbox(ex.message,'Warning','warn');
                        uiwait
                    end
                else
                    % display the line of the error if "RK_ICOP" produced the 
                    % error
                    [found,str] = extractEntryFromErrorStack(ex.stack,'RK_ICOP');
                    if found            
                        errordlg({str;ex.message},'Error Dialog','modal');
                    else
                        errordlg(ex.message,'Error Dialog','modal')
                    end
                    uiwait
                    return
                end
            end
            
            if isempty(sysr)
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg('ICOP did not converge, reduction failed','Error dialog','modal')
                uiwait
                return
            end

        case 1 %IRKA
            
            %Options
            
            if get(handles.pu_mor_krylov_StopCrit,'Value') == 3
                Opts.stopCrit = 'combAny';
            elseif get(handles.pu_mor_krylov_StopCrit,'Value') == 1
                Opts.stopCrit = 's0';
            elseif get(handles.pu_mor_krylov_StopCrit,'Value') == 2
                Opts.stopCrit = 'sysr';
            else
                Opts.stopCrit = 'combAll';
            end
            
            Opts.maxiter=str2double(get(handles.ed_mor_krylov_max,'String'));
            Opts.epsilon=str2double(get(handles.ed_mor_krylov_epsilon,'String'));
            if isnan(Opts.maxiter) || Opts.maxiter<=0 || isnan(Opts.epsilon) || Opts.epsilon<=0
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg({'Please correct the maximum number of iterations and epsilon first!', ...
                        'Both values have to be real numbers greater than zero.'},'Error dialog','modal')
                uiwait
                return
            end
            
            %Reduce
            
            lastwarn('');
            
            try
                if (sys.m > 1 || sys.p > 1) && handles.MimoParam.block == 0  %MIMO
                
                    [sysr,V,W,s0] = irka(sys, s_inp, Rt, Lt, Opts);
                
                else                                %SISO or Block-Krylov
                    
                    [sysr,V,W,s0] = irka(sys, s_inp, Opts);
                    
                end

                error('IRKA:WarningOccured',lastwarn);
            catch ex                
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                
                if strcmp(ex.identifier,'IRKA:WarningOccured')
                    if ~isempty(ex.message)
                        msgbox(ex.message,'Warning','warn');
                        uiwait
                    end
                else
                    % display the line of the error if "irka" produced the 
                    % error
                    [found,str] = extractEntryFromErrorStack(ex.stack,'irka');
                    if found            
                        errordlg({str;ex.message},'Error Dialog','modal');
                    else
                        errordlg(ex.message,'Error Dialog','modal')
                    end
                    uiwait
                    return
                end
            end
        end
    end

    %Write system to workspace
    
    if get(handles.cb_mor_krylov,'Value')==1
        % impose E_r=I?
        sysr = sysr.resolveDescriptor;
    end

    assignin('base',get(handles.ed_mor_sysred,'String'),sysr)
    
    %Write Projection matrices and optimal shifts(IRKA) to workspace
    
    if get(handles.cb_mor_savew,'Value')==1
        assignin('base',get(handles.ed_mor_w,'String'),W)
    end
    
    if get(handles.cb_mor_savev,'Value')==1
        assignin('base',get(handles.ed_mor_v,'String'),V)
    end
    
    if get(handles.cb_mor_saveShifts,'Value')==1 && get(handles.pu_mor_krylov_algorithm,'Value')==1 && ...
            get(handles.pu_mor_method,'Value')==3
        assignin('base',get(handles.ed_mor_saveShifts,'String'),s0); 
    end
    
    if get(handles.cb_mor_saveHsv,'Value')==1 && get(handles.pu_mor_method,'Value')==1
        assignin('base',get(handles.ed_mor_saveHsv,'String'),hsv); 
    end

    %Tell the user that the reduction was successfull
    
    set(handles.figure1,'Pointer','arrow')
    msgbox('Reduction was successful!','Information','modal')
    uiwait
    
    %Update the displays of systems in workspace
    
    pb_refreshsys_Callback(hObject, eventdata, handles)
    pb_refreshlb_Callback(hObject, eventdata, handles)
    pb_mor_krylov_refresh_Callback(handles.pb_mor_krylov_refresh, eventdata, handles)
    
    %Set the optimal shifts as default for vector-import from workspace
    
    if get(handles.cb_mor_saveShifts,'Value')==1 && get(handles.pu_mor_krylov_algorithm,'Value')==1
        name = get(handles.ed_mor_saveShifts,'String');
        x = get(handles.pu_mor_krylov_s0,'String');
        for i = 1:length(x)
           if strcmp(x{i,1},name)
               set(handles.pu_mor_krylov_s0,'Value',i);
               break;
           end
        end
    end
    
    %Suggest names for next reduction
    
    s=get(handles.ed_mor_sysred,'String');
    s=s(~isstrprop(s,'digit'));
    s=suggestVarname(s,handles.ed_mor_sysred);
    suggestVarname(sprintf('%s_w',s),handles.ed_mor_w);
    suggestVarname(sprintf('%s_v',s),handles.ed_mor_v);
    set(hObject,'Enable','on')

function ed_mor_w_Callback(hObject, eventdata, handles)
% check variable name
isvalidvarname(hObject)

function ed_mor_v_Callback(hObject, eventdata, handles)
% check variable name
isvalidvarname(hObject)

function ed_mor_sysred_Callback(hObject, eventdata, handles)
% check variable name
isvalidvarname(hObject)

function ed_mor_saveHsv_Callback(hObject, eventdata, handles)
    % check variable name
    isvalidvarname(hObject)
    
function ed_mor_saveShifts_Callback(hObject, eventdata, handles)
% check variable name
    isvalidvarname(hObject)

function cb_mor_savew_Callback(hObject, eventdata, handles)
if get(hObject,'Value') == 1
    set(handles.ed_mor_w,'Enable','on')
else
    set(handles.ed_mor_w,'Enable','off')
end

function cb_mor_savev_Callback(hObject, eventdata, handles)
if get(hObject,'Value') == 1
    set(handles.ed_mor_v,'Enable','on')
else
    set(handles.ed_mor_v,'Enable','off')
end

function cb_mor_saveHsv_Callback(hObject, eventdata, handles)

    if get(hObject,'Value') == 1
        set(handles.ed_mor_saveHsv,'Enable','on')
    else
        set(handles.ed_mor_saveHsv,'Enable','off')
    end

function cb_mor_saveShifts_Callback(hObject, eventdata, handles)
 
    if get(hObject,'Value') == 1
        set(handles.ed_mor_saveShifts,'Enable','on')
    else
        set(handles.ed_mor_saveShifts,'Enable','off')
    end

%Callbacks for Balacing & Truncation

function updateTBR(hObject, eventdata, handles)

    %Show wait-cursor
    
    set(handles.figure1,'Pointer','watch')
    drawnow

    %Get the selected system

    sys_x=get(handles.pu_mor_systems,'String');
    sysname=sys_x{get(handles.pu_mor_systems,'Value')};
    try
        sys=evalin('base',sysname);
    catch
       errordlg('Selected system could not be found in the workspace.')
       uiwait
       pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles)
       return
    end
        
        
    %Get Axes and desired model order
    
    hr=get(handles.axes_mor_hsv,'UserData');
    q=get(handles.sl_mor_q,'Value');
    
    if mod(q,1) ~= 0
        q = 1;
        set(handles.sl_mor_q,'Value',1);
    end
    
    [hsvStored,hsv] = getHankelSingularValues(sys,sysname,handles);
    
    if get(handles.pu_mor_method,'Value')==1 && hsvStored
        %Calculate the signal-norms H_1 and H_inf
        
        e=2*sum(hsv((q+1):end));
        erel=e/max(hsv);
        set(handles.st_mor_tbr_error,'String',num2str(e, '%.3e'))
        set(handles.st_mor_tbr_relerror,'String',num2str(erel, '%.3e'))
        
        %Plot which shows the Hankel-Singular-Values
            
        set(handles.panel_mor_hsv, 'Visible','on') 
        cla(handles.axes_mor_hsv)
        hold(handles.axes_mor_hsv, 'on')

        if get(handles.rb_mor_tbr_norm,'Value')==1

            maxValue = max(hsv);
            h = plot(handles.axes_mor_hsv,hsv./maxValue);

        else
            h = plot(handles.axes_mor_hsv, hsv);
        end

        % make callback react to click on red HSV line

        set(h,'HitTest','off');

        % legend

        legend(handles.axes_mor_hsv, regexprep(sysname, '_', ' '));

        % set scale

        if get(handles.rb_mor_tbr_log,'Value')==1 || ...
           get(handles.rb_mor_tbr_norm,'Value')==1
            set(handles.axes_mor_hsv,'YScale','log')
        else
            set(handles.axes_mor_hsv,'YScale','linear')                
        end
        
        %Plot the red line
            
        if ishandle(hr)
            set(hr,'XData',[q,q])
            set(hr,'YData',get(handles.axes_mor_hsv,'YLim'));
        else
            hr=plot(handles.axes_mor_hsv,[q,q], get(handles.axes_mor_hsv,'YLim'),'r');
            set(handles.axes_mor_hsv,'UserData',hr)
        end
    end
    
    %Set cursor back to arrow
    
    set(handles.figure1,'Pointer','arrow')

function axes_mor_hsv_ButtonDownFcn(hObject, eventdata, handles)
% refresh red line in HSV axes and update q-slider and q-textfield
x=get(hObject,'CurrentPoint');
% use only x value and round it
q=round(x(1,1));
% control boundary of slider
if q<get(handles.sl_mor_q,'Min')
    q=get(handles.sl_mor_q,'Min');
elseif q>get(handles.sl_mor_q,'Max')
    q=get(handles.sl_mor_q,'Max');
end
set(handles.ed_mor_q,'String',q)
set(handles.sl_mor_q,'Value',q)

updateTBR(hObject, eventdata, handles)

function pb_mor_hsv_Callback(hObject, eventdata, handles)
% calculate and plot Hankel Singular Values

% disable to avoid double call
set(hObject,'Enable','off')
set(handles.figure1,'Pointer','watch')
drawnow

% Get the system from the workspace

try
    [sys, sysname] = getSysFromWs(handles.pu_mor_systems);
catch ex
    set(handles.figure1,'Pointer','arrow')
    set(hObject,'Enable','on')
    if strfind(ex.identifier, 'unassigned')
        errordlg('Please choose a system first','Error Dialog','modal')
        return
    end
    errordlg(ex.message,'Error Dialog','modal')
    uiwait
    return
end

%Calculate the Hankel-Singular-Values

try 
    Opts.rctol      = 1e-3;
    [~,~,~,hsv,R,L] = tbr(sys,1,Opts);
    handles = saveHankelSingularValues(sys,sysname,hsv,R,L,handles);
catch ex %***
    if strcmp(ex.identifier,'MATLAB:nomem')
        errordlg('Out of memory. System is too large to calculate Hankel Singular Values.','Error Dialog','modal')
    elseif strcmp(ex.identifier,'Control:foundation:LyapChol4')
        errordlg('A or (A,E) must have all their eigenvalues in the left-half plane.','Error Dialog','modal')
    else
        [found,str] = extractEntryFromErrorStack(ex.stack,'tbr');
        if found            
            errordlg({str;ex.message},'Error Dialog','modal');
        else
            errordlg(ex.message,'Error Dialog','modal');
        end
    end

    uiwait
    set(handles.figure1,'Pointer','arrow')
    set(hObject,'Enable','on')
    return
end
set(handles.panel_mor_hsv,'Visible','off')

%Plot the Hankel-Singular-Values

updateTBR(hObject,eventdata,handles);

set(handles.figure1,'Pointer','arrow')
set(hObject,'Enable','on')

function bg_mor_plot_SelectionChangeFcn(hObject, eventdata, handles) 

    updateTBR(hObject,eventdata,handles);


function pb_mor_hsv_zoomin_Callback(hObject, eventdata, handles)
% zoom +
if ~isa(handles.zoom,'graphics.zoom')
    % do not zoom logos
    handles.zoom=zoom(handles.figure1);
    setAllowAxesZoom(handles.zoom,handles.logo_about,false)
	setAllowAxesZoom(handles.zoom,handles.logo_tum,false)
end
set(handles.zoom,'Motion','both') % horizontally und vertically
set(handles.zoom,'Direction','in') % zoom in on click
set(handles.zoom,'Enable','on') % activate zoom
guidata(hObject,handles)

function pb_mor_hsv_zoomout_Callback(hObject, eventdata, handles)
%zoom -
if ~isa(handles.zoom,'graphics.zoom')
    % do not zoom logos
    handles.zoom=zoom(handles.figure1);
    setAllowAxesZoom(handles.zoom,handles.logo_about,false)
	setAllowAxesZoom(handles.zoom,handles.logo_tum,false)
end
set(handles.zoom,'Motion','both') % horizontally und vertically
set(handles.zoom,'Direction','out') % zoom out on click
set(handles.zoom,'Enable','on') % activate zoom
guidata(hObject,handles)

function pb_mor_hsv_zoomoff_Callback(hObject, eventdata, handles)
% zoom off
if isa(handles.zoom,'matlab.graphics.interaction.internal.zoom')
    set(handles.zoom,'Enable','off') % deactivate zoom
end
guidata(hObject,handles)


function pb_info_tbr_Callback(hObject, eventdata, handles)

    infoBox({'InfoErrorBound.png'});
    uiwait;

%Callbacks for Modal

function ed_mor_modal_neig_Callback(hObject, eventdata, handles)
% check of double. replace ',' with '.'
x=str2num(get(hObject,'String')); %#ok<ST2NM>
if isempty(x) || length(x)>1
    if length(x)==2
        x=sprintf('%i.%i',round(x(1)),round(x(2)));
        set(hObject,'String',x)
        set(hObject,'UserData',0)
        return
    end
    errordlg('Has to be a real or complex number.','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
else
    set(hObject,'UserData',0)
end

%Callbacks for Krylov

function pu_mor_krylov_algorithm_Callback(hObject, eventdata, handles)

    %Set the values for the elements dependent on SISO or MIMO

    if strcmp(get(handles.panel_mor_krylov_MimoExps,'Visible'),'off')

        if get(hObject,'Value')==1
            %IRKA
            set(handles.cb_mor_krylov_hermite,'Value',1)
            set(handles.cb_mor_krylov_hermite,'Enable','off')
            set(handles.pb_mor_krylov_input,'Visible','off')
            set(handles.pb_mor_krylov_output,'Visible','off')
            set(handles.uitable_mor_krylov_output,'Visible','off')
            set(handles.uitable_mor_krylov,'Visible','on')
        else %explicit moment matching
            set(handles.cb_mor_krylov_hermite,'Enable','on')
            set(handles.cb_mor_krylov_hermite,'Value',1)       
        end
    
    else
       layoutMimoKrylov(handles); 
    end
    
    % if get(hObject,'Value')==2
    %     % ICOP, only sigle expansion point and Arnoldi
    %     set(handles.rb_mor_krylov_single,'Value',1)
    %     set(handles.rb_mor_krylov_single,'Enable','inactive')
    %     set(handles.rb_mor_krylov_multi,'Enable','off')
    %     set(handles.pu_mor_krylov_alg,'Visible','off')
    % else
    %     set(handles.rb_mor_krylov_single,'Enable','on')
    %     set(handles.rb_mor_krylov_multi,'Enable','on')
    %     set(handles.pu_mor_krylov_alg,'Visible','on')
    % end
    
    %Set the values for elements indepentent from SISO or MIMO
    
    if get(hObject,'Value')==2
        % explicit moment matching
        set(handles.rb_mor_krylov_twosided,'Enable','on')
        set(handles.rb_mor_krylov_input,'Enable','on')
        set(handles.rb_mor_krylov_output,'Enable','on')
        
        set(handles.panel_mor_krylov_iter,'Visible','off')
        set(handles.st_mor_krylov_AlgoParams,'Visible','off')
        set(handles.pb_mor_krylov_infoAlgoParams,'Visible','off')
        
        set(handles.st_mor_krylov_points,'String','2. Expansion points')
    else
        %IRKA
        set(handles.rb_mor_krylov_twosided,'Value',1)
        set(handles.rb_mor_krylov_twosided,'Enable','inactive')
        set(handles.rb_mor_krylov_input,'Enable','off')
        set(handles.rb_mor_krylov_output,'Enable','off')
        
        set(handles.panel_mor_krylov_iter,'Visible','on')
        set(handles.st_mor_krylov_AlgoParams,'Visible','on')
        set(handles.pb_mor_krylov_infoAlgoParams,'Visible','on')
        
        set(handles.st_mor_krylov_points,'String','2. Starting points')
    end    
    
    %Update matched Moments and reduced order

    countMatchedMoments(handles);
    
    %Set the panel for saving the optimal shifts visible if irka is
    %selected
    
    if get(hObject,'Value')==1
        set(handles.uipanel_mor_optimalShifts,'Visible','on');
    else
        set(handles.uipanel_mor_optimalShifts,'Visible','off');
    end
    
    %Change the suggested Names for the reducted system

    x = get(handles.pu_mor_systems,'String');
    y = x{get(handles.pu_mor_systems,'Value')};

    suggestNamesMOR(y,handles);  

function bg_mor_krylov_side_SelectionChangedFcn(hObject, eventdata, handles)
    
    if strcmp(get(handles.panel_mor_krylov_MimoExps,'Visible'),'off')

        layoutSisoKrylov(handles);
        countMatchedMoments(handles);    
    else       
        layoutMimoKrylov(handles);        
    end
    
    
    
function ed_mor_krylov_max_Callback(hObject, eventdata, handles)
% maximum iterations
x=str2double(get(hObject,'String'));
x=round(abs(real(x)));
set(hObject,'String',x)

function ed_mor_krylov_epsilon_Callback(hObject, eventdata, handles)
% epsilon
x=str2double(get(hObject,'String'));
x=abs(real(x));
set(hObject,'String',x)



function pb_mor_krylov_add_Callback(hObject, eventdata, handles) 
% add empty row to table
if strcmp(get(handles.uitable_mor_krylov,'Visible'),'on')
    x=get(handles.uitable_mor_krylov,'Data');
    x{size(x,1)+1,1}=[];
    set(handles.uitable_mor_krylov,'Data',x)
else
    x=get(handles.uitable_mor_krylov_output,'Data');
    x{size(x,1)+1,1}=[];
    set(handles.uitable_mor_krylov_output,'Data',x)
end

function pb_mor_krylov_reset_Callback(hObject, eventdata, handles)

if strcmp(get(handles.uitable_mor_krylov,'Visible'),'on')
    set(handles.uitable_mor_krylov,'Data',cell(1,2))
else
    set(handles.uitable_mor_krylov_output,'Data',cell(1,2))
end

countMatchedMoments(handles);

function pb_mor_krylov_AddCompConj_Callback(hObject, eventdata, handles)

%Get the indices of the selected cells

indices = handles.SelectedIndex;

%Add the complex conjugate of the selected cells as new rows

if strcmp(get(handles.uitable_mor_krylov,'Visible'),'on')
    x=get(handles.uitable_mor_krylov,'Data');
    xConj = zeros(size(indices,1),2);
    
    for i = 1:size(indices,1)      
       xConj(i,1) = conj(x{indices(i,1),1});
       xConj(i,2) = x{indices(i,1),2};
    end
    
    x = [x;num2cell(xConj)];
    set(handles.uitable_mor_krylov,'Data',x);
else
    x=get(handles.uitable_mor_krylov_output,'Data');
    xConj = zeros(size(indices,1),2);
    
    for i = 1:size(indices,1)      
       xConj(i,1) = conj(x{indices(i,1),1});
       xConj(i,2) = x{indices(i,1),2};
    end
    
    x = [x;num2cell(xConj)];
    set(handles.uitable_mor_krylov_output,'Data',x);
end

countMatchedMoments(handles);

function pb_mor_krylov_deleteRow_Callback(hObject, eventdata, handles)

%Read out the selected indices from the handles

indices = handles.SelectedIndex;

%Delete selected rows

if strcmp(get(handles.uitable_mor_krylov,'Visible'),'on')
    x=get(handles.uitable_mor_krylov,'Data');
    
    for i = 1:size(indices,1)
       x(indices(i,1)-i+1,:) = []; 
    end
    
    set(handles.uitable_mor_krylov,'Data',x)
else
    x=get(handles.uitable_mor_krylov_output,'Data');
    
    for i = 1:size(indices,1)
       x(indices(i,1)-i+1,:) = []; 
    end
    
    set(handles.uitable_mor_krylov_output,'Data',x)
end

countMatchedMoments(handles);



function pu_mor_krylov_s0_CreateFcn(hObject, eventdata, handles)
% list of vectors in workspace that might be expansion points
set(hObject,'String',listS0InWorkspace)

function ed_mor_krylov_exp_Callback(hObject, eventdata, handles)
% expansion point must be real scalar
x=str2num(get(hObject,'String')); %#ok<ST2NM>
x=real(x);
if length(x)==2
    x=str2double(sprintf('%i.%i',round(x(1)),round(x(2))));
elseif length(x)>2
    x=str2double(get(hObject,'String'));
    x=real(x);
elseif isempty(x)
    x=0;
end
set(hObject,'String',x)

function cb_mor_krylov_hermite_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 0
    
    set(handles.uitable_mor_krylov_output,'Visible','Off');
    set(handles.uitable_mor_krylov,'Visible','On');

    set(handles.pb_mor_krylov_input,'Visible','On');
    set(handles.pb_mor_krylov_output,'Visible','On');

    set(handles.pb_mor_krylov_input,'BackgroundColor',[1;0.843;0]);
    set(handles.pb_mor_krylov_output,'BackgroundColor',[0.94;0.94;0.94]);
    
else
    
    set(handles.uitable_mor_krylov_output,'Visible','Off');
    set(handles.uitable_mor_krylov,'Visible','On');

    set(handles.pb_mor_krylov_input,'Visible','Off');
    set(handles.pb_mor_krylov_output,'Visible','Off');
    
end

countMatchedMoments(handles);

function et_mor_krylov_commandLine_Callback(hObject, eventdata, handles)

    %Execute the typed in command in the workspace

    text = get(handles.et_mor_krylov_commandLine,'String');

    try
        m = evalin('base',text);
    catch er
        errordlg(er.message);
        return;
    end

    %Read out the expension points from the vector

    if (size(m,1) == 1 || size(m,2) == 1) || correctExpensionPointMatrix(m) == 1
        data = addMatrixToTable(m);
    else
        errordlg('The matrix created with the command line has not the correct format.');
        return;
    end

    %Write the data to the table

    if strcmp(get(handles.uitable_mor_krylov,'Visible'),'on')
        set(handles.uitable_mor_krylov,'Data',data);
    else
        set(handles.uitable_mor_krylov_output,'Data',data);
    end
    
    countMatchedMoments(handles);

function pb_mor_krylov_refresh_Callback(hObject, eventdata, handles)
% list of vectors in workspace that might be expansion points
set(handles.pu_mor_krylov_s0,'String',listS0InWorkspace)

function pb_mor_krylov_importVector_Callback(hObject, eventdata, handles)
    % if vector from workspace was selected, add it to table
    contents = cellstr(get(handles.pu_mor_krylov_s0,'String'));
    try
        x=evalin('base',contents{get(handles.pu_mor_krylov_s0,'Value')});
    catch ex
       errordlg('No valid matrix selected');
       pb_mor_krylov_refresh_Callback(handles.pb_mor_krylov_refresh,eventdata,handles);
       return;
    end
    
    tableData = addMatrixToTable(x);
 
    if strcmp(get(handles.uitable_mor_krylov,'Visible'),'on')
        set(handles.uitable_mor_krylov,'Data',tableData);
    else
        set(handles.uitable_mor_krylov_output,'Data',tableData);
    end
    
    countMatchedMoments(handles);   

    
    
function uitable_mor_krylov_CellEditCallback(hObject, eventdata, handles)
    
countMatchedMoments(handles);

function uitable_mor_krylov_CreateFcn(hObject, eventdata, handles)
% adapt size of table
set(hObject,'Data',cell(1,2))
handles.SelectedIndex = 1;
guidata(hObject,handles);

function uitable_mor_krylov_output_CellEditCallback(hObject, eventdata, handles)

countMatchedMoments(handles);

function uitable_mor_krylov_output_CreateFcn(hObject, eventdata, handles)
% adapt size of table
set(hObject,'Data',cell(1,2))
handles.SelectedIndex = 1;
guidata(hObject,handles);

function uitable_mor_krylov_CellSelectionCallback(hObject, eventdata, handles)
%Reads out the selected indices and stores it in a the handles structure
currentCell = eventdata.Indices;
currentCell = currentCell(:,1);
currentCell = unique(currentCell);
handles.SelectedIndex = currentCell;
guidata(hObject,handles);

function uitable_mor_krylov_output_CellSelectionCallback(hObject, eventdata, handles)
%Reads out the selected indices and stores it in a the handles structure
currentCell = eventdata.Indices;
currentCell = currentCell(:,1);
currentCell = unique(currentCell);
handles.SelectedIndex = currentCell;
guidata(hObject,handles);

function pb_mor_krylov_input_Callback(hObject, eventdata, handles)

%Set the right table visible

set(handles.uitable_mor_krylov,'Visible','On');
set(handles.uitable_mor_krylov_output,'Visible','Off');

%Change the color of the buttons to signal witch table is selected

set(handles.pb_mor_krylov_input,'BackgroundColor',[1;0.843;0]);
set(handles.pb_mor_krylov_output,'BackgroundColor',[0.94;0.94;0.94]);

function pb_mor_krylov_output_Callback(hObject, eventdata, handles)

%Set the right table visible

set(handles.uitable_mor_krylov_output,'Visible','On');
set(handles.uitable_mor_krylov,'Visible','Off');

%Change the color of the buttons to signal witch table is selected

set(handles.pb_mor_krylov_output,'BackgroundColor',[1;0.843;0]);
set(handles.pb_mor_krylov_input,'BackgroundColor',[0.94;0.94;0.94]);



function pb_mor_krylov_infoHermite_Callback(hObject, eventdata, handles)

    infoBox({'InfoHermiteInterpolation.png'});
    uiwait;

function pb_mor_krylov_infoInOut_Callback(hObject, eventdata, handles)

    infoBox({'InfoMomentMatching.png'});
    uiwait;

function pb_mor_krylov_infoAlgoParams_Callback(hObject, eventdata, handles)

    infoBox({'InfoStoppingCriterium.png'});
    uiwait;
    
function pb_mor_krylov_infoExpPoints_Callback(hObject, eventdata, handles)

    if get(handles.pu_mor_krylov_algorithm,'Value')==1  %IRKA
       
        infoBox({'InfoStartingPointsIRKA.png'});
        uiwait;
        
    else
        
        infoBox({'InfoExpansionPointsRK.png'});
        uiwait;
        
    end
    
    
    
function pb_mor_krylov_mimo_Callback(hObject, eventdata, handles)

%Set the struct with the parameters

if get(handles.pu_mor_krylov_algorithm,'Value') == 2
    parameter.algorithm = 'rk';
    if get(handles.rb_mor_krylov_input,'Value') == 1
        parameter.side = 'input';
    elseif get(handles.rb_mor_krylov_output,'Value') == 1
        parameter.side = 'output';
    else
        parameter.side = 'twosided';
    end
else
    parameter.algorithm = 'irka';
    parameter.side = 'twosided';
end

x = get(handles.pu_mor_systems,'String');
y = x{get(handles.pu_mor_systems,'Value')};
    
if ~isempty(y)
    try
        parameter.system = evalin('base',y);              
    catch ex
       errordlg('Selected system could not be found in the workspace.')
       uiwait
       pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles)
       return
    end
else
    errordlg('Please select a valid system fist!','Error Dialog','modal')
    uiwait
    return
end
    
%Read out the data from the tables

tableIn = get(handles.uitable_mor_krylov_MimoExps,'Data');
tableOut = get(handles.uitable_mor_krylov_MimoExps_output,'Data');

%Enable Algorithm- and Side-Selection

set(handles.rb_mor_krylov_input,'Enable','off');
set(handles.rb_mor_krylov_output,'Enable','off');
set(handles.rb_mor_krylov_twosided,'Enable','off');
set(handles.pu_mor_krylov_algorithm,'Enable','off');

%Start the Sub-GUI

data = mimoKrylov({parameter,tableIn,tableOut});

%Write the data from the Sub-GUI to the tables

handles.MimoParam.block = 0;

try
    
    wOut = handles.widthTableMimoKrylovOutput;
    wIn = handles.widthTableMimoKrylovInput;
    
    if ~isempty(data.inputData)
       if size(data.inputData,2) == 2
          tableData = data.inputData;
          tableData = [tableData, cell(size(tableData,1),2)];
          
          %Set default values for the input directions
          
          for i = 1:size(tableData,1)
            tableData{i,3} = vec2string((1:parameter.system.m==min(i,parameter.system.m))); 
          end
          
          %Set default values for the output directions
          
          for i = 1:size(tableData,1)
            tableData{i,4} = vec2string((1:parameter.system.p==min(i,parameter.system.p))); 
          end
          
          set(handles.uitable_mor_krylov_MimoExps,'Data',tableData);
          set(handles.uitable_mor_krylov_MimoExps,'ColumnWidth',{round(wIn/2),round(wIn/2),0,0});
          handles.MimoParam.block = 1;
       elseif size(data.inputData,2) == 3
          tableData = data.inputData;
          tableData = [tableData, cell(size(tableData,1),1)];
          
          %Set default values for the output directions
          
          for i = 1:size(tableData,1)
            tableData{i,4} = vec2string((1:parameter.system.p==min(i,parameter.system.p))); 
          end
          
          set(handles.uitable_mor_krylov_MimoExps,'Data',tableData);
          set(handles.uitable_mor_krylov_MimoExps,'ColumnWidth',{round(wIn/3),round(wIn/3),round(wIn/3),0});
       else
          set(handles.uitable_mor_krylov_MimoExps,'Data',data.inputData);
          set(handles.uitable_mor_krylov_MimoExps_output,'Data',...
              [data.inputData(:,1),data.inputData(:,2),data.inputData(:,4)]);
          set(handles.uitable_mor_krylov_MimoExps,'ColumnWidth',{round(wIn/4),round(wIn/4),round(wIn/4),round(wIn/4)});
       end
    end

    if ~isempty(data.outputData)
       if size(data.outputData,2) == 2
          tableData = data.outputData;
          tableData = [tableData, cell(size(tableData,1),1)];
          
          %Set default values for the output directions
          
          for i = 1:size(tableData,1)
            tableData{i,3} = vec2string((1:parameter.system.p==min(i,parameter.system.p))); 
          end
          
          set(handles.uitable_mor_krylov_MimoExps_output,'Data',tableData); 
          set(handles.uitable_mor_krylov_MimoExps_output,'ColumnWidth',{round(wOut/2),round(wOut/2),0});
          handles.MimoParam.block = 1;
       else
          tableData = data.outputData;
          set(handles.uitable_mor_krylov_MimoExps_output,'Data',tableData);
          set(handles.uitable_mor_krylov_MimoExps_output,'ColumnWidth',{round(wIn/3),round(wIn/3),round(wIn/3)});
       end
    end

    %Set the correct table visible

    if get(handles.rb_mor_krylov_output,'Value') == 1
       set(handles.uitable_mor_krylov_MimoExps,'Visible','off'); 
       set(handles.uitable_mor_krylov_MimoExps_output,'Visible','on'); 
    else
       set(handles.uitable_mor_krylov_MimoExps,'Visible','on'); 
       set(handles.uitable_mor_krylov_MimoExps_output,'Visible','off');  
    end

    if data.hermite || get(handles.rb_mor_krylov_output,'Value') == 1 ...
            || get(handles.rb_mor_krylov_input,'Value') == 1
       set(handles.pb_mor_krylov_MimoExps_Input,'Visible','off'); 
       set(handles.pb_mor_krylov_MimoExps_Output,'Visible','off');
    else
       set(handles.pb_mor_krylov_MimoExps_Input,'Visible','on'); 
       set(handles.pb_mor_krylov_MimoExps_Output,'Visible','on');

       set(handles.pb_mor_krylov_MimoExps_Input,'BackgroundColor',[1;0.843;0]);
       set(handles.pb_mor_krylov_MimoExps_Output,'BackgroundColor',[0.94;0.94;0.94]);
    end
    
    %Set the specified values from the subgui in the handles
    
    handles.MimoParam.hermite = data.hermite;
    
    guidata(hObject,handles);
end

%Enable back all disabled controls

set(handles.rb_mor_krylov_twosided,'Enable','on');
set(handles.pu_mor_krylov_algorithm,'Enable','on');

if get(handles.pu_mor_krylov_algorithm,'Value') ~= 1    %not IRKA
    set(handles.rb_mor_krylov_input,'Enable','on');
    set(handles.rb_mor_krylov_output,'Enable','on');
end

function uitable_mor_krylov_MimoExps_CreateFcn(hObject, eventdata, handles)

    set(hObject,'Data',cell(1,4));

function pb_mor_krylov_MimoExps_Input_Callback(hObject, eventdata, handles)

%Set the right table visible

set(handles.uitable_mor_krylov_MimoExps,'Visible','On');
set(handles.uitable_mor_krylov_MimoExps_output,'Visible','Off');

%Change the color of the buttons to signal witch table is selected

set(handles.pb_mor_krylov_MimoExps_Input,'BackgroundColor',[1;0.843;0]);
set(handles.pb_mor_krylov_MimoExps_Output,'BackgroundColor',[0.94;0.94;0.94]);

function pb_mor_krylov_MimoExps_Output_Callback(hObject, eventdata, handles)

%Set the right table visible

set(handles.uitable_mor_krylov_MimoExps,'Visible','Off');
set(handles.uitable_mor_krylov_MimoExps_output,'Visible','On');

%Change the color of the buttons to signal witch table is selected

set(handles.pb_mor_krylov_MimoExps_Input,'BackgroundColor',[0.94;0.94;0.94]);
set(handles.pb_mor_krylov_MimoExps_Output,'BackgroundColor',[1;0.843;0]);
    
    
%--------------------------------------------------------------------------
%                            SYSTEM ANALYSIS
%--------------------------------------------------------------------------

function pu_an_sys1_Callback(hObject, eventdata, handles)
%System has been selected

    set(handles.figure1,'Pointer','watch')
    
    %Set all displays back to default
    
    if eventdata ~= -1
        set(handles.tx_an_sys1_h2,'String','')
        set(handles.tx_an_sys1_dissipativity,'String','')
        set(handles.tx_an_sys1_decaytime,'String','')
        set(handles.tx_an_sys1_hinf,'String','')
        set(handles.tx_an_sys1_stability,'String','')
        set(handles.tx_an_compare_h2_abs,'String','')
        set(handles.tx_an_compare_h2_rel,'String','')
        set(handles.tx_an_compare_hinf_abs,'String','')
        set(handles.tx_an_compare_hinf_rel,'String','')
        drawnow
    end

    %Get the selected system from workspace
    
    try
        sys = getSysFromWs(hObject);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys1_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
            set(handles.pu_an_sys2,'Value',1);
            set(handles.pu_an_sys2,'String',{''});
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys1_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
        uiwait
        set(handles.pu_an_sys2,'Value',1);
        set(handles.pu_an_sys2,'String',{''});
        return
    end

    displaySystemInformation(handles.tx_an_sys1_sysinfo,sys);
    set(handles.tx_an_sys1_sysinfo, 'HorizontalAlignment','left');
    
    %Set the values which belong to all systems
    
    if eventdata ~= -1 
        
        if isa(sys,'sss') && ~isempty(sys.h2Norm)
            set(handles.tx_an_sys1_h2,'String', num2str(sys.h2Norm))
        else
            set(handles.tx_an_sys1_h2,'String','')
        end

        if isa(sys,'sss') && ~isempty(sys.hInfNorm)
            set(handles.tx_an_sys1_hinf,'String',num2str(sys.hInfNorm))
        else
            set(handles.tx_an_sys1_hinf,'String','')
        end

        if isa(sys,'sss') && ~isempty(sys.decayTime)
            set(handles.tx_an_sys1_decaytime,'String',num2str(sys.decayTime))
        else
            set(handles.tx_an_sys1_decaytime,'String','')
        end
    end
    
    %Update the list of selectable systems for the second system (number of
    %inputs and outputs must be consistent
    
    l=systemsInWorkspace();
    
    for i = 1:length(l)      
        try
           sysTemp = evalin('base',l{i,1}); 
           if size(sysTemp.B,2) ~= size(sys.B,2) || size(sysTemp.C,1) ~= size(sys.C,1)
              l{i,1} = []; 
           end
        catch ex
            l{i,1} = [];
        end        
    end
    
    l = l(~cellfun(@isempty, l));
    
    l_alt = get(handles.pu_an_sys2,'String');
    
    if ~isempty(l_alt)

        sOld=l_alt{get(handles.pu_an_sys2,'Value')};
        indexOld=find(strcmp(sOld,l));

        if ~isempty(indexOld)
            set(handles.pu_an_sys2,'Value',indexOld);
        else
            set(handles.pu_an_sys2,'Value',1);
        end

    else
      set(handles.pu_an_sys2,'Value',1);
    end
    
    set(handles.pu_an_sys2,'String',l);
    pu_an_sys2_Callback(handles.pu_an_sys2,-1,handles);
    
    if isempty(l)
        set(handles.pb_an_compare_h2,'Enable','off');
        set(handles.pb_an_compare_hinf,'Enable','off');
    end
    
    set(handles.figure1,'Pointer','arrow')

function pb_an_sys1_calcall_Callback(hObject, eventdata, handles)
%Calculate everything on the panel for system 1

    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys2);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys2_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys2_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end

    pb_an_sys1_stability_Callback(handles.pb_an_sys1_stability,1,handles);
    pb_an_sys1_dissipativity_Callback(handles.pb_an_sys1_stability,1,handles);
    pb_an_sys1_h2_Callback(handles.pb_an_sys1_h2,1,handles);
    pb_an_sys1_hinf_Callback(handles.pb_an_sys1_hinf,1,handles);
    pb_an_sys1_decaytime_Callback(handles.pb_an_sys1_decaytime,1,handles);

function pb_an_sys1_stability_Callback(hObject, eventdata, handles)
%Determine stability of system
    
    %Get the system from workspace
    
    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys1);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys1_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys1_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow
    
    %Check for Stability
    
    try
        stable = isstable(sys);
    catch ex
        % display the line of the error if "isstable" produced the error
        [found,str] = extractEntryFromErrorStack(ex.stack,'isstable');
        if found            
            errordlg({str;ex.message},'Error Dialog','modal');
        else
            errordlg(ex.message,'Error Dialog','modal')
        end
        uiwait
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
        return;
    end

    %Display the solution to the user
    
    if stable == 1
        set(handles.tx_an_sys1_stability,'String','stable')
        set(handles.tx_an_sys1_stability,'ForegroundColor',[0,0.6,0]);
    elseif stable == 0
        set(handles.tx_an_sys1_stability,'String','unstable')       
        set(handles.tx_an_sys1_stability,'ForegroundColor',[0.8,0,0]);
    else
        set(handles.tx_an_sys1_stability,'String','unknown')
        set(handles.tx_an_sys1_stability,'ForegroundColor',[0,0,0]);
    end

    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')

function pb_an_sys1_dissipativity_Callback(hObject, eventdata, handles)
%Determin if the system is strictly dissipative or not    
    
    %Get the system from workspace
    
    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys1);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys1_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys1_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow
    
    %Check for strict dissipativity

    try
        dissipativ = issd(sys);
    catch ex
        % display the line of the error if "issd" produced the error
        [found,str] = extractEntryFromErrorStack(ex.stack,'issd');
        if found            
            errordlg({str;ex.message},'Error Dialog','modal');
        else
            errordlg(ex.message,'Error Dialog','modal')
        end
        uiwait
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
        return;
    end

    %Display the solution to the user
    
    if dissipativ == 1
        set(handles.tx_an_sys1_dissipativity,'String','strictly dissipative')
    elseif dissipativ == 0
        set(handles.tx_an_sys1_dissipativity,'String','not strictly dissipative')
    else
        set(handles.tx_an_sys1_dissipativity,'String','unknown')
    end

    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')

function pb_an_sys1_h2_Callback(hObject, eventdata, handles)
%Determine h2-norm of system

    %Get system from workspace

    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys1);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys1_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys1_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow

    %Get the h2-norm of the system

    if isa(sys,'sss') && ~isempty(sys.h2Norm)

        h2 = sys.h2Norm;
        
    else
        
        try
            h2 = norm(sys, 2);
        catch ex
            if strcmp(ex.identifier,'MATLAB:nomem')
                errordlg('Out of memory, system is too large to solve lyapunov equotation','Error Dialog','modal')
            elseif strcmp(ex.identifier,'Control:foundation:LyapChol4')
                errordlg('A or (A,E) must have all their eigenvalues in the left-half plane','Error Dialog','modal')
            else
                % display the line of the error if "norm" produced the error
                [found,str] = extractEntryFromErrorStack(ex.stack,'norm');
                if found            
                    errordlg({str;ex.message},'Error Dialog','modal');
                else
                    errordlg(ex.message,'Error Dialog','modal')
                end
            end
            uiwait
            set(handles.figure1,'Pointer','arrow')
            set(handles.virtgr_an_red_buttons,'Enable','on')
            return
        end
    end

    %Display the solution to the user

    set(handles.tx_an_sys1_h2,'String',sprintf('%0.4g',h2))

    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')

function pb_an_sys1_hinf_Callback(hObject, eventdata, handles)
%Determine h-infinity-norm of system
    
    %Get system from workspace
    
    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys1);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys1_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys1_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow
    
    %Get hInf-Norm

    if isa(sys,'sss') && ~isempty(sys.hInfNorm)   
        hinf = sys.hInfNorm;
    else
        
        try
            hinf=norm(sys, inf);
        catch ex
            % display the line of the error if "norm" produced the error
            [found,str] = extractEntryFromErrorStack(ex.stack,'norm');
            if found            
                errordlg({str;ex.message},'Error Dialog','modal');
            else
                errordlg(ex.message,'Error Dialog','modal')
            end
            set(handles.figure1,'Pointer','arrow')
            set(handles.virtgr_an_red_buttons,'Enable','on')
            uiwait
            return
        end
    end

    %Display solution to the user
    
    set(handles.tx_an_sys1_hinf,'String',sprintf('%0.4g',hinf))

    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
    
function pb_an_sys1_decaytime_Callback(hObject, eventdata, handles)
%Determin the decay time for the system

    %Get the system from workspace
    
    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys1);     
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys1_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys1_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys1_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow

    %Get the decay Time
    
    if ~isfield(sys,'decayTime') || isempty(sys.decayTime)
        try
            decTime = decayTime(sys);
        catch ex
            % display the line of the error if "decayTime" produced the error
            [found,str] = extractEntryFromErrorStack(ex.stack,'decayTime');
            if found            
                errordlg({str;ex.message},'Error Dialog','modal');
            else
                errordlg(ex.message,'Error Dialog','modal')
            end
            set(handles.figure1,'Pointer','arrow')
            set(handles.virtgr_an_red_buttons,'Enable','on')
            uiwait
            return
        end
    else
        decTime = sys.decayTime;
    end
    
    %Display solution to the user
    
    set(handles.tx_an_sys1_decaytime,'String',decTime)
    
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')

    
    
function pu_an_sys2_Callback(hObject, eventdata, handles)
%System has been selected
    
    %Set all displays back to default
    
    if eventdata ~= -1
        set(handles.tx_an_sys2_h2,'String','')
        set(handles.tx_an_sys2_dissipativity,'String','')
        set(handles.tx_an_sys2_decaytime,'String','')
        set(handles.tx_an_sys2_hinf,'String','')
        set(handles.tx_an_sys2_stability,'String','')
        set(handles.tx_an_compare_h2_abs,'String','')
        set(handles.tx_an_compare_h2_rel,'String','')
        set(handles.tx_an_compare_hinf_abs,'String','')
        set(handles.tx_an_compare_hinf_rel,'String','')
        drawnow
    end

    %Get the selected system from workspace
    
    try
        sys = getSysFromWs(hObject);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys2_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys2_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
        uiwait
        return
    end

    displaySystemInformation(handles.tx_an_sys2_sysinfo,sys);
    set(handles.tx_an_sys2_sysinfo, 'HorizontalAlignment','left');
    
    %Set the values which belong to all systems
    
    if eventdata ~= -1
        
        if isa(sys,'sss') && ~isempty(sys.h2Norm)
            set(handles.tx_an_sys2_h2,'String', num2str(sys.h2Norm))
        else
            set(handles.tx_an_sys2_h2,'String','')
        end

        if isa(sys,'sss') && ~isempty(sys.hInfNorm)
            set(handles.tx_an_sys2_hinf,'String',num2str(sys.hInfNorm))
        else
            set(handles.tx_an_sys2_hinf,'String','')
        end

        if isa(sys,'sss') && ~isempty(sys.decayTime)
            set(handles.tx_an_sys2_decaytime,'String',num2str(sys.decayTime))
        else
            set(handles.tx_an_sys2_decaytime,'String','')
        end
    
    end
    
function pb_an_sys2_calcall_Callback(hObject, eventdata, handles)
%Calculate everything on the panel for system 2
    
    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys2);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys2_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys2_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end

    pb_an_sys2_stability_Callback(handles.pb_an_sys2_stability,1,handles);
    pb_an_sys2_dissipativity_Callback(handles.pb_an_sys2_stability,1,handles);
    pb_an_sys2_h2_Callback(handles.pb_an_sys2_h2,1,handles);
    pb_an_sys2_hinf_Callback(handles.pb_an_sys2_hinf,1,handles);
    pb_an_sys2_decaytime_Callback(handles.pb_an_sys2_decaytime,1,handles);
    
function pb_an_sys2_stability_Callback(hObject, eventdata, handles)
%Determine stability of system
    
    %Get the system from workspace
    
    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys2);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys2_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys2_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow
    
    %Check for Stability
    
    try
        stable = isstable(sys);
    catch ex
        % display the line of the error if "isstable" produced the error
        [found,str] = extractEntryFromErrorStack(ex.stack,'isstable');
        if found            
            errordlg({str;ex.message},'Error Dialog','modal');
        else
            errordlg(ex.message,'Error Dialog','modal')
        end
        uiwait
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
        return;
    end

    %Display the solution to the user
    
    if stable == 1
        set(handles.tx_an_sys2_stability,'String','stable')
        set(handles.tx_an_sys2_stability,'ForegroundColor',[0,0.6,0]);
    elseif stable == 0
        set(handles.tx_an_sys2_stability,'String','unstable')       
        set(handles.tx_an_sys2_stability,'ForegroundColor',[0.8,0,0]);
    else
        set(handles.tx_an_sys2_stability,'String','unknown')
        set(handles.tx_an_sys2_stability,'ForegroundColor',[0,0,0]);
    end

    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
    
function pb_an_sys2_dissipativity_Callback(hObject, eventdata, handles)
%Determin if the system is strictly dissipative or not    
    
    %Get the system from workspace
    
    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys2);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys2_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys2_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow
    
    %Check for strict dissipativity

    try
        dissipativ = issd(sys);
    catch ex
        % display the line of the error if "issd" produced the error
        [found,str] = extractEntryFromErrorStack(ex.stack,'issd');
        if found            
            errordlg({str;ex.message},'Error Dialog','modal');
        else
            errordlg(ex.message,'Error Dialog','modal')
        end
        uiwait
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
        return;
    end

    %Display the solution to the user
    
    if dissipativ == 1
        set(handles.tx_an_sys2_dissipativity,'String','strictly dissipative')
    elseif dissipativ == 0
        set(handles.tx_an_sys2_dissipativity,'String','not strictly dissipative')
    else
        set(handles.tx_an_sys2_dissipativity,'String','unknown')
    end

    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
    
function pb_an_sys2_h2_Callback(hObject, eventdata, handles)
%Determine h2-norm of system

    %Get system from workspace

    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys2);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys2_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys2_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow

    %Get the h2-norm of the system

    if isa(sys,'sss') && ~isempty(sys.h2Norm)
        h2 = sys.h2Norm;
    else
        try
            h2 = norm(sys, 2);
        catch ex
            if strcmp(ex.identifier,'MATLAB:nomem')
                errordlg('Out of memory, system is too large to solve lyapunov equotation','Error Dialog','modal')
            elseif strcmp(ex.identifier,'Control:foundation:LyapChol4')
                errordlg('A or (A,E) must have all their eigenvalues in the left-half plane','Error Dialog','modal')
            else
                % display the line of the error if "norm" produced the error
                [found,str] = extractEntryFromErrorStack(ex.stack,'norm');
                if found            
                    errordlg({str;ex.message},'Error Dialog','modal');
                else
                    errordlg(ex.message,'Error Dialog','modal')
                end
            end
            uiwait
            set(handles.figure1,'Pointer','arrow')
            set(handles.virtgr_an_red_buttons,'Enable','on')
            return
        end
    end

    %Display the solution to the user

    set(handles.tx_an_sys2_h2,'String',sprintf('%0.4g',h2))

    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
    
function pb_an_sys2_hinf_Callback(hObject, eventdata, handles)
%Determine h-infinity-norm of system
    
    %Get system from workspace
    
    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys2);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys2_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys2_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow
    
    %Get hInf-Norm

    if isa(sys,'sss') && ~isempty(sys.hInfNorm)
        hinf = sys.hInfNorm;
    else
        try
            hinf=norm(sys, inf);
        catch ex
            % display the line of the error if "norm" produced the error
            [found,str] = extractEntryFromErrorStack(ex.stack,'norm');
            if found            
                errordlg({str;ex.message},'Error Dialog','modal');
            else
                errordlg(ex.message,'Error Dialog','modal')
            end
            uiwait
            set(handles.figure1,'Pointer','arrow')
            set(handles.virtgr_an_red_buttons,'Enable','on')
            return
        end
    end

    %Display solution to the user
    
    set(handles.tx_an_sys2_hinf,'String',sprintf('%0.4g',hinf))

    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
    
function pb_an_sys2_decaytime_Callback(hObject, eventdata, handles)
%Determin the decay time for the system

    %Get the system from workspace
    
    try
        [sys,sysname] = getSysFromWs(handles.pu_an_sys2);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            set(handles.tx_an_sys2_sysinfo,'String','Please select a system!')
            set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
            errordlg('Please select a system first','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        set(handles.tx_an_sys2_sysinfo,'String','Invalid model')
        set(handles.tx_an_sys2_sysinfo,'HorizontalAlignment','center')
        uiwait
        pb_refreshsys_Callback(handles.pb_refreshsys,eventdata,handles);
        return
    end
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow

    %Get the decay Time
    
    if ~isfield(sys,'decayTime') || isempty(sys.decayTime)
        try
            decTime = decayTime(sys);
        catch ex
            % display the line of the error if "decayTime" produced the error
            [found,str] = extractEntryFromErrorStack(ex.stack,'decayTime');
            if found            
                errordlg({str;ex.message},'Error Dialog','modal');
            else
                errordlg(ex.message,'Error Dialog','modal')
            end
            uiwait
            set(handles.figure1,'Pointer','arrow')
            set(handles.virtgr_an_red_buttons,'Enable','on')
            return
        end
    else
        decTime = sys.decayTime;
    end
    
    %Display solution to the user
    
    set(handles.tx_an_sys2_decaytime,'String',decTime)
    
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')  
    
    
    
function pb_an_compare_h2_Callback(hObject, eventdata, handles)
%Calculate the h2-norm of the difference-system of the two systems

    %Get the h2-norm of the first system

    norm1 = get(handles.tx_an_sys1_h2,'String');

    if isempty(norm1)
       pb_an_sys1_h2_Callback(handles.pb_an_sys1_h2,-1,handles);
       norm1 = get(handles.tx_an_sys1_h2,'String'); 
    end

    norm1 = str2double(norm1);
    
    %Get the selected system from workspace
    
    try
        sys1 = getSysFromWs(handles.pu_an_sys1);
        sys2 = getSysFromWs(handles.pu_an_sys2);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            errordlg('No system selected!','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        uiwait
        return
    end

    %Calculate the h2-norm of the diffence-system
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow
    
    try
       normAbs = norm(sys1-sys2); 
    catch ex
        % display the line of the error if "norm" produced the error
        [found,str] = extractEntryFromErrorStack(ex.stack,'norm');
        if found            
            errordlg({str;ex.message},'Error Dialog','modal');
        else
            errordlg(ex.message,'Error Dialog','modal')
        end
        uiwait
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
        return
    end
    
    %Display the solution to the user
    
    set(handles.tx_an_compare_h2_abs,'String',sprintf('%0.4g',normAbs));
    set(handles.tx_an_compare_h2_rel,'String',sprintf('%0.4g',normAbs/norm1));
    
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')

function pb_an_compare_hinf_Callback(hObject, eventdata, handles)
%Calculate the h-inf-norm of the difference-system of the two systems

    %Get the h-inf-norm of the first system

    norm1 = get(handles.tx_an_sys1_hinf,'String');

    if isempty(norm1)
       pb_an_sys1_hinf_Callback(handles.pb_an_sys1_hinf,-1,handles);
       norm1 = get(handles.tx_an_sys1_hinf,'String'); 
    end

    norm1 = str2double(norm1);
    
    %Get the selected system from workspace
    
    try
        sys1 = getSysFromWs(handles.pu_an_sys1);
        sys2 = getSysFromWs(handles.pu_an_sys2);
    catch ex
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier, 'unassigned')
            errordlg('No system selected!','Error Dialog','modal')
            uiwait
            return
        end
        errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
        uiwait
        return
    end

    %Calculate the h2-norm of the diffence-system
    
    set(handles.figure1,'Pointer','watch')
    set(handles.virtgr_an_red_buttons,'Enable','off')
    drawnow
    
    try
       normAbs = norm(sys1-sys2,inf); 
    catch ex
        % display the line of the error if "norm" produced the error
        [found,str] = extractEntryFromErrorStack(ex.stack,'norm');
        if found            
            errordlg({str;ex.message},'Error Dialog','modal');
        else
            errordlg(ex.message,'Error Dialog','modal')
        end
        uiwait
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
        return
    end
    
    %Display the solution to the user
    
    set(handles.tx_an_compare_hinf_abs,'String',sprintf('%0.4g',normAbs));
    set(handles.tx_an_compare_hinf_rel,'String',sprintf('%0.4g',normAbs/norm1));
    
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on') 

function pb_an_compare_h2_info_Callback(hObject, eventdata, handles)

    infoBox({'InfoDifferenceSystemH2Norm.png'});
    uiwait;

function pb_an_compare_hinf_info_Callback(hObject, eventdata, handles)

    infoBox({'InfoDifferenceSystemHinfNorm.png'});
    uiwait;

    
    
%--------------------------------------------------------------------------
%                               FOOTER
%--------------------------------------------------------------------------

%Callbacks of the panel

function logos_footer_ButtonDownFcn(hObject, eventdata, handles)
    %Web-links to the diverent homepages (Open if the user klicks on a logo)
    %(to convert from pixel units to chararcter units for Windows OS the 
    %values are multiplied with 0.2 for horizontal dimensions and with 
    %0.0769 for vertical dimensions)
    
    p=get(hObject,'CurrentPoint');
    p(1,1) = p(1,1)*handles.PixToCharWidth/0.2;
    p(1,2) = p(1,2)*handles.PixToCharHeight/0.0769;

    if p(1,1)>990 && p(1,1)<1040
        web www.tum.de
    elseif p(1,1)>5 && p(1,1)<80
        web www.rt.mw.tum.de/forschung/forschungsgebiete/modellreduktion/sssmor/
    elseif p(1,1)>650 && p(1,1)<765 
        web www.rt.mw.tum.de
    elseif p(1,1)>305 && p(1,1)<420
        web www.rt.mw.tum.de/forschung/forschungsgebiete/modellreduktion/
    end
    
    
    
%--------------------------------------------------------------------------
%                           UTILITY FUNCTIONS
%--------------------------------------------------------------------------


%Read out the expension-points from the GUI

function [s0,sOut]=getExpansionPoints(handles)
% read expansion point(s) from GUI

    %Set Default values
    
    s0 = [];
    sOut = [];

    %Read the expension points from the first table

    if get(handles.uitable_mor_krylov,'UserData')==1
        error('Please correct expansion points first');
    end
    cellEntries = get(handles.uitable_mor_krylov,'Data');
    if ~checkPointsSisoKrylov(cellEntries)
       error('Table with the shifts and moments has the wrong format'); 
    end
    s0=cell2mat(cellEntries);% Entwicklungspunkte
    if isempty(s0) || size(s0,2)==1 % Tabelle leer, oder nur Momente ohne Punkt
        error('Please choose expansion points and number of matching moments!');
    end
    x = get(handles.pu_mor_systems,'String');
    y=x{get(handles.pu_mor_systems,'Value')};
    s=evalin('base',sprintf('size(%s.A,1);',y)); %Größe Originalsystem
    if s<sum(s0(:,2))
       error('Number of matching moments must be smaller than dimension of original system!');
    end 
    s0=transpose(s0);
    
    %Check if input, output or two-sided is selected

    if get(handles.rb_mor_krylov_input,'Value') == 1    %Input
        
        sIn = s0;
        sOut = [];     
        
    elseif get(handles.rb_mor_krylov_output,'Value') == 1   %Output
        
        sOut = s0;
        sIn = [];
        
    else        %Two sided
        
        %Check whethter hermite-interpolation or not
        
        if get(handles.cb_mor_krylov_hermite,'Value') == 1
           
            sIn = s0;
            sOut = s0;
            
        else
            
            %Read out the expension points from the output table
            
            sIn = s0;
            
            if get(handles.uitable_mor_krylov_output,'UserData')==1
                error('Please correct expansion points first');
            end
            sOut=cell2mat(get(handles.uitable_mor_krylov_output,'Data'));% Entwicklungspunkte
            if isempty(sOut) || size(sOut,2)==1 % Tabelle leer, oder nur Momente ohne Punkt
                error('Please choose expansion points and number of matching moments');
            end
            if sum(sIn(2,:)) ~= sum(sOut(:,2))
               error('Number of matching moments must be equal for input and output');
            end 
            sOut=transpose(sOut);
            
        end       
    end
    
function [sIn, sOut, Rt, Lt] = getMimoExpensionPoints(handles)
%Reads out the expansion points for Mimo systems from the GUI

    tableIn = get(handles.uitable_mor_krylov_MimoExps,'Data');
    tableOut = get(handles.uitable_mor_krylov_MimoExps_output,'Data');
    
    sIn = cell2mat(tableIn(:,1:2))';
    sOut = cell2mat(tableOut(:,1:2))';
    
    %Read out the tangential directions if block-krylov is not selected
    
    if handles.MimoParam.block == 0
    
        for i = 1:size(tableIn,1)
           tableIn{i,3} = evalin('base',tableIn{i,3}); 
        end

        for i = 1:size(tableOut,1)
           tableOut{i,3} = evalin('base',tableOut{i,3});  
        end

        Rt = zeros(size(tableIn{1,3},1),size(tableIn,1));
        Lt = zeros(size(tableOut{1,3},1),size(tableOut,1));

        for i = 1:size(Rt,2)
            for j = 1:size(Rt,1)
                Rt(j,i) = tableIn{i,3}(j,1);
            end
        end

        %Read out the output directions from the forth column if twosided and
        %hermite-interpolation is selected or algorithm is IRKA

        if (get(handles.rb_mor_krylov_twosided,'Value') == 1 && handles.MimoParam.hermite == 1) ||...
                get(handles.pu_mor_krylov_algorithm,'Value') == 1

            for i = 1:size(tableIn,1)
                tableIn{i,4} = evalin('base',tableIn{i,4}); 
            end

            for i = 1:size(Lt,2)
                for j = 1:size(Lt,1)
                    Lt(j,i) = tableIn{i,4}(j,1);
                end
            end

            sOut = sIn;
            
        else

            for i = 1:size(Lt,2)
                for j = 1:size(Lt,1)
                    Lt(j,i) = tableOut{i,3}(j,1);
                end
            end

        end

        Rt = getDirectionMatrix(sIn,Rt);
        Lt = getDirectionMatrix(sOut,Lt);
    
    else
       
        Rt = [];
        Lt = [];
        
    end
    
    if get(handles.rb_mor_krylov_output,'Value') == 1       %Output
        
        sIn = [];
        Rt = [];     
        
    elseif get(handles.rb_mor_krylov_input,'Value') == 1    %Input
        
        sOut = [];
        Lt = [];
        
    else                                                    %Twosided
        
        if handles.MimoParam.hermite == 1       %Hermite Interpolation
            
            sOut = sIn;
        end      
    end
       
function m = getDirectionMatrix(s0,mOld)
%Gives back the matrix with the tangential directions in a format which can
%be used for the reduction algorithms (multiple shifts are taken into account)

    if ~isempty(s0)
       
        numberOfMoments = sum(s0(2,:));
        
        m = zeros(size(mOld,1),numberOfMoments);
        
        counter = 1;
        
        for i = 1:size(s0,2)            
           for j = 1:s0(2,i)              
               for k = 1:size(mOld,1)
                  
                   m(k,counter) = mOld(k,i);
                        
               end
               
               counter = counter + 1;
               
           end            
        end
        
    else
        
        m = [];
        
    end
        
        
%Functions for testing Plot-attributes in the Visualisation-Menue-Point     

function testWidth(hObject)
% check inserted linewidth

h=str2num(get(hObject,'String')); %#ok<ST2NM>
if isempty(h)|| isnan(h)
    errordlg('Line Width must be a real number grater than zero','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
elseif imag(h)~=0
    errordlg('No imaginary numbers allowed','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
elseif length(h)>1
    if length(h)==2
        h=sprintf('%i.%i',round(h(1)),round(h(2)));
        set(hObject,'String',h)
        set(hObject,'UserData',0)
        return
    end
    errordlg('Use ''.'' as decimal seperator','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
elseif h<0
    errordlg('Line Width must not be less than zero','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
else
    set(hObject,'UserData',0)
end

function testMax(hObject,handles)
% check maximum frequency/time in manual mode
h=str2num(get(hObject,'String')); %#ok<ST2NM>
if isempty(h) || isnan(h)
    errordlg('Max must be greater than Min. Use ''.'' as decimal seperator. Don''t use characters','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
    set(hObject,'String','')
elseif imag(h)~=0
    errordlg('No imaginary numbers allowed','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
    set(hObject,'String','')
elseif length(h)>1
    if length(h)==2
        h=sprintf('%i.%i',round(h(1)),round(h(2)));
        set(hObject,'String',h)
        set(hObject,'UserData',0)
        return
    end
    errordlg('Use ''.'' as decimal seperator','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
    set(hObject,'String','')
elseif h<0
    errordlg('Max must not be less than zero','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
    set(hObject,'String','')
elseif h<= str2num(get(handles.ed_min,'String')) %#ok<ST2NM>
    errordlg('Max must be greater than Min','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
    set(hObject,'String','')
else
    set(hObject,'UserData',0)
end

function testMin(hObject,handles)
% check minimum frequency/time in manual mode
h=str2num(get(hObject,'String')); %#ok<ST2NM>
if isempty(h) || isnan(h)
    errordlg('The minimal frequency must be greater than zero. Use ''.'' as decimal seperator. Don''t use characters','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
    set(hObject,'String','')
elseif imag(h)~=0
    errordlg('No imaginary numbers allowed','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
    set(hObject,'String','')
elseif length(h)>1
    if length(h)==2
        h=sprintf('%i.%i',round(h(1)),round(h(2)));
        set(hObject,'String',h)
        set(hObject,'UserData',0)
        return
    end
    errordlg('Use ''.'' as decimal seperator','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
    set(hObject,'String','')
elseif h<0
    errordlg('The minimal frequency must not be less than zero','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
    set(hObject,'String','')
else
    set(hObject,'UserData',0)
end


%Import a matrix with expansion-points from the workspace

function x=listS0InWorkspace
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
        
        if correctExpensionPointMatrix(vec) == 1
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

function t=correctExpensionPointMatrix(m)
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

    
%Auxiliary-functions for krylov
    
function data = addMatrixToTable(m)
%Takes a matrix m with expension points and creates a cell array that can
%be used as Data for an uitable out of the values of m
    
    if size(m,1) == 1 || size(m,2) == 1
        
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
    if size(m,2)==2 && size(m,1)==2
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
    
    data=num2cell(m);   
    
function [] = setKrylovDefaultValues(handles)
       
if get(handles.pu_mor_method,'Value')==3        %Krylov selected
   
    x = get(handles.pu_mor_systems,'String');
    y = x{get(handles.pu_mor_systems,'Value')};
    
    if ~isempty(y)
        
        sys = evalin('base', y);        
    
        if sys.m > 1 || sys.p > 1                       %Mimo system

            valueInput = get(handles.uitable_mor_krylov_MimoExps,'Data');
            valueOutput = get(handles.uitable_mor_krylov_MimoExps_output,'Data');
            
            %Check if the directions in the table are consistent with the
            %number of inputs and outputs of the current system
            
            consistent = 0;
            
            if isempty(cell2mat(valueInput(1,3)))
               if ~isempty(cell2mat(valueInput(1,4)))
                  try
                     vecTemp = evalin('base',valueInput{1,4});
                     if length(vecTemp) == sys.p
                        consistent = 1; 
                     end
                  end
               end                
            elseif isempty(cell2mat(valueInput(1,4)))
                if ~isempty(cell2mat(valueInput(1,3)))
                  try
                     vecTemp = evalin('base',valueInput{1,3});
                     if length(vecTemp) == sys.m
                        consistent = 1; 
                     end
                  end
                end
            else
                try
                    vecTemp = evalin('base',valueInput{1,3});
                    vecTemp2 = evalin('base',valueInput{1,4});
                    if length(vecTemp) == sys.m && length(vecTemp2) == sys.p
                        consistent = 1; 
                    end
                end
            end
            
            %Set the default parameters for MOR with MIMO systems if the
            %selected system has changed
            
            if ~consistent
                handles.MimoParam.hermite = 1;
                handles.MimoParam.block = 0;

                guidata(handles.pu_mor_systems,handles);
            end
            
            %Set the default Values if no values are specified or the
            %system has been changed
            
            if (isempty(cell2mat(valueInput(1,1))) && isempty(cell2mat(valueInput(1,2))) ...
                  && isempty(cell2mat(valueInput(1,3))) && isempty(cell2mat(valueInput(1,4))) ...
                  && isempty(cell2mat(valueOutput(1,1))) && isempty(cell2mat(valueOutput(1,2))) ...
                  && isempty(cell2mat(valueOutput(1,3))) && size(valueOutput,1) == 1 ...
                  && size(valueInput,1) == 1) || ~consistent
                  
                    
                if size(sys.A,1) > 6

                   %Shifts and Moments 

                   value{1,1} = 0;
                   value{1,2} = 2;
                   value{2,1} = 1+i;
                   value{2,2} = 2;
                   value{3,1} = 1-i;
                   value{3,2} = 2;

                   valueInput = value;
                   valueOutput = value;

                   %Tangential Directions

                   valueInput{1,3} = vec2string((1:sys.m==1));
                   valueInput{2,3} = vec2string((1:sys.m==min(2,sys.m)));
                   valueInput{3,3} = vec2string((1:sys.m==min(3,sys.m)));

                   valueInput{1,4} = vec2string((1:sys.p==1));
                   valueInput{2,4} = vec2string((1:sys.p==min(2,sys.p)));
                   valueInput{3,4} = vec2string((1:sys.p==min(3,sys.p)));

                   valueOutput{1,3} = vec2string((1:sys.p==1));
                   valueOutput{2,3} = vec2string((1:sys.p==min(2,sys.p)));
                   valueOutput{3,3} = vec2string((1:sys.p==min(3,sys.p)));

                else

                   %Shifts and Moments

                   value{1,1} = 0;
                   value{1,2} = 1;

                   valueInput = value;
                   valueOutput = value;

                   %Tangential Directions

                   valueInput{1,3} = vec2string((1:sys.m==1));
                   valueInput{1,4} = vec2string((1:sys.p==1));
                   valueOutput{1,3} = vec2string((1:sys.p==1));

                end
            end
            
            set(handles.uitable_mor_krylov_MimoExps,'Data',valueInput);
            set(handles.uitable_mor_krylov_MimoExps_output,'Data',valueOutput);

        else                                        %Siso system

            value = get(handles.uitable_mor_krylov,'Data');

            if isempty(cell2mat(value(1,1))) && isempty(cell2mat(value(1,2)))

                   if size(sys.A,1) > 6
                       value{1,1} = 0;
                       value{1,2} = 2;
                       value{2,1} = 1+i;
                       value{2,2} = 2;
                       value{3,1} = 1-i;
                       value{3,2} = 2;           
                   else
                       value{1,1} = 0;
                       value{1,2} = 1;
                   end
                   set(handles.uitable_mor_krylov,'Data',value);
                   set(handles.uitable_mor_krylov_output,'Data',value);
                   countMatchedMoments(handles);       
            end  
        end
    end
end  
           
function [string] = vec2string(vect)
        
    string = '[';
    for i = 1:length(vect)-1
        string = strcat(strcat(string,num2str(vect(i))),';');
    end
    string = strcat(strcat(string,num2str(vect(length(vect)))),']');
   
function [] = countMatchedMoments(handles)
%Updates the number of matched moemts if for MOR with Krylov

    try
        %Read out the data from the tables for input and output shifts
        
        x=get(handles.uitable_mor_krylov,'Data');
        data1=cell2mat(x(:,2));
        x = get(handles.uitable_mor_krylov_output,'Data');
        data2=cell2mat(x(:,2));
        
        if ~isempty(data1)                  %Input
           data1 = sum(data1(:,1)); 
        else
           data1 = 0;
        end

        if ~isempty(data2)                  %Output
           data2 = sum(data2(:,1)); 
        else
           data2 = 0;
        end      
        
        if isnan(data1)
           data1 = 0; 
        end
        
        if isnan(data2)
           data2 = 0; 
        end
        
        %Update the displays of reduced order and matched moments
        
        redOrder = 0;

        if get(handles.rb_mor_krylov_twosided,'Value') == 1     %Two sided
            if get(handles.pu_mor_krylov_algorithm,'Value') == 2 && ...
                    get(handles.cb_mor_krylov_hermite,'Value') == 0  %Input and Output specifiable

                set(handles.st_mor_krylov_matchedMom,'String',data2+data1);
            else
                set(handles.st_mor_krylov_matchedMom,'String',2*data1);
            end
            
            redOrder = data1;
            
        elseif get(handles.rb_mor_krylov_input,'Value') == 1    %Input
            set(handles.st_mor_krylov_matchedMom,'String',data1);
            redOrder = data1;
        else                                             %Output
            set(handles.st_mor_krylov_matchedMom,'String',data2);
            redOrder = data2;
        end
        
        set(handles.st_mor_krylov_redOrder,'String',redOrder);

        %Check weather the reduced order is bigger than the original order

        x = get(handles.pu_mor_systems,'String');
        y = x{get(handles.pu_mor_systems,'Value')};

        if ~isempty(y)
           sys = evalin('base',y);
           if sys.n < redOrder
              msgbox('Reduced order is bigger than the original order of the system. Please correct that before reducing the system.','Warning','Warn');
              uiwait;
           end
        end
    end
    
function result = checkPointsSisoKrylov(cellArray)
%Cheks if the values for the points and moments specified by the user are
%valid and gives back the corrected Table if possible

    result = 1;

    for i = 1:size(cellArray,1)
      if isempty(cellArray{i,1}) || isnan(cellArray{i,1}) || ...
         isempty(cellArray{i,2}) || isnan(cellArray{i,2}) || ...
         mod(cellArray{i,2},1) ~= 0 
        
         result = 0;
         return;            
      end     
    end
        
           
    
%Different Layouts for SISO- and MIMO-Krylov
    
function [] = layoutMimoKrylov(handles)
%Sets the right table with shifts and directions visible for Mimo-systems 

    wOut = handles.widthTableMimoKrylovOutput;
    wIn = handles.widthTableMimoKrylovInput;

    if get(handles.pu_mor_krylov_algorithm,'Value') == 1        %IRKA
        
        set(handles.uitable_mor_krylov_MimoExps,'Visible','on');
        set(handles.uitable_mor_krylov_MimoExps_output,'Visible','off');
        
        set(handles.pb_mor_krylov_MimoExps_Input,'Visible','off');
        set(handles.pb_mor_krylov_MimoExps_Output,'Visible','off');
        
        set(handles.uitable_mor_krylov_MimoExps,'ColumnWidth',{round(wIn/4) round(wIn/4) round(wIn/4) round(wIn/4)});
        
    else                                                        %RK
        
        if get(handles.rb_mor_krylov_input,'Value') == 1
            
            set(handles.uitable_mor_krylov_MimoExps,'Visible','on');
            set(handles.uitable_mor_krylov_MimoExps_output,'Visible','off');
            
            set(handles.pb_mor_krylov_MimoExps_Input,'Visible','off');
            set(handles.pb_mor_krylov_MimoExps_Output,'Visible','off');
        
            set(handles.uitable_mor_krylov_MimoExps,'ColumnWidth',{round(wIn/3) round(wIn/3) round(wIn/3) 0});
            
        elseif get(handles.rb_mor_krylov_output,'Value') == 1
            
            set(handles.uitable_mor_krylov_MimoExps,'Visible','off');
            set(handles.uitable_mor_krylov_MimoExps_output,'Visible','on');
            
            set(handles.pb_mor_krylov_MimoExps_Input,'Visible','off');
            set(handles.pb_mor_krylov_MimoExps_Output,'Visible','off');           
            
        else                                                    %twosided
            
            if handles.MimoParam.hermite == 1
               
                set(handles.uitable_mor_krylov_MimoExps,'Visible','on');
                set(handles.uitable_mor_krylov_MimoExps_output,'Visible','off');
        
                set(handles.pb_mor_krylov_MimoExps_Input,'Visible','off');
                set(handles.pb_mor_krylov_MimoExps_Output,'Visible','off');
        
                set(handles.uitable_mor_krylov_MimoExps,'ColumnWidth',{round(wIn/4) round(wIn/4) round(wIn/4) round(wIn/4)});              
                
            else
                
                set(handles.uitable_mor_krylov_MimoExps,'Visible','on');
                set(handles.uitable_mor_krylov_MimoExps_output,'Visible','off');
        
                set(handles.pb_mor_krylov_MimoExps_Input,'Visible','on');
                set(handles.pb_mor_krylov_MimoExps_Output,'Visible','on');
                
                set(handles.pb_mor_krylov_MimoExps_Input,'BackgroundColor',[1 0.843 0]);
                set(handles.pb_mor_krylov_MimoExps_Output,'BackgroundColor',[0.94 0.94 0.94]);
        
                set(handles.uitable_mor_krylov_MimoExps,'ColumnWidth',{round(wIn/3) round(wIn/3) round(wIn/3) 0});
                
            end          
        end   
    end    
    
function [] = layoutSisoKrylov(handles)

    if get(handles.pu_mor_krylov_algorithm,'Value')~=1 %RK
        
        if get(handles.rb_mor_krylov_input,'Value') == 1
        
            set(handles.uitable_mor_krylov_output,'Visible','off');
            set(handles.uitable_mor_krylov,'Visible','on');
            
            set(handles.pb_mor_krylov_input,'Visible','off');
            set(handles.pb_mor_krylov_output,'Visible','off');
        
            set(handles.cb_mor_krylov_hermite,'Enable','off');
            
        elseif get(handles.rb_mor_krylov_output,'Value') == 1
            
            set(handles.uitable_mor_krylov_output,'Visible','on');
            set(handles.uitable_mor_krylov,'Visible','off');
            
            set(handles.pb_mor_krylov_input,'Visible','off');
            set(handles.pb_mor_krylov_output,'Visible','off');
        
            set(handles.cb_mor_krylov_hermite,'Enable','off')
            
        else        %twosided
        
            set(handles.uitable_mor_krylov_output,'Visible','off');
            set(handles.uitable_mor_krylov,'Visible','on');
            
            if get(handles.cb_mor_krylov_hermite,'Value') == 0
                
                set(handles.pb_mor_krylov_input,'Visible','on');
                set(handles.pb_mor_krylov_output,'Visible','on');

                set(handles.pb_mor_krylov_input,'BackgroundColor',[1;0.843;0]);
                set(handles.pb_mor_krylov_output,'BackgroundColor',[0.94;0.94;0.94]);
                
            else    %Hermite Interpolation

                set(handles.pb_mor_krylov_input,'Visible','off');
                set(handles.pb_mor_krylov_output,'Visible','off');
                
            end

            set(handles.cb_mor_krylov_hermite,'Enable','on');
        
        end
        
    else                                                %IRKA
        
        set(handles.uitable_mor_krylov_output,'Visible','off');
        set(handles.uitable_mor_krylov,'Visible','on');
        
        set(handles.pb_mor_krylov_input,'Visible','off');
        set(handles.pb_mor_krylov_output,'Visible','off');
        
        set(handles.cb_mor_krylov_hermite,'Enable','off');
    end

    
%Auxiliary-functions for MOR

function [] = displaySystemInformation(object,sys)
%Displays informations about the system sys on the text-object object.
%Prevents the display of the full name of the system which could include
%full path names

    systemName = sys.Name;
    sys.Name = '';
    set(object, 'String', sys.disp);
    sys.Name = systemName;
      
function [] = suggestNamesMOR(sysName,handles)
%Sets the suggests names for the reduced system and the projection-matrices
%in the model-order-reduction menue 

    %Set the name of the system
    
    if ~isempty(sysName)

        splittedString = strsplit(sysName,'_');

        if length(splittedString) >= 2

           name = char(strcat('sysr_',splittedString(1,2)));

        else

            name = char(strcat('sysr_',sysName));

        end

        %Get the selected reduction method 

        if get(handles.pu_mor_method,'Value') == 1

            method = 'tbr';

        elseif get(handles.pu_mor_method,'Value') == 2

            method = 'modal';

        elseif get(handles.pu_mor_krylov_algorithm,'Value') == 1

            method = 'irka';

        else

            method = 'rk';

        end

        name = strcat(name,'_',method);

        %Set the names to the text-fields of the GUI-objects

        suggestVarname(name,handles.ed_mor_sysred);
        suggestVarname(sprintf('%s_w',name),handles.ed_mor_w);
        suggestVarname(sprintf('%s_v',name),handles.ed_mor_v);
        suggestVarname(sprintf('%s_shifts',name),handles.ed_mor_saveShifts);
        suggestVarname(sprintf('%s_hsv',sysName),handles.ed_mor_saveHsv);
    
    else
       
        set(handles.ed_mor_sysred,'String','');
        set(handles.ed_mor_w,'String','');
        set(handles.ed_mor_v,'String','');
        set(handles.ed_mor_saveShifts,'String','');
        
    end
    
function updatedHandles = saveHankelSingularValues(sys,sysName,hsv,R,L,handles)
%Stores the HankelSingularValues hsv to the table saved in handles

    % check if values for this system already exist and override them if
    % the case
    for i = 1:size(handles.storedHsv,1)
       if strcmp(handles.storedHsv{i,1},sysName) && handles.storedHsv{i,2} == size(sys.A,1)
           handles.storedHsv{i,3} = hsv;
           handles.storedHsv{i,4} = R;
           handles.storedHsv{i,5} = L;
           guidata(handles.figure1,handles);
           updatedHandles = handles;
           return;
       end
    end
    
    % add a new row containing the system-name, the system order and the
    % HankelSingularValues to the table
    index = size(handles.storedHsv,1)+1;
    handles.storedHsv{index,1} = sysName;
    handles.storedHsv{index,2} = size(sys.A,1);
    handles.storedHsv{index,3} = hsv;
    handles.storedHsv{index,4} = R;
    handles.storedHsv{index,5} = L;
    updatedHandles = handles;
    guidata(handles.figure1,handles);
    
function [success,hsv,R,L] = getHankelSingularValues(sys,sysName,handles)
%Reads out the HankelSingularValues for the given system from the table
%stored in the handles-structure

    % check if values for this system exist and read them out if the case
    for i = 1:size(handles.storedHsv,1)
       if strcmp(handles.storedHsv{i,1},sysName) && handles.storedHsv{i,2} == size(sys.A,1)
           hsv = handles.storedHsv{i,3};
           R = handles.storedHsv{i,4};
           L = handles.storedHsv{i,5};
           success = 1;
           return;
       end
    end   
    
    % the HankelSingularValues for this system are not stored
    success = 0;
    hsv = [];
    R = [];
    L = [];
    
        
    
%Functions which list variables from workspace    
    
function x = systemsInWorkspace()
% finds and lists all dynamical systems that are contained in workspace
    % read all variable names
    s=evalin('base', 'whos');
    % preallocate memory
    x=cell(length(s),1);
    for i=1:length(s)
        if strcmp(s(i).class,'sss') || strcmp(s(i).class,'ssRed') && ...
            ~strcmp(s(i).name,'load_var_dont_destroy') && ...
            ~strcmp(s(i).name,'GUI_dont_destroy')
            % save name
            x{i}=s(i).name;
        end
    end
    % remove empty (non-system) entries
    x(cellfun(@isempty,x)) = [];   
    
function x=matricesInWorkspace
% Returns all matrices from workspace

s=evalin('base', 'whos');
% preallocate memory
x=cell(length(s),1); %alle, auch quadtratische
for i=1:length(s)   
    if strcmp(s(i).class,'double') && length(s(i).size)==2 && any(s(i).size)
        % save name
        x{i}=s(i).name;
    end
end
% remove empty (non-system) entries
x(cellfun(@isempty,x)) = [];

function [sys, sysname] = getSysFromWs(namehandle)
% imports system from base workspace
% namehandle may be system name or handle to an edit/combo-control

if ishandle(namehandle)
    sysname = get(namehandle,'String');
    if iscell(sysname)
        sysname=sysname{get(namehandle,'Value')};
    end
else
    sysname=namehandle;
end
if isempty(sysname)
    return
end


sys = evalin('base', sysname); 

function x = listClassesInWorkspace(class)
%Finds and lists all objects of the given class from workspace

    % read all variable names
    s=evalin('base', 'whos');
    % preallocate memory
    x=cell(length(s),1);
    for i=1:length(s)
        if strcmp(s(i).class,class) && ...
            ~strcmp(s(i).name,'load_var_dont_destroy') && ...
            ~strcmp(s(i).name,'GUI_dont_destroy')
            % save name
            x{i}=s(i).name;
        end
    end
    
    % remove empty (non-system) entries
    x(cellfun(@isempty,x)) = []; 
    
function success = loadSystemWithLoadSss(filename,path)
% create a system from the matrices in a .mat file by using the function
% loadSss
   try
       lastwarn('');
       success = 1;
       
       % create a name for the system based on the filename
       splittedString = strsplit(filename,'.');
       name = char(strcat('sys_',splittedString(1,1)));

       count = 1;
       sTemp = name;

       % check wheater the name already exists in workspace
       while existInBaseWs(sTemp)~=0
            sTemp = strcat(name,num2str(count));
            count = count+1;
       end

       name = sTemp;
       
       % replace "-" with "_", because variable names are not allowed
       % contain "-"
       name = strrep(name,'-','_');

       % create system using loadSss
       sys = loadSss(strcat(path,filename));         
       assignin('base',name,sys);

       % check whether the system is DAE and warn the user if the case
       if sys.isDae
           msgbox(strcat(filename,': System is DAE. This User-Interface does not fully support systems in DAE-format'),'Warning','Warn');
           uiwait
       end    

       error('loadSss:WarningOccured',lastwarn);

    catch ex
        if strcmp(ex.identifier,'loadSss:WarningOccured')
            if ~isempty(ex.message)
                msgbox(strcat(filename,': ',ex.message),'Warning','warn');
                uiwait
            end
        else
            msgbox({strcat(filename,': '),'Error while evaluating function loadSss.', ...
                    'Try to load the matrices and then compose the model.'},'Error','error');
            success = 0;
        end            
   end
   
function [found,str] = extractEntryFromErrorStack(stack,funcName)
% This function extracts one paricular entry from the list of functions
% execption.stack which represents the chain of function calls that led to
% an error. "funcName" is the name of the function for which more detailed
% information about the error are desired (i.e. 'tbr')

    str = '';
    found = 0;

    for i = 1:length(stack)
        if strcmp(funcName,stack(i).name)
            str = [funcName,' (line ',num2str(stack(i).line),'):'];
            found = 1;
            break;
        end
    end   

    
%Auxiliary-functions for plotting 

function c = plotOutput2cell(data)
%Save the plots from different in- and output-channels in an cell
%array

    if length(size(data)) == 2

        c = cell(1,size(data,2));

        for i = 1:size(data,2)
           c{1,i} = data(:,i); 
        end

    else %Length == 3

        c = cell(size(data,2),size(data,3));

        for i = 1:size(data,2)
           for j = 1:size(data,3)
              c{i,j} = data(:,i,j); 
           end
        end

    end
    
function [legendText] = suggestDefaultLegendText(handles,sysName,isMimo)
%Sets a default legend text if the corresponding text field is empty

    %Read out the system name if it is not given

    if isempty(sysName)
        list = get(handles.lb_PaV_selectedSystems,'String');
        sysName = list{get(handles.lb_PaV_selectedSystems,'Value')};
    end

    if ~isempty(sysName)

        %Find the stored data for this system
        
        data = [];
        
        for i = 1:length(handles.plotData)
           if strcmp(sysName,handles.plotData{i,1}.name)
              data = handles.plotData{i,1};
              break;
           end
        end
        
        %Replace the _'s in the string (because of latex interpreter)

        legendText = strrep(sysName,'_',' ');
        
        %Add the selected input and output channels to the legendText for
        %Mimo systems
        
        if ~isempty(data)
            if ~isempty(data.sizeInputs) || ~isempty(data.sizeOutputs)    %MIMO

                in = data.in;
                out = data.out;

                if ~strcmp(in,'all') && ~strcmp(out,'all')
                   legendText = strcat(legendText,' (In: ',in,' Out: ',out,')');
                end
            end     
        end
    end
    
function [variableName] = suggestPlotDataName(handles,sysName)
%Sets a default variable name for the plot data that should be saved to
%workspace

    %Read out the system name if it is not given

    if isempty(sysName)
        list = get(handles.lb_PaV_selectedSystems,'String');
        sysName = list{get(handles.lb_PaV_selectedSystems,'Value')};
    end

    %Get the selected plot-type

    type = '';
    
    selection = get(handles.plot_type,'Value');
    
    switch selection
        
        case 1      %Bode
            type = 'frd';
        case 2      %Magnitude
            type = 'frd';
        case 3      %Singular
            type = 'frd';
        case 4      %Impulse
            type = 'tf';
        case 5      %Step
            type = 'tf';
        case 6      %Pole-Zero-Map
            type = 'zpk';
    end
    
    variableName = strcat(sysName,'_',type);
    
    %Check if the constructed name is a valid variable name
    
    newName = variableName;
    counter = 1;
    
    while existInBaseWs(newName)
        
        newName = strcat(variableName,num2str(counter));
        counter = counter + 1;      
    end
        
    variableName = newName;
   
function [variableName] = suggestPlotDataNameTimeseries(handles,sysName)
%Sets a default variable name for the plot data that should be saved to
%workspace (for the timeseries data)

    %Read out the system name if it is not given

    if isempty(sysName)
        list = get(handles.lb_PaV_selectedSystems,'String');
        sysName = list{get(handles.lb_PaV_selectedSystems,'Value')};
    end

    %Add plot-type to the name for easier identification
    
    if get(handles.plot_type,'Value') == 4      %Impulse
        variableName = strcat(sysName,'_impulseData');
    else                                        %Step
        variableName = strcat(sysName,'_stepData');
    end
    
    %Check if the constructed name is a valid variable name
    
    newName = variableName;
    counter = 1;
    
    while existInBaseWs(newName)
        
        newName = strcat(variableName,num2str(counter));
        counter = counter + 1;      
    end
        
    variableName = newName;
    
function handles = savePlotData(handles)
%This function saves the values for the plot options of the currently
%selected system in the list (handles.lb_PaV_selectedSystems)

    %Get the name of the currently selected system
    
    list = get(handles.lb_PaV_selectedSystems,'String');
    systemName = list{get(handles.lb_PaV_selectedSystems,'Value'),1};
    
    if ~isempty(systemName)
        
       %Read out the values from the controls
       
       data.name = systemName;
       
       if strcmp(get(handles.panel_intoout,'Visible'),'on')
          listTemp = get(handles.pu_in,'String');
          data.in = listTemp{get(handles.pu_in,'Value'),1};
          listTemp = get(handles.pu_out,'String');
          data.out = listTemp{get(handles.pu_out,'Value'),1};
       else
          data.in = [];
          data.out = [];
       end
       
       data.save = get(handles.cb_PaV_SaveData,'Value');
       data.variableName = get(handles.et_PaV_saveData,'String');
       
       data.saveTimeseries =  get(handles.cb_PaV_SaveTimeseries,'Value');
       data.variableNameTimeseries = get(handles.et_PaV_saveTimeseries,'String');
        
       data.color = get(handles.pu_PaV_color,'Value');
       data.lineStyle = get(handles.pu_PaV_lineStyle,'Value');
       data.legendText = get(handles.ed_legend,'String');
       data.markerType = get(handles.pu_PaV_plotStyle_marker,'Value');
       
       if get(handles.rb_PaV_plotStyle_manual,'Value')
          data.plotStyle = 'manual'; 
       else
          data.plotStyle = 'auto';
       end
       
       if strcmp(get(handles.cb_PaV_SaveData,'Visible'),'on')   %System
          if get(handles.rb_manual,'Value')
               data.resolution = 'manual'; 

               if get(handles.rb_distlin,'Value')
                  data.distribution = 'linear';
               else
                  data.distribution = 'logarithmic'; 
               end

               data.min = get(handles.ed_min,'String');
               data.max = get(handles.ed_max,'String');
               data.steps = get(handles.sl_steps,'Value');
          else
             data.resolution = 'auto';
          end
       else                                                     %Object
          data.resolution = []; 
          
          data.save = [];
          data.variableName = [];
          
          data.saveTimeseries = [];
          data.variableNameTimeseries = [];
          
          data.distribution = []; 
          data.min = [];
          data.max = [];
          data.steps = [];   
       end
       
       %Check if there exists stored data for this system
       
       if isempty(handles.plotData)
          handles.plotData = {data}; 
       else
          for i = 1:length(handles.plotData)
              if strcmp(systemName,handles.plotData{i,1}.name)
                  
                 data.sizeInputs = handles.plotData{i,1}.sizeInputs;
                 data.sizeOutputs = handles.plotData{i,1}.sizeOutputs;
                 data.isSystem = handles.plotData{i,1}.isSystem;
                 
                 if strcmp(data.plotStyle,'auto')
                    data.color =  handles.plotData{i,1}.color;
                    data.markerType = handles.plotData{i,1}.markerType;
                    data.lineStyle = handles.plotData{i,1}.lineStyle;
                    data.legendText = handles.plotData{i,1}.legendText;
                 end
                 
                 if strcmp(data.resolution,'auto')
                    data.distribution = handles.plotData{i,1}.distribution; 
                    data.min = handles.plotData{i,1}.min;
                    data.max = handles.plotData{i,1}.max;
                    data.steps = handles.plotData{i,1}.steps;
                 end
                 
                 handles.plotData{i,1} = data;
                 guidata(handles.lb_PaV_selectedSystems,handles);
                 return;
              end
          end
          
          handles.plotData{end+1,1} = data;
          
       end
       
       guidata(handles.lb_PaV_selectedSystems,handles)
        
    end
    
function x = removeObjectsFromList(list,class)
%Removes all objects that are from the given class from the given list

    if ~isempty(list)

        if size(list,1) < size(list,2)
           list = list'; 
        end

        for i = 1:length(list) 
            
            try
               object = evalin('base',list{i,1});
               
               if isa(object,class)
                    list{i,1} = [];
               end
            end     
        end

        x = list;

        x(cellfun(@isempty,x)) = [];
    
    else
       x = list; 
    end

function [] = addRelativePaths()
%Add the subfolder from the GUI-folder-structure to the search-paths
    
    guiPath = which('sssMOR_App.m');
    pathArray = strsplit(guiPath,filesep);
    path = '';
    
    for i = 1:(size(pathArray,2)-1)
       
        if i > 1
           path = strcat(path,filesep,pathArray{1,i});
        else
           path = strcat(path,pathArray{1,i}); 
        end       
    end
    
    addpath(genpath(path));
    
    
%Auxiliary functions for operations on imported pictures
    
function [out] = imresize(im, out_dims)
% Implements the bilinear interpolation algorithm to resize an image. This
% function is used instead of the "imresize" function from Matlab, because
% otherwise it would be necessary to have the "Image acquisition" toolbox
% installed to start the App

    % Get some necessary variables first
    in_rows = size(im,1);
    in_cols = size(im,2);
    out_rows = out_dims(1);
    out_cols = out_dims(2);

    % Let S_R = R / R'        
    S_R = in_rows / out_rows;
    % Let S_C = C / C'
    S_C = in_cols / out_cols;

    % Define grid of co-ordinates in our image
    % Generate (x,y) pairs for each point in our image
    [cf, rf] = meshgrid(1 : out_cols, 1 : out_rows);

    % Let r_f = r'*S_R for r = 1,...,R'
    % Let c_f = c'*S_C for c = 1,...,C'
    rf = rf * S_R;
    cf = cf * S_C;

    % Let r = floor(rf) and c = floor(cf)
    r = floor(rf);
    c = floor(cf);

    % Any values out of range, cap
    r(r < 1) = 1;
    c(c < 1) = 1;
    r(r > in_rows - 1) = in_rows - 1;
    c(c > in_cols - 1) = in_cols - 1;

    % Let delta_R = rf - r and delta_C = cf - c
    delta_R = rf - r;
    delta_C = cf - c;

    % Final line of algorithm
    % Get column major indices for each point we wish
    % to access
    in1_ind = sub2ind([in_rows, in_cols], r, c);
    in2_ind = sub2ind([in_rows, in_cols], r+1,c);
    in3_ind = sub2ind([in_rows, in_cols], r, c+1);
    in4_ind = sub2ind([in_rows, in_cols], r+1, c+1);       

    % Now interpolate
    % Go through each channel for the case of colour
    % Create output image that is the same class as input
    out = zeros(out_rows, out_cols, size(im, 3));
    out = cast(out, class(im));

    for idx = 1 : size(im, 3)
        chan = double(im(:,:,idx)); %// Get i'th channel
        % Interpolate the channel
        tmp = chan(in1_ind).*(1 - delta_R).*(1 - delta_C) + ...
                       chan(in2_ind).*(delta_R).*(1 - delta_C) + ...
                       chan(in3_ind).*(1 - delta_R).*(delta_C) + ...
                       chan(in4_ind).*(delta_R).*(delta_C);
        out(:,:,idx) = cast(tmp, class(im));
    end


