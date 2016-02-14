function varargout = sssMOR_GUI(varargin)
% SSSMOR_GUI - nes at end of demos 
% ------------------------------------------------------------------
% USAGE:  TODO
%
% See also SSSMOR_GETTINGSTARTED, SSS.
%
% ------------------------------------------------------------------
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                    -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Heiko Panzer, Sylvia Cremer, Maria Cruz Varona, 
%               Alessandro Castagnotto
% Last Change:  11 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------


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
                   'gui_OpeningFcn', @sssMOR_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @sssMOR_GUI_OutputFcn, ...
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
        
        
function sssMOR_GUI_OpeningFcn(hObject, eventdata, handles, varargin)  %#ok<*INUSL>

    %Make latex code possible for all static text-fields whose tag names start
    %with "latex"  
    
    lbls = findall(hObject);
    
    for i=1:length(lbls)
        
          l = lbls(i);
          t = get(l,'tag');
          
          if strfind(t,'latex_') == 1
              
             % Get current text, position and tag
             t = 'nothing'; 
             set(l,'units','pixels');
             s = get(l,'string');
             p = get(l,'position');
             parent = get(l,'parent');
             
             % Remove the UICONTROL
             delete(l);
             
             % Replace it with an axis with an text object on it
             handles.(t) = axes('parent',parent,'units','pixels','position',[p(1) p(2) p(3) p(4)],'visible','off');
             handles.(t) = text(0,0.3,s,'interpreter','latex');            
          end
    end
    
    %Add the subfolder from the GUI-folder-structure to the search-paths
    
    guiPath = which('sssMOR_GUI.m');
    pathArray = strsplit(guiPath,'\');
    path = '';
    
    for i = 1:(size(pathArray,2)-1)
       
        if i > 1
           path = strcat(path,'\',pathArray{1,i});
        else
           path = strcat(path,pathArray{1,i}); 
        end       
    end
    
    %Add path with all subfolders
    
    addpath(genpath(path));
    
    %Set default-values for the variables in saved in handles
    
    handles.splash = splash('splash.jpg');
    handles.output = hObject;
    handles.virtgr_an_red_buttons=[handles.pb_an_stability,handles.pb_an_h2,handles.pb_an_hinf,handles.pb_an_simtime,handles.pb_an_h2error,handles.pb_an_hinferror,handles.pb_an_simtimeobt,handles.pb_an_simtimeorig];
    handles.allbuttons=[handles.pb_save,handles.pb_load,handles.pb_create,handles.pb_readascii];
    handles.zoom=[];
    handles.letzterpfad='*';

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
    evalin('base','stop(timerfind),delete(timerfind)')


function varargout = sssMOR_GUI_OutputFcn(hObject, eventdata, handles)
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

    p = get(hObject,'CurrentPoint');

    if p(1,1)>915 && p(1,1)<1044 && p(1,2)>730 && p(1,2)<1000
        %TUM-Logo (Header)
        set(gcf,'Pointer','hand'); 
    elseif p(1,1)> 5 && p(1,1)<75 && p(1,2)>730 && p(1,2)<1000 
        %Lehrstuhl-Logo (Header)
        set(gcf,'Pointer','hand');
    elseif p(1,1)>5 && p(1,1)<170 && p(1,2)>0 && p(1,2)<70
        %sssMOR-Logo (Footer)
        set(gcf,'Pointer','hand');
    elseif p(1,1)>915 && p(1,1)<1044 && p(1,2)>0 && p(1,2)<70
        %TUM-Logo (Footer)
        set(gcf,'Pointer','hand');
    elseif p(1,1)>250 && p(1,1)<480 && p(1,2)>0 && p(1,2)<70
        %MorLab-Logo (Footer)
        set(gcf,'Pointer','hand');
    elseif p(1,1)>580 && p(1,1)<850 && p(1,2)>0 && p(1,2)<70
        %Lehrstuhl-Logo (Footer)
        set(gcf,'Pointer','hand');
    else
        %Not over any logo
        set(gcf,'Pointer','arrow');
    end




%--------------------------------------------------------------------------
%                             MAIN-MENUE
%--------------------------------------------------------------------------

%Callbacks of the panel

function logo_tum_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
    % load TUM logo
    A=imread('pictures\Heading.png');
    h=image(A);
    set(h,'HitTest','off')
    set(hObject,'YDir','reverse');

function logo_tum_ButtonDownFcn(hObject, eventdata, handles)
    % link to web page
    p=get(hObject,'CurrentPoint');
    if p(1,1)>915 && p(1,1)<1044 && p(1,2)>5 && p(1,2)<80
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
    A=imread('pictures\MOR.jpg');
    image(A);
    set(hObject,'YDir','reverse');

function text_about_weblink_ButtonDownFcn(hObject, eventdata, handles)
    % Open a weblink if the mouse is klicked over the text element
    web https://www.rt.mw.tum.de/?sssMOR

    
    
%--------------------------------------------------------------------------
%                    POSTPROCESSING AND VISUALISATION
%--------------------------------------------------------------------------

function pu_in_Callback(hObject, eventdata, handles)
% MIMO input to be analyzed
% refresh list of suitable open figures
list_open_figures(handles)

function pu_out_Callback(hObject, eventdata, handles)
% MIMO output to be analyzed
% refresh list of suitable open figures
list_open_figures(handles)

function syschoice_Callback(hObject, eventdata, handles)
% occurs when system has been chosen

% get and convert system
try
    [sys, sysname] = get_sys_from_ws(hObject);
catch ex
    if strfind(ex.identifier, 'unassigned')
        set(handles.ed_legend,'String','');    
        return
    end
    errordlg(ex.message,'Error Dialog','modal')
    uiwait
    return
end

% x = get(hObject,'String');
% y = x{get(hObject,'Value')};
% if isempty(y)
%     set(handles.sysinfo, 'String', 'Please choose a system!')
%      return % no system chosen
% end
% % load system from workspace
% sys = evalin('base', y);
% if ~strcmp(class(sys), 'ss') && ~strcmp(class(sys), 'sss')
%     errordlg('Variable is not a valid state space model.','Error Dialog','modal')
%     set(handles.sysinfo,'String','Invalid model')
%     uiwait
%     return
% end

% show system information
set(handles.sysinfo, 'String', sys.disp);

% adapt pu_in and pu_out for MIMO systems
if sys.isMimo
    if get(handles.pu_in,'Value') > sys.m
        set(handles.pu_in,'Value', 1)
    end
    if get(handles.pu_out,'Value') > sys.p
        set(handles.pu_out,'Value', 1)
    end
    set(handles.panel_intoout,'Visible','on')
    in={num2str((1:sys.m)'),'all'};
    out={num2str((1:sys.p)'),'all'};
    set(handles.pu_in,'String',in)
    set(handles.pu_out,'String',out)
else
    % invisible for SISO
    set(handles.panel_intoout,'Visible','off')
end
% replay _ by space
% s= regexprep(sysname, '_', ' ');
% if ~isempty(sys.mor_info)
%     s=sprintf('%s; %s', s ,sys.mor_info.method);
% end
% suggest legend
%set(handles.ed_legend,'String',s);
% refresh list of open figures
list_open_figures(handles)


function plot_type_Callback(hObject, eventdata, handles)
if (get(hObject,'Value')==3 || get(hObject,'Value')==5)&&get(handles.figure,'Value')==1
    % make logarithmic checkbox visible for bode and frequency
    set(handles.st_scale,'Visible','on')
    set(handles.bg_x,'Visible','on')
    set(handles.bg_y,'Visible','on')
else
    % othervise invisible
    set(handles.st_scale,'Visible','off')
    set(handles.bg_x,'Visible','off')
    set(handles.bg_y,'Visible','off')
end
if get(hObject,'Value')==4
    % pzmap
    set(handles.panel_marker,'Visible','off')
    set(handles.rb_auto,'Value',1)
    set(handles.rb_auto,'Enable','inactive')
    set(handles.panel_manual,'Visible','off')
    set(handles.bg_distribution,'Visible','off')
    set(handles.rb_manual,'Enable','off')
else
    set(handles.panel_marker,'Visible','on')
    set(handles.rb_auto,'Enable','on')
    set(handles.rb_manual,'Enable','on')
    try
        x = get(handles.syschoice,'String');
        % system chosen for postprocessing
        y=x{get(handles.syschoice,'Value')};
        if ~isempty(y)
            % for MIMO, in and out must be able to be chosen
            sys = evalin('base', y);
            if ~strcmp(class(sys), 'sss')
                set(handles.sysinfo, 'String', 'Not an sss model')
                errordlg('Variable is not a valid state space model.')
                uiwait
                return
            end
            set(handles.sysinfo, 'String', sys.disp);
            if sys.isMimo
                if get(handles.pu_in,'Value')>size(sys.B,2)+1
                    set(handles.pu_in,'Value',1)
                end
                if get(handles.pu_out,'Value')>size(sys.C,1)+1
                    set(handles.pu_out,'Value',1)
                end
                set(handles.panel_intoout,'Visible','on')
                in={num2str((1:size(sys.B,2))'),'all'};
                out={num2str((1:size(sys.C,1))'),'all'};
                set(handles.pu_in,'String',in)
                set(handles.pu_out,'String',out)
            else
                set(handles.panel_intoout,'Visible','off')
            end
        else
            set(handles.sysinfo,'String','Please choose system')
        end
    catch ex %#ok<NASGU>
        set(handles.panel_intoout,'Visible','off')
    end
end
% refresh list of open figures
list_open_figures(handles)


function figure_Callback(hObject, eventdata, handles)
% handle of chosen figure
global fighand
userdata=get(hObject,'UserData');
fighand=userdata(get(hObject,'Value'));
if get(handles.plot_type,'Value')==3 ||get(handles.plot_type,'Value')==5
    if fighand~=0
        % figure already exists -> no choice of logarithmic
        set(handles.st_scale,'Visible','off')
        set(handles.bg_y,'Visible','off')
        set(handles.bg_x,'Visible','off')
    else
        % new figure -> log can be chosen
        set(handles.st_scale,'Visible','on')
        set(handles.bg_y,'Visible','on')
        set(handles.bg_x,'Visible','on')
    end
end
list_open_figures(handles); % aktualisieren, falls figures geschlossen wurden
if fighand ~=0
    figure(fighand)
end

function rb_auto_Callback(hObject, eventdata, handles)
% automatic choice of time/frequency
set(handles.panel_manual,'Visible','off')
set(handles.bg_distribution,'Visible','off')

function rb_manual_Callback(hObject, eventdata, handles)
% manual choice of time/frequency
set(handles.panel_manual,'Visible','on')
set(handles.bg_distribution,'Visible','on')

function ed_min_Callback(hObject, eventdata, handles)
% check user input
test_min(hObject,handles) % if ok, save 0 to userdata, otherwise 1
if get(handles.rb_manual,'Value')==0 % manual selected but inactive
    set(handles.rb_manual,'Value',1)
    test_max(handles.ed_max,handles)
end


function ed_max_Callback(hObject, eventdata, handles)
% see ed_min_Callback
test_max(hObject,handles)
if get(handles.rb_manual,'Value')==0
    set(handles.rb_manual,'Value',1)
    test_min(handles.ed_min,handles)
end


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
set(handles.st_curstep,'String',h)
% check min and max
if get(handles.rb_manual,'Value')==0
    set(handles.rb_manual,'Value',1)
    test_max(handles.ed_max,handles)
    test_min(handles.ed_min,handles)
end

function rb_colorst_Callback(hObject, eventdata, handles)
% reset inserted colors if standard color is chosen
set(handles.ed_r,'String','')
set(handles.ed_g,'String','')
set(handles.ed_b,'String','')
set(handles.ed_r,'UserData','0')
set(handles.ed_g,'UserData','0')
set(handles.ed_b,'UserData','0')


function rb_colorvek_Callback(hObject, eventdata, handles)
% standard: black
set(handles.ed_r,'String','0.0')
set(handles.ed_g,'String','0.0')
set(handles.ed_b,'String','0.0')

function ed_r_Callback(hObject, eventdata, handles)
% check input
test_colors(handles,hObject)
% if radio button for RGB is not selected, set it and unset other colors
if get(handles.rb_colorvek,'Value')==0
    set(handles.rb_colorvek,'Value',1)
    set(handles.ed_g,'UserData',1)
    set(handles.ed_b,'UserData',1)
end

function ed_g_Callback(hObject, eventdata, handles)
%see ed_r_Callback
test_colors(handles,hObject)
if get(handles.rb_colorvek,'Value')==0
    set(handles.rb_colorvek,'Value',1)
    set(handles.ed_r,'UserData',1)
    set(handles.ed_b,'UserData',1)
end


function ed_b_Callback(hObject, eventdata, handles)
% see ed_r_Callback
test_colors(handles,hObject)
if get(handles.rb_colorvek,'Value')==0
    set(handles.rb_colorvek,'Value',1)
    set(handles.ed_g,'UserData',1)
    set(handles.ed_r,'UserData',1)
end

function ed_post_markersize_Callback(hObject, eventdata, handles)
% check input data of linewidth
test_width(hObject)

function pb_plot_Callback(hObject, eventdata, handles)
% plot
set(hObject,'String','Busy')
set(hObject,'Enable','off') 
global fighand
% check whether a system was selected
if get(handles.syschoice,'Value')==1
    errordlg('Please choose a system.','Error Dialog','modal')
    set(hObject,'String','Plot')
    set(hObject,'Enable','on')
    uiwait
    return
end

% get and convert system
try
    [sys, sysname] = get_sys_from_ws(handles.syschoice);
catch ex
    set(hObject,'String','Plot')
    set(hObject,'Enable','on') 
    errordlg(ex.message,'Error Dialog','modal')
    uiwait
    return
end

% if time vector is given, check min, max
if get(handles.rb_auto,'Value')==1
    W=[];
elseif get(handles.rb_manual,'Value')==1
    if get(handles.ed_min,'UserData')==1
        errordlg('Please correct minimal frequency/time first','Error Dialog','modal')
        uiwait
        set(hObject,'String','Plot')
        set(hObject,'Enable','on')
        return
    elseif get(handles.ed_max,'UserData')==1
        errordlg('Please correct maximal frequency/time first','Error Dialog','modal')
        uiwait
        set(hObject,'String','Plot')
        set(hObject,'Enable','on')
        return
    end
    minimum=str2double(get(handles.ed_min,'String'));
    maximum=str2double(get(handles.ed_max,'String'));
    steps=get(handles.sl_steps,'Value');
    if get(handles.rb_distlin,'Value')==1
        W=linspace(minimum,maximum,steps);
    elseif get(handles.rb_distlog,'Value')==1
        %bei logspace muss der kleinste Wert>0 sein
        if str2double(get(handles.ed_min,'String'))==0
            errordlg('Minimum frequency/time must be greater than zero if you choose logarithmical distribution!','Error Dialog','modal')
            uiwait
            set(hObject,'String','Plot')
            set(hObject,'Enable','on')
            return
        end
        W=logspace(log10(minimum),log10(maximum),steps);
    else
        errordlg('Please select distribution of steps first','Error Dialog','modal')
        uiwait
        set(hObject,'String','Plot')
        set(hObject,'Enable','on')
        return
    end
else
    errordlg('Please choose resolution method first','Error Dialog','modal')
    uiwait
    set(hObject,'String','Plot')
    set(hObject,'Enable','on')
    return
end
% input color: standard colors are part of colorlist. selected color: col
if get(handles.rb_colorst,'Value')==1
    % standard color
    temp=get(handles.colorlist,'Value');
    list=get(handles.colorlist,'UserData');
    col=list(temp,:);
elseif get(handles.rb_colorvek,'Value')==1
    % rgb color
    if get(handles.ed_r,'UserData')==1 || ...
       get(handles.ed_g,'UserData')==1 || ...
       get(handles.ed_b,'UserData')==1
        errordlg('Please correct colorvector first','Error Dialog','modal')
        uiwait
        set(hObject,'Enable','on')
        set(hObject,'String','Plot')
        return
    end
    r=str2double(get(handles.ed_r,'String'));
    g=str2double(get(handles.ed_g,'String'));
    b=str2double(get(handles.ed_b,'String'));
    col=[r g b];
else
    errordlg('Please select color first','Error Dialog','modal')
    uiwait
    set(hObject,'String','Plot')
    set(hObject,'Enable','on')
    return
end
% check whether line width is correct: width
if get(handles.ed_width,'UserData')==1
    errordlg('Please correct Line Width','Error Dialog','modal')
    uiwait
    set(hObject,'String','Plot')
    set(hObject,'Enable','on')
    return
else
    width=str2double(get(handles.ed_width,'String'));
end
% check marker
if get(handles.pu_post_markerstyle,'Value')>1
    if get(handles.ed_post_markersize,'UserData')==1
        errordlg('Please correct Marker Size','Error Dialog','modal')
        uiwait
        set(hObject,'String','Plot')
        set(hObject,'Enable','on')
        return
    end
    val=get(handles.pu_post_markerstyle,'Value');
    k=get(handles.pu_post_markerstyle,'String');
    j=k(val);
    % first char of string is passed on to plot function
    l=regexp(j, '\ ', 'split');
    markerstyle=l{1}{1,1};
    val=get(handles.pu_post_markercolor,'Value');
    k=get(handles.pu_post_markercolor,'UserData');
    markercol=k(val,:);
    markersize=str2double(get(handles.ed_post_markersize,'String'));
else
    markerstyle='none';
    markercol=[0 0 0];
    markersize=5;
end
% save line style
switch get(handles.style,'Value')
    case 1
        line='-';
    case 2
        line='--';
    case 3
        line=':';
    case 4
        line='-.';
end
% select transfer channel
if strcmp(get(handles.panel_intoout,'Visible'),'on')
    %mimo
    input=get(handles.pu_in,'String');
    input=str2double(input{get(handles.pu_in,'Value')});
    output=get(handles.pu_out,'String');
    output=str2double(output{get(handles.pu_out,'Value')});
else
    input=1;
    output=1;
end
set(handles.figure1,'Pointer','watch')
drawnow
switch get(handles.plot_type,'Value')
case 1 %impulse response
    try
        [y, t]=impulse(sys(input, output),W);
        set(handles.ed_min,'String',t(1))
        set(handles.ed_min,'UserData',0)
        set(handles.ed_max,'String',t(end))
        set(handles.ed_max,'UserData',0)
    catch ex
        uiwait
        set(hObject,'Enable','on')
        set(hObject,'String','Plot')
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier,'nomem')
            errordlg('System is too big','Error Dialog','modal')
        elseif strfind(ex.message,'singular')
            errordlg('The system contains algebraic states.','Error Dialog','modal')
        else
            errordlg(ex.message, 'Error Dialog', 'modal')
        end
        uiwait
        return
    end
    if get(handles.figure,'Value')==1 || fighand==0
        fighand=figure;
    end
    figure(fighand)
    assignin('base',sysname,sys);
    UserData=get(fighand,'UserData');
    % use existing figure
    if get(fighand,'Tag')=='i'
        axes(UserData(1,1))
        hold on
        % merge new and old legend
        leg2=get(legend,'String');
        leg1={get(handles.ed_legend,'String')};
        leg=[leg2,leg1];
    else
        hold off
        leg=get(handles.ed_legend,'String');
    end
    % plot respective impulse responses in figure
    if ~isempty(UserData)
        for i=1:size(UserData,1)
            for j=1:size(UserData,2)
                subplot(UserData(i,j))
                hold on
                y_plot=y{i,j};
                plot(t,y_plot,'Color',col,'Linestyle',line,'LineWidth',width,'Marker',...
                    markerstyle,'MarkerFaceColor',markercol,'MarkerSize',...
                    markersize)
            end
        end
    else
        % plot into new figure
        UserData=zeros(size(y,1),size(y,2));
        k=1;
        for i=1:size(UserData,1)
            for j=1:size(UserData,2)
                UserData(i,j)=subplot(size(y,1),size(y,2),k);
                y_plot=y{i,j};
                plot(t,y_plot,'Color',col,'Linestyle',line,'LineWidth',width,'Marker',...
                    markerstyle,'MarkerFaceColor',markercol,'MarkerSize',...
                    markersize)
                k=k+1;
                if j==1 % label "To Out"
                    y_lab=sprintf('To Out(%i)',i);
                    ylabel(y_lab)
                end
                if i==1 % label "To In"
                    x_lab=sprintf('From In(%i)',j);
                    title(x_lab)
                end
            end
        end
        set(fighand,'UserData',UserData)
    end
    h=axes('position',[0,0,1,1],'Visible','off'); % axes unsichtbar in den hintergrund legen zur beschriftung
    text(0.5,0.98,'Impulse Response');
    text(0.5,0.02,'Time [s]')
    text(0.01,0.5,'Amplitude','Rotation',90)
    set(h,'HandleVisibility','off') % um zoomen usw zu vermeiden wird das handle unsichtbar gemacht
    hold on
    for i=1:size(UserData,1) % damit man auf die einzelnen subplots zugreifen kann werden alle nochmal angesprochen
        for j=1:size(UserData,2)
            axes(UserData(i,j))
        end
    end
    axes(UserData(1,1))
    legend(leg)
    set(fighand,'Tag','i')
case 2 %Step Response
    try % Sprungantwort berechnen, siehe auch impulsantwort
       [y, t]=step(sys(input, output),W);
       set(handles.ed_min,'String',t(1))
       set(handles.ed_min,'UserData',0)
       set(handles.ed_max,'String',t(end))
       set(handles.ed_max,'UserData',0)
    catch ex
        set(hObject,'Enable','on')
        set(hObject,'Strin','Plot')
        set(handles.figure1,'Pointer','arrow')
        if strfind(ex.identifier,'nomem')
            errordlg('System is too big','Error Dialog','modal')
        elseif strfind(ex.message,'singular')
            errordlg('The system contains algebraic states.','Error Dialog','modal')
        else
            errordlg(ex.message,'Error Dialog','modal')
        end
        uiwait
        return
    end
    if get(handles.figure,'Value')==1 || fighand==0
        figure
        fighand=gcf;
    end
    figure(fighand)
    assignin('base',sysname,sys); % neue systemdaten in den workspace schreiben
    UserData=get(fighand,'UserData');    
    if get(fighand,'Tag')=='s' % vorhandenes figure nutzen
        axes(UserData(1,1)) 
        hold on
        leg2=get(legend,'String');
        leg1={get(handles.ed_legend,'String')};
        leg=[leg2,leg1];
    else
    hold off
    leg=get(handles.ed_legend,'String');
    end
    if ~isempty(UserData)
         for i=1:size(UserData,1)
             for j=1:size(UserData,2)
                 subplot(UserData(i,j))
                 hold on
                 y_plot=y{i,j};
                 plot(t,y_plot,'Color',col,'Linestyle',line,'LineWidth',width,'Marker',...
                     markerstyle,'MarkerFaceColor',markercol,'MarkerSize',...
                     markersize)
             end
         end
    else
        UserData=zeros(size(y,1),size(y,2));
        k=1;
        for i=1:size(UserData,1)
            for j=1:size(UserData,2)
                UserData(i,j)=subplot(size(y,1),size(y,2),k);
                y_plot=y{i,j};
                plot(t,y_plot,'Color',col,'Linestyle',line,'LineWidth',width,'Marker',...
                    markerstyle,'MarkerFaceColor',markercol,'MarkerSize',...
                    markersize)
                k=k+1;
                if j==1
                    y_lab=sprintf('To Out(%i)',i);
                    ylabel(y_lab)
                end
                if i==1
                    x_lab=sprintf('From In(%i)',j);
                    title(x_lab)
                end
            end
        end
        set(fighand,'UserData',UserData)
    end
    h=axes('position',[0,0,1,1],'Visible','off');
    text(0.4,0.98,'Step Response');
    text(0.5,0.02,'Time [s]')
    text(0.01,0.5,'Amplitude','Rotation',90)
    set(h,'HandleVisibility','off')
    hold on
    for i=1:size(UserData,1)
        for j=1:size(UserData,2)
            axes(UserData(i,j))
        end
    end
    axes(UserData(1,1))
    legend(leg)
    set(fighand,'Tag','s')

case 3 %Bode
    try
        [mag, phase, omega]=bode(sys(input, output),W);
        set(handles.ed_min,'String',omega(1))
        set(handles.ed_min,'UserData',0)
        set(handles.ed_max,'String',omega(end))
        set(handles.ed_max,'UserData',0)
    catch ex
        errordlg(ex.message,'Error Dialog','modal')
        uiwait
        set(hObject,'String','Plot')
        set(handles.figure1,'Pointer','arrow')
        set(hObject,'Enable','on')
        return
    end
    if get(handles.figure,'Value')==1 || fighand==0
        figure
        fighand=gcf;
    end
    figure(fighand)
    assignin('base',sysname,sys);
    UserData=get(fighand,'UserData');
    if ~isempty(strfind(get(fighand,'Tag'),'b'))
        axes(UserData(1,1)) 
        hold on
        leg2=get(legend,'String');
        leg1={get(handles.ed_legend,'String')};
        leg=[leg2,leg1];
    else
        hold off
        leg=get(handles.ed_legend,'String');
    end
    if ~isempty(UserData)% mag und phase plotten in vorhandene fig plotten
        if ~isempty(strfind(get(fighand,'Tag'),'O'))
            mag=cellfun(@(x) 20*log10(x),mag,'UniformOutput',false);
        end
         for i=1:size(UserData,1)
             for j=1:size(UserData,2)
                 subplot(UserData(i,j))
                 hold on
                 if rem(i,2)==1 %amplitude plotten
                     y_plot=mag{ceil(i/2),j};
                     plot(omega,y_plot,'Color',col,'Linestyle',line,'LineWidth',width,'Marker',...
                         markerstyle,'MarkerFaceColor',markercol,'MarkerSize',...
                         markersize)
                 else %phase plotten
                     y_plot=phase{i/2,j};
                     plot(omega,y_plot,'Color',col,'Linestyle',line,'LineWidth',width,'Marker',...
                         markerstyle,'MarkerFaceColor',markercol,'MarkerSize',...
                         markersize)
                 end
             end
         end
    else
        if get(handles.rb_ylog,'Value')==1 %y logarithmisch in dB
            mag=cellfun(@(x) 20*log10(x),mag,'UniformOutput',false);
            set(fighand,'Tag','bO') % O für logarithmische skalierung
        else
            set(fighand,'Tag','bI') % I für lineare skalierung
        end
        UserData=zeros(2*size(mag,1),size(mag,2));
        k=1;
        for i=1:size(UserData,1)
            for j=1:size(UserData,2)
                UserData(i,j)=subplot(2*size(mag,1),size(mag,2),k);
                if rem(i,2)==1 %amplitude plotten
                    y_plot=mag{ceil(i/2),j};
                    plot(omega,y_plot,'Color',col,'Linestyle',line,'LineWidth',width,'Marker',...
                        markerstyle,'MarkerFaceColor',markercol,'MarkerSize',...
                        markersize)
                    if get(handles.rb_xlog,'Value')==1
                        set(gca, 'XScale', 'log');
                    end
                else %phase plotten
                    y_plot=phase{i/2,j};
                    plot(omega,y_plot,'Color',col,'Linestyle',line,'LineWidth',width,'Marker',...
                        markerstyle,'MarkerFaceColor',markercol,'MarkerSize',...
                        markersize)
                    if get(handles.rb_xlog,'Value')==1
                        set(gca, 'XScale', 'log');
                    end
                end
                if j==1
                    y_lab=sprintf('To Out(%i)',ceil(i/2));
                    ylabel(y_lab)
                end
                if i==1
                    x_lab=sprintf('From In(%i)',j);
                    title(x_lab)
                end
                k=k+1;
            end
        end
        set(fighand,'UserData',UserData)
    end
    h=axes('position',[0,0,1,1],'Visible','off');
    text(0.4,0.98,'Bode Diagram');
    text(0.5,0.02,'Frequency [rad/sec]')
    if get(handles.rb_ylog,'Value')==1
        text(0.01,0.5,'Magnitude [dB] ; Phase [deg]','Rotation',90)
    else
        text(0.01,0.5,'Magnitude (abs) ; Phase [deg]','Rotation',90)
    end
    set(h,'HandleVisibility','off')
    hold on
    for i=1:size(UserData,1)
        for j=1:size(UserData,2)
            axes(UserData(i,j))
        end
    end
    axes(UserData(1,1))
   legend(leg)

case 4 %Pole-Zero Map

   try
       [p,z]=pzmap(sys,input,output);
   catch ex
        if strfind(ex.identifier,'nomem')
            errordlg('Out of memory, system is too big','Error Dialog','modal')
            uiwait
            set(handles.figure1,'Pointer','arrow')
            set(hObject,'String','Plot')
            set(hObject,'Enable','on')
            return
        elseif strfind(ex.message,'singular')
            errordlg(ex.message,'Error Dialog','modal')
            uiwait
            set(handles.figure1,'Pointer','arrow')
            set(hObject,'String','Plot')
            set(hObject,'Enable','on')
            return
        else
            throw(ex)
        end
   end
    assignin('base',sysname,sys);
    if get(handles.figure,'Value')==1 || fighand==0
        figure
        fighand=gcf;
    end
    figure(fighand)
    UserData=get(fighand,'UserData');

    if ~isempty(strfind(get(fighand,'Tag'),'p'))
        axes(UserData(1,1)) 
        hold on
        leg2=get(legend,'String');
        leg1={get(handles.ed_legend,'String')};
        leg=[leg2,leg1];
    else
        hold off
        leg=get(handles.ed_legend,'String');
    end
    if ~isempty(UserData) % plot all pzmaps in open figure
        for i=1:size(UserData,1)
            for j=1:size(UserData,2)
                subplot(UserData(i,j))
                hold on
                z_plot=z{i,j};
                plot(real(z_plot),imag(z_plot),'--wo','MarkerEdgeColor',col)
                temp=plot(real(p),imag(p),'--wx','MarkerEdgeColor',col);
                hAnnotation = get(temp,'Annotation');
                hLegendEntry = get(hAnnotation','LegendInformation');
                set(hLegendEntry,'IconDisplayStyle','off')
            end
        end
    else
        UserData=zeros(size(z,1),size(z,2)); %plot in new figure
        k=1;
        for i=1:size(UserData,1)
            for j=1:size(UserData,2)
                UserData(i,j)=subplot(size(z,1),size(z,2),k);
                hold on
                z_plot=z{i,j};
                plot(real(z_plot),imag(z_plot),'--wo','MarkerEdgeColor',col)
                temp=plot(real(p),imag(p),'--wx','MarkerEdgeColor',col);
                hAnnotation = get(temp,'Annotation');
                hLegendEntry = get(hAnnotation','LegendInformation');
                set(hLegendEntry,'IconDisplayStyle','off')
                k=k+1;
                uy=min([min(imag(p)),min(imag(z_plot))]);
                oy=max([max(imag(p)),max(imag(z_plot))]);
                ux=min([min(real(p)),min(real(z_plot))]);
                ox=max([max(real(p)),max(real(z_plot))]);
                temp=plot([0 0],[1.1*uy 1.1*oy],':k');
                hAnnotation = get(temp,'Annotation');
                hLegendEntry = get(hAnnotation','LegendInformation');
                set(hLegendEntry,'IconDisplayStyle','off')
                temp=plot([1.1*ux 1.1*ox],[0 0],':k');
                hAnnotation = get(temp,'Annotation');
                hLegendEntry = get(hAnnotation','LegendInformation');
                set(hLegendEntry,'IconDisplayStyle','off')
                if j==1 % an die äußere Seite jeweils To Out schreiben
                    y_lab=sprintf('To Out(%i)',i);
                    ylabel(y_lab)
                end
                if i==1 % an die obere Seite jeweils To In schreiben
                    x_lab=sprintf('From In(%i)',j);
                    title(x_lab)
                end
            end
        end
        set(fighand,'UserData',UserData) % handles der axes in UserData der figure speichern
    end
    h=axes('position',[0,0,1,1],'Visible','off'); % axes unsichtbar in den hintergrund legen zur beschriftung
    text(0.4,0.98,'Pole-Zero Map');
    text(0.5,0.02,'Real Axis')
    text(0.01,0.5,'Imaginary Axis','Rotation',90)
    set(h,'HandleVisibility','off') % um zoomen usw zu vermeiden wird das handle unsichtbar gemacht
    hold on
    for i=1:size(UserData,1) % damit man auf die einzelnen subplots zugreifen kann werden alle nochmal angesprochen
        for j=1:size(UserData,2)
            axes(UserData(i,j))
        end
    end
    axes(UserData(1,1))
    legend(leg)
    set(fighand,'Tag','p')
case 5 %Frequency response
    try
        [mag, phase, omega]=bode(sys(input, output),W);
        set(handles.ed_min,'String',omega(1))
        set(handles.ed_min,'UserData',0)
        set(handles.ed_max,'String',omega(end))
        set(handles.ed_max,'UserData',0)
    catch ex
        errordlg(ex.message,'Error Dialog','modal')
        uiwait
        set(hObject,'String','Plot')
        set(handles.figure1,'Pointer','arrow')
        return
    end
    if get(handles.figure,'Value')==1 || fighand==0
        figure
        fighand=gcf;
    end
    figure(fighand)
    assignin('base',sysname,sys);
    UserData=get(fighand,'UserData');
    if ~isempty(strfind(get(fighand,'Tag'),'f'))
        axes(UserData(1,1)) 
        hold on
        leg2=get(legend,'String');
        leg1={get(handles.ed_legend,'String')};
        leg=[leg2,leg1];
    else
        hold off
        leg=get(handles.ed_legend,'String');
   end
   if ~isempty(UserData)% mag und phase plotten in vorhandene fig plotten
       if ~isempty(strfind(get(fighand,'Tag'),'O'))
           mag=cellfun(@(x) 20*log10(x),mag,'UniformOutput',false);
       end
        for i=1:size(UserData,1)
            for j=1:size(UserData,2)
                subplot(UserData(i,j))
                hold on
                y_plot=mag{i,j};
                plot(omega,y_plot,'Color',col,'Linestyle',line,'LineWidth',width,'Marker',...
                    markerstyle,'MarkerFaceColor',markercol,'MarkerSize',...
                    markersize)
            end
        end
    else
       if get(handles.rb_ylog,'Value')==1 %y logarithmisch in dB
           mag=cellfun(@(x) 20*log10(x),mag,'UniformOutput',false);
           set(fighand,'Tag','fO') % O für logarithm skalierung
       else
           set(fighand,'Tag','fI') % I für lineare skalierung
       end
        UserData=zeros(size(mag));
        k=1;
        for i=1:size(UserData,1)
            for j=1:size(UserData,2)
                UserData(i,j)=subplot(size(mag,1),size(mag,2),k);
                    y_plot=mag{i,j};
                    plot(omega,y_plot,'Color',col,'Linestyle',line,'LineWidth',width,'Marker',...
                        markerstyle,'MarkerFaceColor',markercol,'MarkerSize',...
                        markersize)
                    if get(handles.rb_xlog,'Value')==1
                        set(gca, 'XScale', 'log');
                    end
                if j==1
                    y_lab=sprintf('To Out(%i)',ceil(i/2));
                    ylabel(y_lab)
                end
                if i==1
                    x_lab=sprintf('From In(%i)',j);
                    title(x_lab)
                end
                k=k+1;
            end
        end
        set(fighand,'UserData',UserData)
    end
    h=axes('position',[0,0,1,1],'Visible','off');
    text(0.4,0.98,'Frequency Response');
    text(0.5,0.02,'Frequency [rad/sec]')
    if get(handles.rb_ylog,'Value')==1 %y logarithmisch in dB
        text(0.01,0.4,'Magnitude [dB]','Rotation',90)
    else
        text(0.01,0.4,'Magnitude (abs)','Rotation',90)
    end
    set(h,'HandleVisibility','off')
    hold on
    for i=1:size(UserData,1)
        for j=1:size(UserData,2)
            axes(UserData(i,j))
        end
    end
    axes(UserData(1,1))
    legend(leg)
end
set(hObject,'String','Plot')
set(hObject,'Enable','on')
set(handles.figure1,'Pointer','arrow')
list_open_figures(handles)
figure(fighand)


function colorlist_Callback(hObject, eventdata, handles)
% if wrong radio butten active, select right one
if get(handles.rb_colorst,'Value')==0
    rb_colorst_Callback(handles.rb_colorst,eventdata, handles)
    set(handles.rb_colorst,'Value',1)
end

function ed_width_Callback(hObject, eventdata, handles)
% check line width
test_width(hObject)





%--------------------------------------------------------------------------
%                     LOADING AND SETTING UP MODELS
%--------------------------------------------------------------------------

%Callbacks of the buttons

function pb_load_Callback(hObject, eventdata, handles)

    % disable to avoid double call of uigetfile
    set(handles.allbuttons,'Enable','off')
    filename=sprintf('%s.mat',handles.letzterpfad);
    [filename,path]=uigetfile(filename);
    set(handles.allbuttons,'Enable','on')
    if filename==0
         return
    end
    handles.letzterpfad=path;
    guidata(hObject, handles);
    % only mat-files can be loaded
    if isempty(strfind(filename,'.mat'))
        errordlg('Only .mat files allowed.','Error Dialog','modal')
        uiwait
        return
    end 
    
    loadingSuccessfull = 1;
    
    %Check which option is selected
    
    if get(handles.rb_loadOptions_matrices,'Value') == 1
       
        %Load matrices to workspace
        
        evalin('base',sprintf('load(''%s%s'');',path,filename));
        
    else
        
        %Create system with the function loadSss
        
        try
            
           %Create a name for the system based on the filename
            
           splittedString = strsplit(filename,'.');
           name = char(strcat('sys_',splittedString(1,1)));
           
           count = 1;
           sTemp = name;
                    
           %Check wheater the name already exists in workspace
           
           while existInBaseWs(sTemp)~=0
                sTemp = strcat(name,num2str(count));
                count = count+1;
           end
           
           name = sTemp;
           
           %Create system using loadSss
          
           evalin('base',strcat(name,sprintf(' = loadSss(''%s%s'');',path,filename)));
           
           %Check whether the system is DAE and warn the user if the case
           
           sys = evalin('base', name);
           
           if sys.isDae
               msgbox('System is DAE. This User-Interface does not fully support systems in DAE-format','Warning','Warn');
               uiwait
           end           
            
        catch
            
            msgbox('Error while evaluating function loadSss','Error','error');
            loadingSuccessfull = 0;
            
        end      
    end
    
    %Refresh the display of the variables in workspace
    
    set(handles.lb_systems,'Value',[])
    set(handles.lb_systems,'String',systemsInWorkspace)
    set(handles.lb_matrixes,'Value',[])
    set(handles.lb_matrixes,'String',matricesInWorkspace)
    
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
    matname=entermatrixname(suggestion);
    if isempty(matname)
        return
    end
    % import matrix
    try
        A=readMatrixMarket([path filename]);
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
    assignin('base',matname,A)
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

    infoBox({'pictures\InfoLoadOptions.png'});
    uiwait;

    
%Callbacks of the list-boxes
    
function lb_systems_Callback(hObject, eventdata, handles)
    set(hObject,'String',systemsInWorkspace)

function lb_systems_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String',systemsInWorkspace)

function lb_matrixes_Callback(hObject, eventdata, handles)
    set(hObject,'String',matricesInWorkspace)

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
    filename=sprintf('save(''%s%s'',temp{:});',path,filename);
    % save variable list in base workspace
    assignin('base','temp', variables)    %*** ??
    % store selected variables
    evalin('base',filename)
    % delete temp variables from base workspace
    evalin('base', 'clear temp')  %*** ??
    msgbox('Saving was successful.','Information','modal')
    uiwait
end

function pb_delete_Callback(hObject, eventdata, handles)
% delete variables

string=get(handles.lb_systems,'String'); % all variables in listbox
value=get(handles.lb_systems,'Value'); % all selected variables
variables=cell(length(value));
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

sys=evalin('base',y);
if ~strcmp(class(sys), 'sss')
    errordlg('Variable is not a valid state space model.','Error Dialog','modal')
    uiwait
    return
end
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
        return
    end

    sys = evalin('base', y);
    if ~strcmp(class(sys), 'sss')
        try
            sys = sss(sys);
        catch ex %#ok<NASGU>
            set(handles.pb_mor_reduce,'Enable','off')
            set(handles.panel_mor_hsv,'Visible','off')
            set(handles.st_mor_sysinfo,'String','Invalid model')
            errordlg('Variable is not a valid state space model.','Error Dialog','modal')
            uiwait
            return
        end
    end

    % set max of slider to system order
    if get(handles.sl_mor_q,'Value') > sys.n || ... 
       str2double(get(handles.ed_mor_q,'String')) > sys.n
        if get(handles.rb_mor_krylov_multi,'Value')==1
            warndlg('Order of reduced system must not be greater than order of original system.','Warning','modal')
            uiwait
        end
        set(handles.sl_mor_q,'Value',1)
        set(handles.ed_mor_q,'String','1')
    end
    
    displaySystemInformation(handles.st_mor_sysinfo,sys);
    
    set(handles.sl_mor_q,'Max',size(sys.A,1));
    set(handles.pb_mor_reduce,'Enable','on');
    
    % suggest names for reduced system and projection matrices
    
    suggestNamesMOR(y,handles);

    % hide plot of Hankel Singular Values
    set(handles.panel_mor_hsv,'Visible','off')
    
    %Set new default values for expension points for krylov
    
    value = get(handles.uitable_mor_krylov,'Data');
   
    if size(value,1) == 1 && (isempty(cell2mat(value(1,1))) && isempty(cell2mat(value(1,2))) || ...
            value{1,1} == 0)
        
        x = get(handles.pu_mor_systems,'String');
        y = x{get(handles.pu_mor_systems,'Value')};
        
        if ~isempty(y)
           sys = evalin('base', y);
           redOrder = round(size(sys.A,1)/3); 
           value{1,1} = 0;
           value{1,2} = redOrder;
           set(handles.uitable_mor_krylov,'Data',value); 
           set(handles.uitable_mor_krylov_output,'Data',value); 
           uitable_mor_krylov_CellEditCallback(handles.uitable_mor_krylov, eventdata, handles);
        end        
    end      

function pb_refreshsys_Callback(hObject, eventdata, handles)

    %Refresh list of open figures and systems
    
    %Get all systems from workspace
    
    l=systemsInWorkspace();
    
    %Check wether the previous selected systems exist in workspace
    
    
    %Syschoisce
    
    l_alt = get(handles.syschoice,'String');

    if ~isempty(l_alt)

        sOld=l_alt{get(handles.syschoice,'Value')};
        indexOld=find(strcmp(sOld,l));

        if ~isempty(indexOld)
            set(handles.syschoice,'Value',indexOld);
        else
            set(handles.syschoice,'Value',1);
        end

    else
      set(handles.syschoice,'Value',1);
    end
    
      
    %System Analysis

    l_alt = get(handles.pu_an_sys,'String');

    if ~isempty(l_alt)

        sOld=l_alt{get(handles.pu_an_sys,'Value')};
        indexOld=find(strcmp(sOld,l));

        if ~isempty(indexOld)
            set(handles.pu_an_sys,'Value',indexOld);
        else
            set(handles.pu_an_sys,'Value',1);
        end

    else
      set(handles.pu_an_sys,'Value',1);
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
    
    
    %If there are no systems in workspace, set the list of selectable
    %options to the empty string, else set it to the list of systems

    if ~isempty(l) && size(l,1) >= 1

        set(handles.syschoice, 'String', l);  
        set(handles.pu_mor_systems,'String', l);
        set(handles.pu_an_sys,'String',l);
        
        %Display informations about the selected system
        
        pu_mor_systems_Callback(handles.pu_mor_systems, eventdata, handles);
        
%         sys = evalin('base',l{get(handles.pu_mor_systems,'Value')});
%         
%         if ~strcmp(class(sys), 'sss')
%             
%             try
%                 sys = sss(sys);
%                 %set(handles.st_mor_sysinfo,'String',sys.disp)
%                 displaySystemInformation(handles.st_mor_sysinfo,sys);
%             catch
%                 set(handles.st_mor_sysinfo,'String','Invalid model')
%             end
%             
%         else
%             %set(handles.st_mor_sysinfo,'String',sys.disp)
%             displaySystemInformation(handles.st_mor_sysinfo,sys);
%         end
        
    else
        set(handles.pu_mor_systems,'Value',1)

        set(handles.syschoice, 'String', [{''}; l]);  
        set(handles.pu_mor_systems,'String', [{''}; l]);
        set(handles.pu_an_sys,'String',[{''}; l]);

        set(handles.pb_mor_reduce,'Enable','off');
        set(handles.panel_mor_hsv,'Visible','off');
    end

    %refresh list of open figures
    
    list_open_figures(handles);
    
    %list vectors in workspace that might be s0
    
    set(handles.pu_mor_krylov_s0,'String',listS0InWorkspace);




return


% systems from list in workspace
x = get(handles.syschoice,'String');
try
    % system chosen for postprocessing
    y=x{get(handles.syschoice,'Value')};
    if ~isempty(y)
        sys = evalin('base', y);
        if ~strcmp(class(sys), 'sss')
            errordlg('Variable is not a valid state space model.')            
        else 
            set(handles.sysinfo, 'String', sys.disp);
            if sys.isMimo
                if get(handles.pu_in,'Value') > sys.m + 1
                    set(handles.pu_in,'Value', 1)
                end
                if get(handles.pu_out,'Value') > sys.p + 1
                    set(handles.pu_out,'Value', 1)
                end
                set(handles.panel_intoout,'Visible','on')
                in={num2str((1:size(sys.B,2))'),'all'};
                out={num2str((1:size(sys.C,1))'),'all'};
                set(handles.pu_in,'String',in)
                set(handles.pu_out,'String',out)
            else
                set(handles.panel_intoout,'Visible','off')
            end
        end
    else
        set(handles.sysinfo,'String','Please choose system')
        set(handles.panel_intoout,'Visible','off')
    end
catch ex  %#ok<NASGU>
    set(handles.panel_intoout,'Visible','off')
    set(handles.sysinfo,'String','')
end

% system chosen for MOR
try
    sys = get_sys_from_ws(handles.pu_mor_systems);
    set(handles.st_mor_sysinfo, 'String', sys.disp);
catch ex
    if strfind(ex.identifier, 'unassigned')
        set(handles.st_mor_sysinfo,'String','Please choose a system.')
        return
    end
    set(handles.st_mor_sysinfo,'String','Variable is not a sparse state space model.')
end

% system chosen for Analysis
try
    sys = get_sys_from_ws(handles.pu_an_sys);
    set(handles.tx_an_sysinfo, 'String', sys.disp);
catch ex
    if strfind(ex.identifier, 'unassigned')
        set(handles.tx_an_sysinfo,'String','Please choose a system.')
        return
    end
    set(handles.tx_an_sysinfo,'String','Variable is not a sparse state space model.')
end
pu_an_sys_Callback(handles.pu_an_sys, eventdata, handles)

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

%If Krylov is selected, the reduced order can't be specified by the slider 
%anymore

if get(hObject,'Value')==3
    set(handles.ed_mor_q,'Enable','off')
    set(handles.sl_mor_q,'Enable','off')
else
    set(handles.ed_mor_q,'Enable','on')
    set(handles.sl_mor_q,'Enable','on')
end

%If Krylov is selected, set the default values for the expension points

if get(hObject,'Value')==3
   
    value = get(handles.uitable_mor_krylov,'Data');
   
    if isempty(cell2mat(value(1,1))) && isempty(cell2mat(value(1,2)))
        
        x = get(handles.pu_mor_systems,'String');
        y = x{get(handles.pu_mor_systems,'Value')};
        
        if ~isempty(y)
           sys = evalin('base', y);
           redOrder = round(size(sys.A,1)/3); 
           value{1,1} = 0;
           value{1,2} = redOrder;
           set(handles.uitable_mor_krylov,'Data',value);
           set(handles.uitable_mor_krylov_output,'Data',value);
           uitable_mor_krylov_CellEditCallback(handles.uitable_mor_krylov, eventdata, handles);
        end        
    end    
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

    if (get(handles.ed_mor_w,'UserData')==1 && get(handles.cb_mor_savew,'Value')==1)...
     ||(get(handles.ed_mor_v,'UserData')==1 && get(handles.cb_mor_savev,'Value')==1)
        errordlg('Please correct names for projection matrices first','Error Dialog','modal')
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
    set(handles.figure1,'Pointer','watch')
    set(hObject,'Enable','off')
    drawnow % this makes the watch appear
    sysname = get(handles.pu_mor_systems,'String');
    sysname = sysname{get(handles.pu_mor_systems,'Value')};
    sys=evalin('base',sysname);

    % convert to sss
    if ~strcmp(class(sys), 'sss')
        try
            sys=sss(sys);
        catch ex
            set(hObject,'String','Plot')
            set(hObject,'Enable','on') 
            errordlg(['Original system is not a valid state space model: ' ex.message],'Error Dialog','modal')
            uiwait
            return
        end
    end

    % Reduce
    switch get(handles.pu_mor_method,'Value')
        
    case 1 %TBR
        if isempty(sys.HankelSingularValues) || isempty(sys.TBal) || isempty(sys.TBalInv)
            set(handles.figure1,'Pointer','arrow')
            set(hObject,'Enable','on')
            errordlg('Please calculate Hankel Singular Values first.','Error Dialog','modal')
            uiwait
            return
        end
        
        try
        
            if get(handles.rb_mor_tbrtruncate,'Value')==1
                % truncation

                Opts.type = 'tbr';

                [sysr, V, W] = tbr(sys, q, Opts);

                %TODO check what happend to that code
                %sysr.morInfo = struct('time', clock, 'method', 'TBR', 'orgsys', sysname);
            else
                % match DC gain

                Opts.type  = 'matchDcGain';

                [sysr, V, W] = tbr(sys, q, Opts);

            end
        
        catch ex
            set(handles.figure1,'Pointer','arrow')
            set(hObject,'Enable','on')
            errordlg(ex.message,'Error Dialog','balancing & truncation')
            uiwait
            return
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

        try
            [sysr, V, W] = modalMor(sys, q, opts);
            %sysr.mor_info=struct('time', clock, 'method', 'Modal Truncation', 'orgsys', sysname);
        catch ex
            set(handles.figure1,'Pointer','arrow')
            set(hObject,'Enable','on')
            errordlg(ex.message,'Error Dialog','modal')
            uiwait
            return
        end
        
        
    case 3 % Krylov
        
        
        %Get expansion points. on error s0:=[]
        
        [s_inp,s_out] = getExpansionPoints(handles);
        if isempty(s_inp) && isempty(s_out)
            set(handles.figure1,'Pointer','arrow')
            set(hObject,'Enable','on')
            return
        end
        
        test = get(handles.pu_mor_krylov_algorithm,'Value');
        
        switch get(handles.pu_mor_krylov_algorithm,'Value')
            
        case 3 % explicit moment matching
            
            if get(handles.rb_mor_krylov_input,'Value')==1
                % input krylov
               if sys.m>1 % multiple input?
                   if sys.p==1
                        set(handles.figure1,'Pointer','arrow')
                        set(hObject,'Enable','on')
                        errordlg('Your system has only one output, use output-Krylov subspace!','Error Dialog','modal')
                        uiwait
                        return  
                   else
                        set(handles.figure1,'Pointer','arrow')
                        set(hObject,'Enable','on')
                        errordlg('MIMO is not supported','Error Dialog','modal')
                        uiwait
                        return
                   end
               end
                method='Rational Krylov, Input Krylov Subspace';
            elseif get(handles.rb_mor_krylov_twosided,'Value')==1
                % output krylov
                if sys.p>1 % multiple output?
                   if sys.m==1
                        set(handles.figure1,'Pointer','arrow')
                        set(hObject,'Enable','on')
                        errordlg('Your system has only one input, use input-Krylov subspace!','Error Dialog','modal')
                        uiwait
                        return  
                   else
                        set(handles.figure1,'Pointer','arrow')
                        set(hObject,'Enable','on')
                        errordlg('MIMO is not supported','Error Dialog','modal')
                        uiwait
                        return
                   end
                end
                method='Rational Krylov, Output Krylov Subspace';
            else
                % two-sided
                if sys.isMimo
                    if sys.p>1 && sys.m==1
                        set(handles.figure1,'Pointer','arrow')
                        set(hObject,'Enable','on')
                        errordlg('Your system has only one input, use input-Krylov subspace!','Error Dialog','modal')
                        uiwait
                        return
                    elseif sys.m>1 && sys.p==1
                        set(handles.figure1,'Pointer','arrow')
                        set(hObject,'Enable','on')
                        errordlg('Your system has only one output, use output-Krylov subspace!','Error Dialog','modal')
                        uiwait
                        return  
                    else
                        set(handles.figure1,'Pointer','arrow')
                        set(hObject,'Enable','on')
                        errordlg('MIMO is not supported yet','Error Dialog','modal')
                        uiwait
                        return
                    end
                end
                method='Rational Krylov, two-sided';
            end
            try
                if get(handles.rb_mor_krylov_input,'Value')==1
                    [sysr, V, W] = rk(sys, s_inp);
                else
                    [sysr, V, W] = rk(sys, s_inp, s_out);
                end
                %sysr.mor_info=struct('time', clock, 'method', method, 'orgsys', sysname);
            catch ex
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg(['Reduction Failed: ' ex.message],'Error Dialog','modal')
                uiwait
                assignin('base','LastError',ex)
                return
            end

        case 2 %ICOP
            
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
            try
                [sysr,V,W,alpha_opt] = RK_ICOP(sys, s0,'maxiter',maxiter,'epsilon',epsilon,'projection',projection );
                %sysr.mor_info=struct('time', clock, 'method', 'Krylov, ICOP', 'orgsys', sysname);
                %sysr.mor_info.alpha_opt=alpha_opt;
            catch ex
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg('Failed','Error Dialog','modal')
                uiwait
                disp(ex.message)
                return
            end
            if isempty(sysr)
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg('ICOP did not converge, reduction failed','Error dialog','modal')
                uiwait
                return
            end

        case 1 %IRKA
            
            if sys.isMimo % MIMO?
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg('MIMO is not supported yet.','Error Dialog','modal')
                uiwait
                return
            end
            Opts.maxiter=str2double(get(handles.ed_mor_krylov_max,'String'));
            Opts.epsilon=str2double(get(handles.ed_mor_krylov_epsilon,'String'));
            if isempty(Opts.maxiter) || Opts.maxiter==0 || isempty(Opts.epsilon) %|| epsilon==0 ??
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg('Please choose the maximum number of iterations and epsilon first','Error dialog','modal')
                uiwait
                return
            end
            try
                [sysr,V,W,s0] = irka(sys, s_inp, Opts);
                %sysr.mor_info=struct('time', clock, 'method', 'Krylov, IRKA', 'orgsys', sysname);
                %sysr.mor_info.s0=s0;
            catch ex
                set(handles.figure1,'Pointer','arrow')
                set(hObject,'Enable','on')
                errordlg(ex.message,'Error Dialog','modal')
                uiwait
                return
            end
        end
    end

    % write system to workspace
    if get(handles.cb_mor_krylov,'Value')==1
        % impose E_r=I?
        sysr.resolve_dae
    end

    assignin('base',get(handles.ed_mor_sysred,'String'),sysr)
    if get(handles.cb_mor_savew,'Value')==1
        assignin('base',get(handles.ed_mor_w,'String'),W)
    end
    if get(handles.cb_mor_savev,'Value')==1
        assignin('base',get(handles.ed_mor_v,'String'),V)
    end
    set(handles.figure1,'Pointer','arrow')
    msgbox('Reduction was successful!','Information','modal')
    uiwait
    pb_refreshsys_Callback(hObject, eventdata, handles)
    pb_refreshlb_Callback(hObject, eventdata, handles)
    s=get(handles.ed_mor_sysred,'String');
    % suggest names for next reduction
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

%Callbacks for Balacing & Truncation

function updateTBR(hObject, eventdata, handles)

    %Show wait-cursor
    
    set(handles.figure1,'Pointer','watch')
    drawnow

    %Get the selected system

    sys_x=get(handles.pu_mor_systems,'String');
    sysname=sys_x{get(handles.pu_mor_systems,'Value')};
    sys=evalin('base',sysname);
    
    %Get Axes and desired model order
    
    hr=get(handles.axes_mor_hsv,'UserData');
    q=get(handles.sl_mor_q,'Value');
    
    if get(handles.pu_mor_method,'Value')==1 && ~isempty(sys.HankelSingularValues)
        
        %Calculate the signal-norms H_1 and H_inf
        
        e=2*sum(sys.HankelSingularValues((q+1):end));
        erel=e/max(sys.HankelSingularValues);
        set(handles.st_mor_tbr_error,'String',num2str(e, '%.3e'))
        set(handles.st_mor_tbr_relerror,'String',num2str(erel, '%.3e'))
        
        %Plot which shows the Hankel-Singular-Values
            
        set(handles.panel_mor_hsv, 'Visible','on') 
        cla(handles.axes_mor_hsv)
        hold(handles.axes_mor_hsv, 'on')

        if get(handles.rb_mor_tbr_norm,'Value')==1

            maxValue = max(sys.HankelSingularValues);
            h = plot(handles.axes_mor_hsv,sys.HankelSingularValues./maxValue);

        else
            h = plot(handles.axes_mor_hsv, sys.HankelSingularValues);
        end

        % make callback react to click on red HSV line

        set(h,'HitTest','off');

        % legend

        legend(handles.axes_mor_hsv, regexprep(sysname, '_', ' '));

        % set scale

        if get(handles.rb_mor_tbr_log,'Value')==1
            set(handles.axes_mor_hsv,'YScale','log')
        else
            set(handles.axes_mor_hsv,'YScale','linear')                
        end
        
        %Plot the red line
            
        if ishandle(hr)
            set(hr,'XData',[q,q])
            set(hr,'YData',get(handles.axes_mor_hsv,'YLim'));
        else
            %hr=plot(handles.axes_mor_hsv, [q,q],sys.HankelSingularValues([end,1]),'r');
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
    [sys, sysname] = get_sys_from_ws(handles.pu_mor_systems);
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
    tbr(sys,1);
    assignin('base', sysname, sys)
catch ex %***
    if strcmp(ex.identifier,'MATLAB:nomem')
        errordlg('Out of memory. System is too large to calculate Hankel Singular Values.','Error Dialog','modal')
    elseif strcmp(ex.identifier,'Control:foundation:LyapChol4')
        errordlg('A or (A,E) must have all their eigenvalues in the left-half plane.','Error Dialog','modal')
    else
        errordlg(ex.message,'Error Dialog','modal')
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
if isa(handles.zoom,'graphics.zoom')
    set(handles.zoom,'Enable','off') % deactivate zoom
end
guidata(hObject,handles)


function pb_info_tbr_Callback(hObject, eventdata, handles)

    infoBox({'pictures\InfoErrorBound.png'});
    uiwait;

%Callbacks for Modal

function ed_mor_modal_neig_Callback(hObject, eventdata, handles)
% check of double. replace ',' with '.'
x=str2num(get(hObject,'String')); %#ok<ST2NM>
if isempty(x) || length(x)>1 || imag(x)~=0
    if length(x)==2
        x=sprintf('%i.%i',round(x(1)),round(x(2)));
        set(hObject,'String',x)
        set(hObject,'UserData',0)
        return
    end
    errordlg('Has to be a real number.','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
else
    set(hObject,'UserData',0)
end

%Callbacks for Krylov

function pu_mor_krylov_algorithm_Callback(hObject, eventdata, handles)

    if get(hObject,'Value')==1
        % IRKA, nur twosided möglich
        %andere radio buttons(rbs) deaktivieren und Texte anpassen
        set(handles.rb_mor_krylov_twosided,'Value',1)
        set(handles.rb_mor_krylov_twosided,'Enable','inactive')
        set(handles.rb_mor_krylov_input,'Enable','off')
        set(handles.rb_mor_krylov_twosided,'Enable','off')
        set(handles.st_mor_krylov_points,'String','2. Starting points')
    else % sonst rbs aktivieren
        set(handles.rb_mor_krylov_twosided,'Enable','on')
        set(handles.rb_mor_krylov_input,'Enable','on')
        set(handles.rb_mor_krylov_output,'Enable','on')
        set(handles.st_mor_krylov_points,'String','2. Expansion points')
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
    
    if get(hObject,'Value')==3
        % explicit moment matching
        set(handles.panel_mor_krylov_iter,'Visible','off')
    else
        set(handles.panel_mor_krylov_iter,'Visible','on')
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

function pb_mor_krylov_reset_Callback(hObject, eventdata, handles)

if strcmp(get(handles.uitable_mor_krylov,'Visible'),'on')
    set(handles.uitable_mor_krylov,'Data',cell(1,2))
else
    set(handles.uitable_mor_krylov_output,'Data',cell(1,2))
end

function pb_mor_krylov_AddCompConj_Callback(hObject, eventdata, handles)

if strcmp(get(handles.uitable_mor_krylov,'Visible'),'on')
    x=get(handles.uitable_mor_krylov,'Data');
    xConj = conj(cell2mat(x));
    x = [x;mat2cell(xConj,ones(size(xConj,1),1),[1,1])];
    set(handles.uitable_mor_krylov,'Data',x);
else
    x=get(handles.uitable_mor_krylov_output,'Data');
    xConj = conj(cell2mat(x));
    x = [x;mat2cell(xConj,ones(size(xConj,1),1),[1,1])];
    set(handles.uitable_mor_krylov_output,'Data',x);
end

uitable_mor_krylov_CellEditCallback(handles.uitable_mor_krylov, eventdata, handles)

function uitable_mor_krylov_CellEditCallback(hObject, eventdata, handles)
    
try
    x=get(hObject,'Data');
    %x=addExpensionPoints(x,eventdata);
    %set(hObject,'Data',x)
    q=cell2mat(x);
    if size(q,2)==2
        q=sum(q(:,2));
        set(handles.st_mor_krylov_redOrder,'String',q);

        if get(handles.rb_mor_krylov_twosided,'Value') == 1    %Two sided
            set(handles.st_mor_krylov_matchedMom,'String',2*q);
        else                                               %One sided
            set(handles.st_mor_krylov_matchedMom,'String',q);
        end
        
        %Check weather the reduced order is bigger than the original order
        
        x = get(handles.pu_mor_systems,'String');
        y = x{get(handles.pu_mor_systems,'Value')};
        
        if ~isempty(y)
           sys = evalin('base',y);
           if sys.n < q
              msgbox('Reduced order is bigger than the original order of the system. Please correct that before reducing the system.','Warning','Warn');
              uiwait;
           end
        end
        
    end
end

function uitable_mor_krylov_CreateFcn(hObject, eventdata, handles)
% adapt size of table
set(hObject,'Data',cell(1,2))

function uitable_mor_krylov_output_CellEditCallback(hObject, eventdata, handles)

function uitable_mor_krylov_output_CreateFcn(hObject, eventdata, handles)
% adapt size of table
set(hObject,'Data',cell(1,2))

function bg_mor_krylov_side_SelectionChangedFcn(hObject, eventdata, handles)

    if get(handles.rb_mor_krylov_twosided,'Value')==1
        
        set(handles.uitable_mor_krylov_output,'Visible','Off');
        set(handles.uitable_mor_krylov,'Visible','On');
        
        if get(handles.rb_mor_krylov_hermite,'Value') == 0
        
            set(handles.pb_mor_krylov_input,'Visible','On');
            set(handles.pb_mor_krylov_output,'Visible','On');

            set(handles.pb_mor_krylov_input,'BackgroundColor',[1;0.843;0]);
            set(handles.pb_mor_krylov_output,'BackgroundColor',[0.94;0.94;0.94]);
        
        end
        
        set(handles.rb_mor_krylov_hermite,'Enable','On');
    else
        
        set(handles.uitable_mor_krylov_output,'Visible','Off');
        set(handles.uitable_mor_krylov,'Visible','On');
        
        set(handles.pb_mor_krylov_input,'Visible','Off');
        set(handles.pb_mor_krylov_output,'Visible','Off');
        
        set(handles.rb_mor_krylov_hermite,'Enable','Off');
    end

    uitable_mor_krylov_CellEditCallback(handles.uitable_mor_krylov, eventdata, handles);

function pb_mor_krylov_importVector_Callback(hObject, eventdata, handles)
    % if vector from workspace was selected, add it to table
    contents = cellstr(get(handles.pu_mor_krylov_s0,'String'));
    try
        x=evalin('base',contents{get(handles.pu_mor_krylov_s0,'Value')});
    catch ex
       errordlg('No valid matrix selected');
       return;
    end
    
    if size(x,2)~=2
        x=x';
    end
    if size(x,2)==2 && size(x,1)==2
       if ~(mod(x(1,2),1)==0 && mod(x(2,2),1)==0)
          if mod(x(1,1),1)==0 && mod(x(2,1),1)==0
              xTemp = x;
              x(1,1) = xTemp(1,2);
              x(2,1) = xTemp(2,2);
              x(1,2) = xTemp(1,1);
              x(2,2) = xTemp(2,1);
          elseif mod(x(1,1),1)==0 && mod(x(1,2),1)==0
              xTemp = x;
              x(1,1) = xTemp(2,1);
              x(2,1) = xTemp(2,2);
              x(1,2) = xTemp(1,1);
              x(2,2) = xTemp(1,2);
          elseif mod(x(2,1),1)==0 && mod(x(2,2),1)==0
              xTemp = x;
              x(1,2) = xTemp(2,1);
              x(2,1) = xTemp(1,2);
          end
       end        
    end
    x=num2cell(x);
    %x=addExpensionPoints(x,[]);
 
    if strcmp(get(handles.uitable_mor_krylov,'Visible'),'on')
        set(handles.uitable_mor_krylov,'Data',x);
    else
        set(handles.uitable_mor_krylov_output,'Data',x);
    end
    uitable_mor_krylov_CellEditCallback(handles.uitable_mor_krylov,eventdata,handles);    

function pb_mor_krylov_input_Callback(hObject, eventdata, handles)

%Set the right table visible

set(handles.uitable_mor_krylov,'Visible','On');
set(handles.uitable_mor_krylov_output,'Visible','Off');

%Change the color of the buttons to signal witch table is selected

set(handles.pb_mor_krylov_input,'BackgroundColor',[1;0.843;0]);
set(handles.pb_mor_krylov_output,'BackgroundColor',[0.94;0.94;0.94]);

function pb_mor_krylov_infoHermite_Callback(hObject, eventdata, handles)

    infoBox({'pictures\InfoHermiteInterpolation.png'});
    uiwait;

function pb_mor_krylov_output_Callback(hObject, eventdata, handles)

%Set the right table visible

set(handles.uitable_mor_krylov_output,'Visible','On');
set(handles.uitable_mor_krylov,'Visible','Off');

%Change the color of the buttons to signal witch table is selected

set(handles.pb_mor_krylov_output,'BackgroundColor',[1;0.843;0]);
set(handles.pb_mor_krylov_input,'BackgroundColor',[0.94;0.94;0.94]);

function rb_mor_krylov_hermite_Callback(hObject, eventdata, handles)

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

%--------------------------------------------------------------------------
%                            SYSTEM ANALYSIS
%--------------------------------------------------------------------------

function pb_an_stability_Callback(hObject, eventdata, handles)
% determine stability of system
set(handles.figure1,'Pointer','watch')
set(handles.virtgr_an_red_buttons,'Enable','off')
drawnow
x = get(handles.pu_an_sys,'String');
sysname=x{get(handles.pu_an_sys,'Value')};
if isempty(sysname)
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
    errordlg('Please choose system first!','Error Dialog','modal')
    uiwait
    return
end
sys=evalin('base',sysname);
if ~isfield(sys,'Stability')
    try
        sys=stability(sys);
    catch ex
        assignin('base','LastError',ex)
        throw(ex)
    end
end
set(handles.tx_an_stability,'String',sys.Stability)
assignin('base',sysname,sys)
if isempty(eventdata) || eventdata~=1
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
end

function pb_an_h2_Callback(hObject, eventdata, handles)
set(handles.figure1,'Pointer','watch')
set(handles.virtgr_an_red_buttons,'Enable','off')
drawnow
x = get(handles.pu_an_sys,'String');
sysname=x{get(handles.pu_an_sys,'Value')};
if isempty(sysname)
    if isempty(eventdata) || eventdata~=1
        errordlg('Please choose system first!','Error Dialog','modal')
        uiwait
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
    end
    return
end

sys=evalin('base', sysname);
if ~strcmp(class(sys), 'sss')
    try
        sys = sss(sys);
        assignin('base', sysname, sys);
    catch ex %#ok<NASGU>
        set(handles.pb_mor_reduce,'Enable','off')
        set(handles.panel_mor_hsv,'Visible','off')
        set(handles.st_mor_sysinfo,'String','Invalid model')
        errordlg('Variable is not a valid state space model.','Error Dialog','modal')
        uiwait
        return
    end
end

try
    h2 = norm(sys, 2);
    assignin('base', sysname, sys);
catch ex
    if strcmp(ex.identifier,'MATLAB:nomem')
        errordlg('Out of memory, system is too large to solve lyapunov equotation','Error Dialog','modal')
    elseif strcmp(ex.identifier,'Control:foundation:LyapChol4')
        errordlg('A or (A,E) must have all their eigenvalues in the left-half plane','Error Dialog','modal')
    else
        errordlg(ex.message)
    end
    uiwait
    if isempty(eventdata) || eventdata~=1
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
    end
    return
end
set(handles.tx_an_h2,'String',num2str(h2))
if isempty(eventdata) || eventdata~=1
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
end

function pb_an_hinf_Callback(hObject, eventdata, handles)
set(handles.figure1,'Pointer','watch')
drawnow
x = get(handles.pu_an_sys,'String');
sysname=x{get(handles.pu_an_sys,'Value')};
set(handles.virtgr_an_red_buttons,'Enable','off')
if isempty(sysname)
    if isempty(eventdata) || eventdata~=1
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
    end
    return
end

sys=evalin('base', sysname);
if ~strcmp(class(sys), 'sss')
    try
        sys = sss(sys);
    catch ex %#ok<NASGU>
        set(handles.pb_mor_reduce,'Enable','off')
        set(handles.panel_mor_hsv,'Visible','off')
        set(handles.st_mor_sysinfo,'String','Invalid model')
        errordlg('Variable is not a valid state space model.','Error Dialog','modal')
        uiwait
        return
    end
end

try
    hinf=norm(sys, inf);
catch ex
    errordlg(ex.message,'Error Dialog','modal')
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
    uiwait
    return
end

set(handles.tx_an_hinf,'String',num2str(hinf))
if isempty(eventdata) || eventdata~=1
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
end

function pb_an_h2error_Callback(hObject, eventdata, handles)
set(handles.figure1,'Pointer','watch')
set(handles.virtgr_an_red_buttons,'Enable','off')
drawnow
x = get(handles.pu_an_sys,'String');
sysname=x{get(handles.pu_an_sys,'Value')};
if isempty(sysname)
    if isempty(eventdata) || eventdata~=1
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
    end
    return
end
sysr=evalin('base',sysname);
if ~isfield(sysr,{'H_2_error_abs','H_2_error_rel_'})
    try
        syso=evalin('base',sysr.originalsysname);
    catch ex
        if strcmp('MATLAB:UndefinedFunction',ex.identifier)
            %Originalsystem liegt nicht im workspace
            errordlg('Originalsystem is not in the workspace','Error Dialog','modal')
            uiwait
            return
        else
            assignin('base','lastError',ex)
            throw(ex)
        end
    end
    
    % Fehler berechnen
    try
        h2e=norm(syso-sysr);
    catch ex
        if strcmp(ex.identifier,'MATLAB:nomem')
            errordlg('Out of memory, errorsystem is too big to calculate the lyapunov equotation!','Error Dialog','modal')
            uiwait
            if isempty(eventdata) || eventdata~=1
                    set(handles.figure1,'Pointer','arrow')
                    set(handles.virtgr_an_red_buttons,'Enable','on')
            end
            return
         elseif strcmp(ex.identifier,'Control:foundation:LyapChol4')
            errordlg('A or (A,E) must have all their eigenvalues in the left-half plane!','Error Dialog','modal')
            uiwait
            if isempty(eventdata) || eventdata~=1
                set(handles.figure1,'Pointer','arrow')
                set(handles.virtgr_an_red_buttons,'Enable','on')
            end
            return
        else
            throw(ex)
        end
    end
    if isempty(syso.H_2_norm)
        try
            h2o=norm(syso);
            assignin('base',sysr.originalsysname,syso)
        catch ex
            if strcmp(ex.identifier,'MATLAB:nomem')
                errordlg('Out of memory, errorsystem is too big to calculate the lyapunov equotation!','Error Dialog','modal')
                uiwait
                if isempty(eventdata) || eventdata~=1
                    set(handles.figure1,'Pointer','arrow')
                    set(handles.virtgr_an_red_buttons,'Enable','on')
                end
                return
            elseif strcmp(ex.identifier,'Control:foundation:LyapChol4')
                errordlg('A or (A,E) must have all their eigenvalues in the left-half plane!','Error Dialog','modal')
                uiwait
                if isempty(eventdata) || eventdata~=1
                    set(handles.figure1,'Pointer','arrow')
                    set(handles.virtgr_an_red_buttons,'Enable','on')
                end
                return
            else
                throw(ex)
            end
        end
    else
        h2o=syso.H_2_norm;
    end
    erel=h2e/h2o*100;
    sysr.H_2_error_abs=h2e;
    sysr.H_2_error_rel=erel;
    assignin('base',sysname,sysr)
else
    xe=sysr.H_2_error_abs;
    xrel=sysr.H_2_error_rel;
end
s=sprintf('%s %%',num2str(xrel));
set(handles.tx_an_h2errorabs,'String',xe)
set(handles.tx_an_h2errorrel,'String',s)
if isempty(eventdata) || eventdata~=1
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
end



function pb_an_hinferror_Callback(hObject, eventdata, handles)
set(handles.figure1,'Pointer','watch')
set(handles.virtgr_an_red_buttons,'Enable','off')
drawnow
x = get(handles.pu_an_sys,'String');
sysname=x{get(handles.pu_an_sys,'Value')};
if isempty(sysname)
    if isempty(eventdata) || eventdata~=1
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
    end
    return
end
sysr=evalin('base',sysname);
if ~isfield(sysr,{'H_inf_error_abs','H_inf_error_rel'})
    try
        syso=evalin('base',sysr.originalsysname);
    catch ex
        if strcmp('MATLAB:UndefinedFunction',ex.identifier)
            % original system is not in workspace
            errordlg('Original system is not in the workspace.','Error Dialog','modal')
            uiwait
            return
        else
            assignin('base','LastError',ex)
            throw(ex)
        end
    end
    %Fehler berechnen

    % Fehlersystem syse
    Ee=sparse([syso.E,zeros(length(syso.E),length(sysr.E));zeros(length(sysr.E),length(syso.E)),sysr.E]);
    Ae=sparse([syso.A,zeros(length(syso.A),length(sysr.A));zeros(length(sysr.A),length(syso.A)),sysr.A]);
    Be=sparse([syso.B;sysr.B]);
    Ce=sparse([syso.C,-sysr.C]);
    De=0;
    syse=struct('A',Ae,'B',Be,'C',Ce,'D',De,'E',Ee);
    xe=norm(syse,'inf');
    if ~isfield(syso,'H_inf_norm')
        
        [xo,syso]=norm(syso,'inf');
        assignin('base',sysr.originalsysname,syso)
    else
        xo=syso.H_inf_norm;
    end
    xrel=xe/xo*100;
    sysr.H_inf_error_abs=xe;
    sysr.H_inf_error_rel=xrel;
    assignin('base',sysname,sysr)
else
    xe=sysr.H_inf_error_abs;
    xrel=sysr.H_inf_error_rel;
end
s=sprintf('%s %%',num2str(xrel));
set(handles.tx_an_hinferrorabs,'String',xe)
set(handles.tx_an_hinferrorrel,'String',s)
if isempty(eventdata) || eventdata~=1
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
end

function pb_an_simtimeobt_Callback(hObject, eventdata, handles)
set(handles.figure1,'Pointer','watch')
set(handles.virtgr_an_red_buttons,'Enable','off')
drawnow
x = get(handles.pu_an_sys,'String');
sysname=x{get(handles.pu_an_sys,'Value')};
if isempty(sysname)
    if isempty(eventdata) || eventdata~=1
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
    end
    return
end
sys=evalin('base',sysname);
if isfield(sys,'obtained_simulation_time')
    s=sprintf('%s %%',num2str(sys.obtained_simulation_time));
    set(handles.tx_an_simtimeobt,'String',s)
    if isempty(eventdata) || eventdata~=1
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
    end
    return
end
if ~isfield(sys,'originalsysname')
    errordlg('Original system unknown','Error Dialog','modal')
    uiwait
    if isempty(eventdata) || eventdata~=1
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
    end
    return
end
if ~existInBaseWs(sys.originalsysname)
    errordlg('Original system is not in the base workspace','Error Dialog','modal')
    uiwait
    if isempty(eventdata) || eventdata~=1
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
    end
    return
end
origsys=evalin('base',sys.originalsysname);
if ~isfield(origsys,'Simulationtime')
    try
        residue(origsys);
        assignin('base',sys.originalsysname,origsys)    
    catch ex
        if strfind(ex.identifier,'nomem')
            errordlg('Out of memory. Original system is too large.','Error Dialog','modal')
            uiwait
            if isempty(eventdata) || eventdata~=1
                set(handles.figure1,'Pointer','arrow')
                set(handles.virtgr_an_red_buttons,'Enable','on')
            end
            return
        elseif strfind(ex.message,'singular')
            errordlg(ex.message,'Error Dialog','modal')
            uiwait
            set(handles.figure1,'Pointer','arrow')
            set(hObject,'String','Plot')
            set(hObject,'Enable','on')
            return
        else
            set(handles.figure1,'Pointer','arrow')
            set(handles.virtgr_an_red_buttons,'Enable','on')
            throw(ex)
        end
    end
end
if ~isfield(sys,'Simulationtime')
    try
        residue(sys);
        assignin('base',sysname,sys)
    catch ex
        if strfind(ex.identifier,'nomem')
            errordlg('System is too big, out of memory','Error Dialog','modal')
            uiwait
            if isempty(eventdata) || eventdata~=1
                 set(handles.figure1,'Pointer','arrow') % *** ??
                 set(handles.virtgr_an_red_buttons,'Enable','on')
            end
            return
        elseif strfind(ex.message,'singular')
                errordlg('The system contains algebraic states.','Error Dialog','modal')
                uiwait
                set(hObject,'Enable','on')
                set(hObject,'String','Plot')
                set(handles.figure1,'Pointer','arrow')
                return
        else
            set(handles.figure1,'Pointer','arrow')
            set(handles.virtgr_an_red_buttons,'Enable','on')
            throw(ex)
        end
    end
end
set(handles.tx_an_simtime,'String',num2str(sys.Simulationtime))
sys.obtained_simulation_time=100-sys.Simulationtime*100/origsys.Simulationtime;
assignin('base',sysname,sys)
s=sprintf('%s %%',num2str(sys.obtained_simulation_time));
set(handles.tx_an_simtimeobt,'String',s)
if isempty(eventdata) || eventdata~=1
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
end


function pb_an_simtime_Callback(hObject, eventdata, handles)
set(handles.figure1,'Pointer','watch')
set(handles.virtgr_an_red_buttons,'Enable','off')
drawnow

% get and convert system
try
    [sys, sysname] = get_sys_from_ws(handles.pu_an_sys);
catch ex
    set(hObject,'String','Plot')
    set(hObject,'Enable','on') 
    errordlg(ex.message,'Error Dialog','modal')
    uiwait
    return
end

x = get(handles.pu_an_sys,'String');
sysname=x{get(handles.pu_an_sys,'Value')};
if isempty(sysname)
    if isempty(eventdata) || eventdata~=1
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
    end
    return
end
sys=evalin('base',sysname);
if ~ismpty(sys.SimulationTime)
    try
        residue(sys);
        assignin('base',sysname,sys)
    catch ex
        if strfind(ex.identifier,'nomem')
            errordlg('System is too big, out of memory','Error Dialog','modal')
            uiwait
            if isempty(eventdata) || eventdata~=1
                set(handles.figure1,'Pointer','arrow')  %*** ??
                set(handles.virtgr_an_red_buttons,'Enable','on')
            end
            return
        elseif strfind(ex.message,'singular')
            errordlg('The system contains algebraic states.','Error Dialog','modal')
            uiwait
            set(hObject,'Enable','on')
            set(hObject,'String','Plot')
            set(handles.figure1,'Pointer','arrow')
            return
        else
            set(handles.figure1,'Pointer','arrow')
            set(handles.virtgr_an_red_buttons,'Enable','on')
            errordlg(ex.message,'Error Dialog','modal')
        end
    end
end
set(handles.tx_an_simtime,'String',num2str(sys.SimulationTime))
if isempty(eventdata) || eventdata~=1
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
end


function pu_an_sys_Callback(hObject, eventdata, handles)
set(handles.figure1,'Pointer','watch')
set(handles.panel_an_sysred,'Visible','off')
set(handles.tx_an_h2,'String','')
set(handles.tx_an_simtime,'String','')
set(handles.tx_an_hinf,'String','')
set(handles.tx_an_stability,'String','')
drawnow

try
    sys = get_sys_from_ws(hObject);
catch ex
    set(handles.figure1,'Pointer','arrow')
    if strfind(ex.identifier, 'unassigned')
        return
    end
    errordlg(['Variable is not a valid state space model. ' ex.message],'Error Dialog','modal')
    set(handles.tx_an_sysinfo,'String','Invalid model')
    uiwait
    return
end

set(handles.tx_an_sysinfo, 'String', sys.disp);
set(handles.tx_an_sysinfored, 'String', sys.disp_mor_info);
%*** system_information(sys,handles.tx_an_sysinfo,handles.tx_an_sysinfored)

if ~isempty(sys.mor_info)
    set(handles.panel_an_sysred,'Visible','on')
    if isfield(sys,'H_2_error_rel')
        s=sprintf('%s %%',num2str(sys.H_2_error_rel));
        set(handles.tx_an_h2errorrel,'String',s)
    else
        set(handles.tx_an_h2errorrel,'String','')
    end
    if isfield(sys,'H_2_error_abs')
        set(handles.tx_an_h2errorabs,'String',sys.H_2_error_abs)
    else
        set(handles.tx_an_h2errorabs,'String','')
    end
    if isfield(sys,'H_inf_error_rel')
        s=sprintf('%s %%',num2str(sys.H_inf_error_rel));
        set(handles.tx_an_hinferrorrel,'String',s)
    else
        set(handles.tx_an_hinferrorrel,'String','')
    end
    if isfield(sys,'H_inf_error_abs')
        set(handles.tx_an_hinferrorabs,'String',sys.H_inf_error_abs)
    else
        set(handles.tx_an_hinferrorabs,'String','')
    end
    if isfield(sys,'obtained_simulation_time')
        set(handles.tx_an_simtimeobt,'String',num2str(sys.obtained_simulation_time))
    else
        set(handles.tx_an_simtimeobt,'String','')
    end
    if isfield(sys,'Simulationtime_original_system')
        set(handles.tx_an_simtimeorig,'String',num2str(sys.Simulationtime_original_system))
    else
        set(handles.tx_an_simtimeorig,'String','')
    end
else
    set(handles.panel_an_sysred,'Visible','off')
end
if isfield(sys, 'H_2_norm')
    set(handles.tx_an_h2,'String', num2str(sys.H_2_norm))
else
    set(handles.tx_an_h2,'String','')
end
if isfield(sys,'Stability')
    set(handles.tx_an_stability,'String',sys.Stability)
else
    set(handles.tx_an_stability,'String','')
end
if isfield(sys, 'H_inf_norm')
    set(handles.tx_an_hinf,'String',num2str(sys.H_inf_norm))
else
    set(handles.tx_an_hinf,'String','')
end
if isfield(sys,'Simulationtime')
    set(handles.tx_an_simtime,'String',num2str(sys.Simulationtime))
else
    set(handles.tx_an_simtime,'String','')
end
set(handles.figure1,'Pointer','arrow')

function pb_an_simtimeorig_Callback(hObject, eventdata, handles)
set(handles.figure1,'Pointer','watch')
set(handles.virtgr_an_red_buttons,'Enable','off')
drawnow
x = get(handles.pu_an_sys,'String');
sysname=x{get(handles.pu_an_sys,'Value')};
if isempty(sysname)
    if isempty(eventdata) || eventdata~=1
        set(handles.figure1,'Pointer','arrow')
        set(handles.virtgr_an_red_buttons,'Enable','on')
    end
    return
end
sys=evalin('base',sysname);
if ~isfield(sys,'Simulationtime_original_system')
    if ~isfield(sys,'originalsysname')
        errordlg('Original system unknown','Error Dialog','modal')
        uiwait
        if isempty(eventdata) || eventdata~=1
            set(handles.figure1,'Pointer','arrow')
            set(handles.virtgr_an_red_buttons,'Enable','on')
        end
        return
    end
    if ~existInBaseWs(sys.originalsysname)
        errordlg('Original system is not in the base workspace','Error Dialog','modal')
        uiwait
        if isempty(eventdata) || eventdata~=1
            set(handles.figure1,'Pointer','arrow')
            set(handles.virtgr_an_red_buttons,'Enable','on')
        end
        return
    end
    origsys=evalin('base',sys.originalsysname);
    if ~isfield(origsys,'Simulationtime')
        try
            residuen(origsys);
            assignin('base',sys.originalsysname,origsys)
        catch ex
            if strfind(ex.identifier,'nomem')
                errordlg('Original system is too big, out of memory','Error Dialog','modal')
                uiwait
                if isempty(eventdata) || eventdata~=1
                    set(handles.figure1,'Pointer','arrow')  % *** ??
                    set(handles.virtgr_an_red_buttons,'Enable','on')
                end
                return
            elseif strfind(ex.message,'singular')
                errordlg('The system contains algebraic states.','Error Dialog','modal')
                uiwait
                set(hObject,'Enable','on')
                set(hObject,'String','Plot')
                set(handles.figure1,'Pointer','arrow')
                return
            else
                set(handles.figure1,'Pointer','arrow')
                set(handles.virtgr_an_red_buttons,'Enable','on')
                throw(ex)
            end
        end
    end
    sys.Simulationtime_original_system=origsys.Simulationtime;
    assignin('base',sysname,sys)
end
set(handles.tx_an_simtimeorig,'String',sys.Simulationtime_original_system)
if isempty(eventdata) || eventdata~=1
    set(handles.figure1,'Pointer','arrow')
    set(handles.virtgr_an_red_buttons,'Enable','on')
end

function pb_an_calcall_Callback(hObject, eventdata, handles)
% calculate everything
set(handles.figure1,'Pointer','watch')
drawnow
x = get(handles.pu_an_sys,'String');
sysname=x{get(handles.pu_an_sys,'Value')};
if isempty(sysname)
    set(handles.figure1,'Pointer','arrow')
    errordlg('Please choose a system first!','Error Dialog','modal')
    uiwait
    return
end
sys=evalin('base',sysname);
% show warning if system order is above 500
if size(sys.A,1)>500
    s=sprintf('%s is a system with high order, this may take some time.\nDo you want to continue?',sysname);
    k=stqd('String',s,'Title','Question Dialog');
    if ~isempty(k) && strcmp(k,'No')
        set(handles.figure1,'Pointer','arrow')
        return
    end
end
set(handles.virtgr_an_red_buttons,'Enable','off')
if strcmp(get(handles.panel_an_sysred,'Visible'),'on')
    % reduced system
    pb_an_simtimeorig_Callback(handles.pb_an_simtimeorig, 1, handles)
    pb_an_simtimeobt_Callback(handles.pb_an_simtimeobt, 1, handles)
    pb_an_h2error_Callback(handles.pb_an_h2error, 1, handles)
    pb_an_hinferror_Callback(handles.pb_an_hinferror, 1, handles)
end
try
    pb_an_stability_Callback(handles.pb_an_stability,1, handles)
    pb_an_h2_Callback(handles.pb_an_h2,1, handles)
    pb_an_hinf_Callback(handles.pb_an_hinf, 1, handles)
    pb_an_simtime_Callback(handles.pb_an_simtime,1, handles)
catch ex
    assignin('base','lastError',ex)
    errordlg(ex.message,'Error Dialog','modal')
    uiwait
end
set(handles.virtgr_an_red_buttons,'Enable','on')
set(handles.figure1,'Pointer','arrow')







%--------------------------------------------------------------------------
%                               FOOTER
%--------------------------------------------------------------------------

%Callbacks of the panel

function logos_footer_CreateFcn(hObject, eventdata, handles)
    %Load the bitmap with the logos

    A=imread('Pictures\Footer.png');
    h=image(A);
    set(h,'HitTest','off');


function logos_footer_ButtonDownFcn(hObject, eventdata, handles)
    %Web-links to the diverent homepages (Open if the user klicks on a logo)

    p=get(hObject,'CurrentPoint');

    if p(1,1)>915 && p(1,1)<1044
        web www.tum.de
    elseif p(1,1)>5 && p(1,1)<170
        web www.rt.mw.tum.de/forschung/forschungsgebiete/modellreduktion/sssmor/
    elseif p(1,1)>580 && p(1,1)<850 
        web www.rt.mw.tum.de
    elseif p(1,1)>250 && p(1,1)<480
        web www.rt.mw.tum.de/forschung/forschungsgebiete/modellreduktion/
    end
    
    
    
    
    
%--------------------------------------------------------------------------
%                           UTILITY FUNCTIONS
%--------------------------------------------------------------------------

function list_open_figures(handles)
%global fighand
val=get(handles.figure,'Value');
if val>1 
    % UserData contains handles to figures
    userdata=get(handles.figure,'UserData');
    % save handle to selected figure
    k=userdata(val);
end
% get handles of all open figures. their order changes each time!
% figures are focussed and thereby can overlap the sssMOR_GUI

openfig=get(0,'Children');
c=cell({'New Figure'});
try
    h=0;
    in=1; % for SISO
    out=1; % for SISO
    if strcmp(get(handles.panel_intoout,'Visible'),'on')
        % MIMO system selected
        x=get(handles.pu_in,'String');
        in=str2double(x{get(handles.pu_in,'Value')}); % chosen input
        x=get(handles.pu_out,'String');
        out=str2double(x{get(handles.pu_out,'Value')}); % chosen output
        x=get(handles.syschoice,'String');
        sys=x{get(handles.syschoice,'Value')};
        % remove figure with wrong number of subplots
        if isnan(in)
            % 'all' has been selected
            in=evalin('base',sprintf('size(%s.B,2)',sys));
        else
            % number of input does not matter, only number of subplots
            in=1; 
        end
        if isnan(out)
            out=evalin('base',sprintf('size(%s.C,1)',sys));
        else
            out=1;
        end
    end
    switch get(handles.plot_type,'Value')
        case 1 %Impulse Response
            tag='i';
        case 2 % Step Response
            tag='s';
        case 3
            % Bode Diagram, requires two subplots per channel
            tag='b';
            out=2*out; 
        case 4 % Pole-Zero map
            tag='p';
        case 5 % Frequency Response
            tag='f';
    end
    
    for i=1:length(openfig)
        if openfig(i)==handles.figure1
            % sssMOR_GUI
            continue
        elseif strcmp(get(openfig(i),'Tag'),'composeModel')
            continue
        end
        temp=get(openfig(i),'Name');
         if ~isempty(strfind(get(openfig(i),'Tag'),tag))
             % include only figures of same type
             if in*out>1
                 % right number of subplots?
                 % store handles to subplots in UserData
                 x=get(openfig(i),'UserData');
                 if size(x,1)~=out || size(x,2)~=in
                     continue
                 end
             end
             if isempty(temp)
                 % figure has no name, so its number is used instead 
                temp=sprintf('Figure %i',openfig(i));
             end
             
             h=[h, openfig(i)]; %#ok<AGROW>
             c=[c, {temp}]; %#ok<AGROW> % names
         end
    end

    if size(c, 2) < val % number of entries
        % selected figure number is above number of entries, use 1 instead
        set(handles.figure,'Value',1)
    end
    set(handles.figure,'String', c); % set name
    set(handles.figure,'UserData', h); % save handles in UserData
    % if figure was selected originally, its handle is k. reselect it
    if exist('k','var')==1
        for i=1:length(h)
            if h(i)==k
                set(handles.figure,'Value',i)
            end
        end
    else
        set(handles.figure,'Value',1)
    end
    % bring sssMOR_GUI back to front
    figure(handles.figure1)
catch %#ok<CTCH>
    % if an error occurs, restore default setting
    set(handles.figure,'Value',1) 
    set(handles.figure,'String','New Figure');
    set(handles.figure,'UserData',0);
end


function [s0,sOut]=getExpansionPoints(handles)
% read expansion point(s) from GUI

    %Read the expension points from the first table

    if get(handles.uitable_mor_krylov,'UserData')==1
        errordlg('Please correct expansion points first','Error Dialog','modal')
        uiwait
        s0=[];
        return
    end
    s0=cell2mat(get(handles.uitable_mor_krylov,'Data'));% Entwicklungspunkte
    if isempty(s0) || size(s0,2)==1 % Tabelle leer, oder nur Momente ohne Punkt
        errordlg('Please choose expansion points and number of matching moments','Error Dialog','modal')
        uiwait
        s0=[];
        return
    end
    x = get(handles.pu_mor_systems,'String');
    y=x{get(handles.pu_mor_systems,'Value')};
    s=evalin('base',sprintf('size(%s.A,1);',y)); %Größe Originalsystem
    if s<sum(s0(:,2))
       errordlg('Number of matching moments must be smaller than dimension of original system!','Error Dialog','modal')
       uiwait
       s0=[];
       return
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
        
        if get(handles.rb_mor_krylov_hermite,'Value') == 1
           
            sIn = s0;
            sOut = s0;
            
        else
            
            %Read out the expension points from the output table
            
            sIn = s0;
            
            if get(handles.uitable_mor_krylov_output,'UserData')==1
                errordlg('Please correct expansion points first','Error Dialog','modal')
                uiwait
                sOut=[];
                return
            end
            sOut=cell2mat(get(handles.uitable_mor_krylov_output,'Data'));% Entwicklungspunkte
            if isempty(sOut) || size(sOut,2)==1 % Tabelle leer, oder nur Momente ohne Punkt
                errordlg('Please choose expansion points and number of matching moments','Error Dialog','modal')
                uiwait
                sOut=[];
                return
            end
            if sum(sIn(2,:)) ~= sum(sOut(:,2))
               errordlg('Number of matching moments must be equal for input and output','Error Dialog','modal')
               uiwait
               sOut=[];
               return
            end 
            sOut=transpose(sOut);
            
        end       
    end

function test_colors(handles,hObject)
% check inserted colors
h=str2num(get(hObject,'String')); %#ok<ST2NM>
if isempty(h)
    errordlg('R,G,B must be elements of the interval [0;1]','Error Dialog','modal')
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
        if str2double(h)>1
            errordlg('R,G,B must not be greater than one','Error Dialog','modal')
            uiwait
            set(hObject,'UserData',1)
            return
        end
        set(hObject,'UserData',0)
        return
    end
    errordlg('Use ''.'' as decimal seperator','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
elseif h<0
    errordlg('R,G,B must not be less than zero','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
elseif h>1
    errordlg('R,G,B must not be greater than one','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
else
    set(hObject,'UserData','0')
end


function test_width(hObject)
% check inserted linewidth

h=str2num(get(hObject,'String')); %#ok<ST2NM>
if isempty(h)
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
    set(hObject,'UserData','0')
end


function test_max(hObject,handles)
% check maximum frequency/time in manual mode
h=str2num(get(hObject,'String')); %#ok<ST2NM>
if isempty(h)
    errordlg('Max must be greater than Min. Use ''.'' as decimal seperator. Don''t use characters','Error Dialog','modal')
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
    errordlg('Max must not be less than zero','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
elseif h<= str2num(get(handles.ed_min,'String')) %#ok<ST2NM>
    errordlg('Max must be greater than Min','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
else
    set(hObject,'UserData','0')
end

function test_min(hObject,handles)
% check minimum frequency/time in manual mode
h=str2num(get(hObject,'String')); %#ok<ST2NM>
if isempty(h)
    errordlg('The minimal frequency must be greater than zero. Use ''.'' as decimal seperator. Don''t use characters','Error Dialog','modal')
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
    errordlg('The minimal frequency must not be less than zero','Error Dialog','modal')
    uiwait
    set(hObject,'UserData',1)
else
    set(hObject,'UserData','0')
end


function x=listS0InWorkspace
% list all vectors from workspace that might be expansion points
s=evalin('base', 'whos');
% preallocate memory
x=cell(length(s),1);
for i=1:length(s)
    if strcmp(s(i).class,'double') && length(s(i).size)==2 && (s(i).size(1)==2 || s(i).size(2)==2)
        %check wheather one dimension is integer
        vec = evalin('base',s(i).name);
        isInt1 = 1;
        isInt2 = 1;
        if size(vec,1) == 2 && size(vec,2) == 2
            isInt1 = 0;
            isInt2 = 0;
            if (mod(vec(1,1),1)==0 && mod(vec(2,1),1)==0) || ...
               (mod(vec(1,2),1)==0 && mod(vec(2,2),1)==0)
           
                isInt1 = 1;
            end
            if (mod(vec(1,1),1)==0 && mod(vec(1,2),1)==0) || ...
               (mod(vec(2,1),1)==0 && mod(vec(2,2),1)==0)
           
                isInt2 = 1;
            end
     
        elseif size(vec,1) == 2
            for j = 1:size(vec,2)
               if mod(vec(1,j),1)~=0
                  isInt1 = 0; 
               end
               if mod(vec(2,j),1)~=0
                  isInt2 = 0; 
               end
            end
            
        else
            for j = 1:size(vec,1)
               if mod(vec(j,1),1)~=0
                  isInt1 = 0; 
               end
               if mod(vec(j,2),1)~=0
                  isInt2 = 0; 
               end
            end
        end
        
        % save name
        
        if isInt1 == 1 || isInt2 == 1
            x{i}=s(i).name;
        end
    end
end
% remove empty (non-system) entries
x(cellfun(@isempty,x)) = [];
if isempty(x)
   x = {''}; 
end



function [sys, sysname] = get_sys_from_ws(namehandle)
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

% convert to sss
if ~strcmp(class(sys), 'sss')
    sys=sss(sys);
end



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
        
    else
        
        method = 'krylov';
        
    end
    
    name = strcat(name,'_',method);
    
    %Set the names to the text-fields of the GUI-objects
    
    suggestVarname(name,handles.ed_mor_sysred);
    suggestVarname(sprintf('%s_w',name),handles.ed_mor_w);
    suggestVarname(sprintf('%s_v',name),handles.ed_mor_v);
        
function x = systemsInWorkspace()
% finds and lists all dynamical systems that are contained in workspace
    % read all variable names
    s=evalin('base', 'whos');
    % preallocate memory
    x=cell(length(s),1);
    for i=1:length(s)
        if strcmp(s(i).class,'ss') || strcmp(s(i).class,'sss') && ...
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

