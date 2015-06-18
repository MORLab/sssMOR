function [p, z] = pzmap(sys, varargin)
% Computes the poles and zeros of an LTI system
% ------------------------------------------------------------------
% [p, z] = pzmap(sys, varargin)
% Inputs:       * sys: an sss-object containing the LTI system
%    [optional] * plot options. see <a href="matlab:help plot">PLOT</a>
% Outputs:      * p, z: vectors containing poles and zeros
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer
% Last Change:  19 Jan 2012
% ------------------------------------------------------------------

% are poles already available?
if ~isempty(sys.poles)
    p=sys.poles;
else
    % no, they are not. solve for eigenvalues of system
    if sys.is_dae
        p = 1./eig(full(sys.E), full(sys.A));
    else
        p = eig(full(sys.A));
    end
    % ensure column vector
    if size(p,1)<size(p,2)
        p=transpose(p);
    end
    sys.poles=p;
    
    % store system in caller workspace
    if inputname(1)
        assignin('caller', inputname(1), sys);
    end
end

% are zeros already available?
if ~isempty(sys.invariant_zeros)
    z=sys.invariant_zeros;
else
    % no, they are not. solve for gen. eigenvalues of Rosenbrock matrix
    z = cell(sys.p,sys.m);
    for i=1:sys.m
        for j=1:sys.p
            z{j,i}=eig(full([sys.A,sys.B(:,i);sys.C(j,:),sys.D(j,i)]), ...
                       [full(sys.E),zeros(sys.n,1);zeros(1,sys.n),0]);
            % ensure column vector
            if size(z{j,i},1)<size(z{j,i},2)
                z{j,i}=transpose(z{j,i});
            end
        end
    end

    % remove zeros at infinity
    z=cellfun(@(x) x(~isinf(x)), z, 'UniformOutput', false);
    sys.invariant_zeros=z;
    
    % store system in caller workspace
    if inputname(1)
        assignin('caller', inputname(1), sys);
    end
end

if nargout>0
    return
end

% --------------- PLOT ---------------

options=varargin;
fig_handle=gcf;
axes_handle=zeros(sys.p,sys.m);

% set random color if figure is not empty
if isempty(options)
    if ~isempty(get(fig_handle, 'Children'))
        c=rand(3,1); c=c/norm(c);
        options = {'Color', c};
    end
end

for i=1:sys.m
    for j=1:sys.p
        
        axes_handle(i,j)=subplot(sys.p,sys.m,i*sys.p+j-1);
        box on
        
        % plot o for zeros
        plot_handle=plot(real(z{j,i}), imag(z{j,i}), options{:});
        set(plot_handle, 'Marker', 'o', 'LineStyle', 'none');
        hold on
        % plot x for poles, remove legend entry
        plot_handle=plot(real(p), imag(p), options{:});
        set(plot_handle, 'Marker', 'x', 'LineStyle', 'none');
        hAnnotation = get(plot_handle,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')

        % determine x and y boundary
        mni=min(imag([p;z{j,i}])); mxi=max(imag([p;z{j,i}]));
        if mni*mxi<0 
            limy = [mni-(mxi-mni)/20 mxi+(mxi-mni)/20];
        elseif mni*mxi>0 
            limy = sort([0 1.05*max(abs([mni mxi]))]*sign(mni));
        else
            limy = [-1 1];
        end
        mnr=min(real([p;z{j,i}])); mxr=max(real([p;z{j,i}]));
        if mnr*mxr<0 
            limx = [mnr-(mxr-mnr)/20 mxr+(mxr-mnr)/20];
        elseif mnr*mxr>0
            limx = sort([0 1.05*max(abs([mnr mxr]))]*sign(mnr));
        else
            limx = [-1 1];
        end
        
        % plot dashed axes through origin, remove legend entry
        plot_handle=plot([0 0],limy,':k');
        hAnnotation = get(plot_handle,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')
        plot_handle=plot(limx,[0 0],':k');
        hAnnotation = get(plot_handle,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')
        set(axes_handle(i,j), 'XLim', limx, 'YLim', limy);

        % label input / output number
        if j==1 && sys.p>1
            y_lab=sprintf('To Out(%i)',i);
            ylabel(y_lab)
        end
        if i==1  && sys.m>1
            x_lab=sprintf('From In(%i)',j);
            title(x_lab)
        end
    end
end

% create invisible background plot for labelling
% h=axes('position',[0,0,1,1],'Visible','off');
% text(0.4,0.98,'Pole-Zero Map');
% text(0.5,0.02,'Real Axis')
% text(0.01,0.5,'Imaginary Axis','Rotation',90)
% set(h,'HandleVisibility','off')
% hold on
% % bring subplots to front
% for i=1:size(axes_handle,1)
%     for j=1:size(axes_handle,2)
% %         set(fig_handle, 'CurrentAxes', axes_handle(i,j));
%         axes(axes_handle(i,j))
%     end
% end
% set(fig_handle, 'CurrentAxes', axes_handle(1,1));


% % make subplots content of the figure
set(fig_handle,'UserData',axes_handle)

% avoid output
clear p z
