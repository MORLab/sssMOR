
function  [mag, phase, omega] = bode(sys, varargin)
% Plots the bode diagram of an LTI system
% ------------------------------------------------------------------
% [mag, phase, omega] = bode(sys, omega, in, out, options)
% Inputs:       * sys: an sss-object containing the LTI system
%    [optional] * vector of imaginary frequencies to plot at
%               * plot options. see <a href="matlab:help plot">PLOT</a>
% Outputs:      * vector of complex frequency response values
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer, Rudy Eid,
%               Alessandro Castagnotto
% Last Change:  31 March 2015
% ------------------------------------------------------------------

options = {};
% in=[]; out=[];

% --------------- EVALUATE OPTIONS ---------------

if nargin>1
    options = varargin;
    if strcmp(class(varargin{1}), 'double')
        omega = varargin{1};
        options(1) = [];
    end
%     if length(options)>=2
%         if strcmp(class(options{1}), 'double') && ...
%            strcmp(class(options{2}), 'double')
%             % I/O
%             if ~isempty(options{1}) && ~isnan(options{1}) 
%                 %input number
%                 in=options{1};
%             end
%             if ~isempty(options{2}) && ~isnan(options{2})
%                 %output number
%                 out=options{2};
%             end
%             options(1:2) = [];
%         end
%     end
end

% --------------- CALCULATE FREQUENCY RESPONSE ---------------

if exist('omega', 'var') && ~isempty(omega)
    % --------- frequency values given ---------
    if size(omega,1)>1
        omega=transpose(omega);
    end
    m = freqresp(sys, 1j*abs(omega));
else
    % --------- frequency values need to be chosen ---------
    dc = freqresp(sys,0);    % G(0)=DCgain
    ft = freqresp(sys,inf);  % G(inf)=feedthrough
  
    %determine minimum frequency
    if any(any(cellfun(@isinf,dc))) || any(any(cellfun(@isnan,dc))) % pole at s=0
        wmin = 0;   %***
        dc = num2cell(ones(size(dc)));
    elseif any(any(cellfun(@abs,dc)<1e-14))   % transfer zero at s=0
        wmin = 0;   %***
        dc = num2cell(ones(size(dc)));
    else
        wmin=0; t = freqresp(sys, 10^wmin);
        while cellfun(@(x,y) norm(x-y)/norm(y),t,dc) > 1e-2
            wmin=wmin-1; t = freqresp(sys, 1i*10^wmin);
        end
        while cellfun(@(x,y) norm(x-y)/norm(y),t,dc) < 1e-2
            wmin=wmin+1; t = freqresp(sys, 1i*10^wmin);
        end
        wmin=wmin-1;
    end
    
    %determine maximum frequency
    wmax=0; t = freqresp(sys, 10^wmax);
    while cellfun(@(x,y,z) norm(x-y)/norm(z),t,ft,dc) > 1e-6
        wmax=wmax+1; t = freqresp(sys, 1i*10^wmax);
    end
    while cellfun(@(x,y,z) norm(x-y)/norm(z),t,ft,dc) < 1e-6
        wmax=wmax-1; t = freqresp(sys, 1i*10^wmax);
    end
    wmax=wmax+1;

    delta = (wmax-wmin)/19; % initial resolution (insert odd number only!)
    omega = 10.^(wmin:delta:wmax);
    m = freqresp(sys, 1i* omega);
    while(1)
        % increase plot density until vertical change per step is below 1%
        for i=1:length(omega)-1
            if cellfun(@(x) abs(abs(x(i)) - abs(x(i+1)))/(abs(x(i)) + abs(x(i+1))),m) > 0.01
                break
            end
        end
        if i==length(omega)-1
            break
        end
        % do not refine above 2000 points
        if length(omega)>1000
            break
        end
        delta = delta/2;
        omega = 10.^(wmin:delta:wmax);

        % calculate new values of frequency response
        temp=freqresp(sys, 1i*omega(2:2:length(omega)));
        % update array of results (insert new values)
        m=cellfun(@(x,y) [reshape([x(1:length(x)-1);y],1,2*length(x)-2),x(end)],m,temp,'UniformOutput',false);
    end
end

% determine magnitude and phase from complex frequency response
mag = cellfun(@abs, m, 'UniformOutput',false);
phase = cellfun(@(x) unwrap(angle(x))*180/pi,m, 'UniformOutput',false);

% determine H_inf-norm (maximum magnitude)
[a,b]=cellfun(@max, mag);
[a,c]=max(a);
[H_inf_norm,d]=max(a);
H_inf_peakfreq=omega(b(c(d),d));
sys.H_inf_norm = max([sys.H_inf_norm, H_inf_norm]);
sys.H_inf_peakfreq = max([sys.H_inf_peakfreq, H_inf_peakfreq]);
% store system in caller workspace
if inputname(1)
    assignin('caller', inputname(1), sys);
end

% --------------- PLOT ---------------
if nargout>0
    return
end

if isempty(options)
    options={'Color','b'};
end

mag=cellfun(@(x) 20*log10(x), mag, 'UniformOutput',false);
[m, p] = size(mag); % inputs, outputs
plot_handles=zeros(2*m,p);

for i=1:m %for each input
    for j=1:p %for each output
        plot_handles(2*i-1,j)=subplot(2*m,p,2*((i-1)*p)+j);
        hold on        
%         if j==1
%             y_lab=sprintf('To Out(%i)',ceil(i/2));
%             ylabel(y_lab)
%         end
        if i==1 && m>1
            x_lab=sprintf('From In(%i)',j);
            title(x_lab)
        end
      
        % amplitude
        y_plot=mag{i,j};
        plot(omega,y_plot,options{:})
        set(gca, 'XScale', 'log');
        set(gca, 'XLim', [min(omega) max(omega)]);
        box on
        mx=max(mag{i,j}); mn=min(mag{i,j});
        set(gca, 'YLim', [mn-(mx-mn)/20, mx+(mx-mn)/20]);
        ylabel('Magnitude [dB]') 

        plot_handles(2*i,j)=subplot(2*m,p,(2*i-1)*p+j);
        hold on
        box on
        % phase
        y_plot=phase{i,j};
        plot(omega,y_plot,options{:})
        set(gca, 'XScale', 'log');
        set(gca, 'XLim', [min(omega) max(omega)]);
        ylabel('Phase [deg]');
        xlabel('Frequency [rad/sec]');
    end
end

% figure and axes properties
set(gcf,'UserData',plot_handles)
set(gcf,'Color', [1 1 1])

% link the axes for simultaneous zooming
linkaxes(reshape(plot_handles,2*m*p,1),'x');

% h=axes('position',[0,0,1,1],'Visible','off');
% text(0.4,0.98,'Bode Diagram');
% text(0.5,0.02,'Frequency [rad/sec]')
% text(0.01,0.5,'Magnitude [dB] ; Phase [deg]','Rotation',90)
% set(h,'HandleVisibility','off')

% for i=1:size(plot_handles,1)
%     for j=1:size(plot_handles,2)
%         axes(plot_handles(i,j))
%     end
% end
% axes(plot_handles(1,1))

if nargout==0
    clear mag phase w
end

