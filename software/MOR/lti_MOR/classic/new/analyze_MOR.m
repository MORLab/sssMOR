function h = analyze_MOR(sys,sysr,h)

%   This function generates several plots that compare the original, high
%   fidelity model "sys" to the reduced order model sysr
%
%   The function is programmed to accept sparse state-space (sss) models
%   only. An extension to "dense" state-space models (ss) can be easily
%   derived.

%%  Parse the input
if ~exist('h','var') || isempty(h)
    % Generate the line handles
    h = struct('bmsys',[],'bmsysr',[],'bpsys',[],'bpsysr',[],'bmsyse',[],...
        'bpsyse',[],'psys',[],'psysr',[]); %empty structure to store the line handles
    
    newplot = 1;
else
    newplot = 0;
end

if newplot
    %% Generate the figure window
    
    %   get the new figure on the second monitor (if available) in full size
    if strcmp(getenv('USERNAME'),'castagnotto') && strcmp(getenv('computername'),'RT44');
        pos = [1281 1 1280 948]; %correct coordinates for castagnotto
    else
        pos = get(0,'MonitorPositions');
        sz = size(pos);
        if (sz(1) > 1)
            pos = pos(2,:);
        end
    end

    
    fh = figure('Name','Analyze MOR','Color','white','Position',pos);

    ax = zeros(6,1); %preallocating memory for the axis handle
    %%  Bode diagram of sys and sysr
    [mag, phase, omega] = get_bode(sys);
    [magr, phaser, omegar]  = get_bode(sysr);
    
    ax(1) = subplot(2,3,1);
            % amplitude
            h.bmsys = plot(omega,mag,'Color',TUM_Blau); hold on
            h.bmsysr = plot(omegar,magr,'Color',TUM_Orange);

            set(gca, 'XScale', 'log');
            set(gca, 'XLim', [min(omega) max(omega)]);
            box on
            mx=max(mag); mn=min(mag);
            set(gca, 'YLim', [mn-(mx-mn)/20, mx+(mx-mn)/20]);
            ylabel('Magnitude [dB]') 

            title('Bode plot');
            legend('HFM','ROM','Location','NorthEast');

    ax(2) = subplot(2,3,4);        
            % phase
            h.bpsys = plot(omega,phase,'Color',TUM_Blau);hold on
            h.bpsysr = plot(omegar,phaser,'Color',TUM_Orange);
            set(gca, 'XScale', 'log');
            set(gca, 'XLim', [min(omega) max(omega)]);
            ylabel('Phase [deg]');
            xlabel('Frequency [rad/sec]');

    % linkaxes(ax(1:2),'x');
    clear mag phase omega

    %%  Bode diagram of error system
    syse = error_system(sys,sysr);
    [mag, phase, omega]     = bode(ss(syse));
   
    if iscell(mag)
        %replace code her
    else
        mag = {squeeze(mag)};
        phase = {squeeze(phase)};
    end
        
    
    ax(3) = subplot(2,3,2);
            % amplitude
            h.bmsyse = plot(omega,mag{1},'Color',TUM_Rot);

            set(gca, 'XScale', 'log');
            set(gca, 'XLim', [min(omega) max(omega)]);
            box on
            mx=max(mag{1}); mn=min(mag{1});
            set(gca, 'YLim', [mn-(mx-mn)/20, mx+(mx-mn)/20]);
            title('Error system');
            ylabel('Magnitude [dB]') 

    ax(4) = subplot(2,3,5);        
            % phase
            h.bpsyse = plot(omega,phase{1},'Color',TUM_Rot);
            set(gca, 'XScale', 'log');
            set(gca, 'XLim', [min(omega) max(omega)]);
            ylabel('Phase [deg]');
            xlabel('Frequency [rad/sec]');

    linkaxes(ax(1:4),'x');
    clear mag phase omega
    %% Compare the poles

    ax(5) = subplot(2,3,3);
        sys_poles = get_poles(sys);
        sysr_poles = get_poles(sysr);

        h.psys = plot(real(sys_poles),imag(sys_poles),'*','Color',TUM_Blau);hold on
        h.psysr = plot(real(sysr_poles),imag(sysr_poles),'o','Color',TUM_Orange);
        title('Pole configuration');
    
%      h.ax = ax;

%     %% Compute the H2 and Hinf norm of the error
%     H2 = norm(sys-sysr);
%     Hinf = norm(sys-sysr,'inf');
%     
%     ax(6) =subplot(2,3,6);
%         h.h2 = plot(H2,'d','Color',TUM_Rot);hold on
%         h.hinf = plot(Hinf,'p','Color',TUM_Rot);
%     legend('H_2','H_{inf}');
else
        %%  Bode diagram of sys and sysr
    [mag, phase, omega]     = get_bode(sys);
    [magr, phaser, omegar]  = get_bode(sysr);

            % amplitude
            set(h.bmsys,'XData',omega,'YData',mag{1});
            set(h.bmsysr,'XData',omegar,'YData',magr{1});
      
            % phase
            set(h.bpsys,'XData',omega,'YData',phase{1});
            set(h.bpsysr,'XData',omegar,'YData',phaser{1});

    clear mag phase omega

    %%  Bode diagram of error system

    [mag, phase, omega]     = bode(sys-sysr);

            % amplitude
            set(h.bmsyse,'XData',omega,'YData',mag{1});
      
            % phase
            set(h.bpsyse,'XData',omega,'YData',phase{1});

    clear mag phase omega
    %% Compare the poles

        sys_poles = get_poles(sys);
        sysr_poles = get_poles(sysr);
        
        set(h.psys,'XData',real(sys_poles),'YData',imag(sys_poles));
        set(h.psysr,'XData',real(sysr_poles),'YData',imag(sysr_poles));
        
%     %% Compute the H2 and Hinf norm of the error
%         H2 = norm(sys-sysr);
%         Hinf = norm(sys-sysr,'inf');
%         xdat = [get(h.h2,'XData'),legth(get(h.h2,'XData'))+1];
%     
%         set(h.h2,'XData',xdat,...
%             'YData',...
%             [get(h.h2,'YData'),H2]);
%         set(h.hinf,'XData',xdat,...
%             'YData',...
%             [get(h.hinf,'YData'),Hinf]);
end
    
function [mag, phase, omega] = get_bode(sys)
%   gets the relevant information if sys is a structure containing the
%   information, or else it computes it from scratch.

if isstruct(sys) && isfield(sys,'bode')
        mag = sys.bode.mag;
        phase = sys.bode.phase;
        omega = sys.bode.omega;
    else
        [mag, phase, omega]     = bode(sys);
        mag = mag{1};
        phase = phase{1};
end

function p = get_poles(sys)

%   gets the relevant information if sys is a structure containing the
%   information, or else it computes it from scratch.

if isstruct(sys) && isfield(sys,'eigs')
    p = [];
    fnames = fieldnames(sys.eigs);
    for ii=1:length(fnames)
        p = [p;sys.eigs.(fnames{ii})];
    end
else
    p = eig(sys);
end



function syse = error_system(sys,sysr)

if isstruct(sys)
    HFM = sys.sys;
else
    HFM = sys;
end
if isstruct(sysr)
    ROM = sysr.sysr;
else
    ROM = sysr;
end

syse = HFM-ROM;