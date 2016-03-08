function lineH = analyzeMor(sys,sysr,Opts)
%   This function generates several plots that compare the original, high
%   fidelity model "sys" to the reduced order model sysr
%
%   The function is programmed to accept sparse state-space (sss) models
%   only. An extension to "dense" state-space models (ss) can be easily
%   derived.

%%  Parse the input
Def.mode = 'red'; %reduced analysis ('full');
Def.lineH = struct('bmsys',[],'bmsysr',[],'bpsys',[],'bpsysr',[],'bmsyse',[],...
        'bpsyse',[],'psys',[],'psysr',[]); %empty structure to store the line handles
    
if ~isfield('Opts','lineH') || isempty(Opts.lineH)
    newplot = 1;
else
    newplot = 0;
end

if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
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
    
    
    mag = mag2db(mag); magr = mag2db(magr);
    switch Opts.mode
        case 'full'
            ax(1) = subplot(2,3,1);
        case 'red'
            ax(1) = subplot(2,1,1);
    end
            % amplitude
            Opts.lineH.bmsys = plot(omega,mag,'Color',TUM_Blau); hold on
            Opts.lineH.bmsysr = plot(omegar,magr,'Color',TUM_Orange);

            set(gca, 'XScale', 'log');
            set(gca, 'XLim', [min(omega) max(omega)]);
            box on
            mx=max(mag); mn=min(mag);
            set(gca, 'YLim', [mn-(mx-mn)/20, mx+(mx-mn)/20]);
            ylabel('Magnitude [dB]') 

            title('Bode plot');
            legend('HFM','ROM','Location','NorthEast');
    switch Opts.mode
        case 'full'
            ax(2) = subplot(2,3,4);  
        case 'red'
            ax(2) = subplot(2,1,2);
    end       
            % phase
            Opts.lineH.bpsys = plot(omega,phase,'Color',TUM_Blau);hold on
            Opts.lineH.bpsysr = plot(omegar,phaser,'Color',TUM_Orange);
            set(gca, 'XScale', 'log');
            set(gca, 'XLim', [min(omega) max(omega)]);
            ylabel('Phase [deg]');
            xlabel('Frequency [rad/sec]');

    % linkaxes(ax(1:2),'x');
    clear mag phase omega
    if strcmp(Opts.mode,'full')
        %%  Bode diagram of error system
        syse = error_system(sys,sysr);
        [mag, phase, omega] = bode(ss(syse));

        if iscell(mag)
            %replace code her
        else
            mag = {squeeze(mag)};
            phase = {squeeze(phase)};
        end

        ax(3) = subplot(2,3,2);
                % amplitude
                Opts.lineH.bmsyse = plot(omega,mag{1},'Color',TUM_Rot);

                set(gca, 'XScale', 'log','YScale','log');
                set(gca, 'XLim', [min(omega) max(omega)]);
                box on
                mx=max(mag{1}); mn=min(mag{1});
                title('Error system');
                ylabel('Magnitude [dB]') 

        ax(4) = subplot(2,3,5);        
                % phase
                Opts.lineH.bpsyse = plot(omega,phase{1},'Color',TUM_Rot);
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

            Opts.lineH.psys = plot(real(sys_poles),imag(sys_poles),'*','Color',TUM_Blau);hold on
            Opts.lineH.psysr = plot(real(sysr_poles),imag(sysr_poles),'o','Color',TUM_Orange);
            title('Pole configuration');

    %      Opts.lineH.ax = ax;

    %     %% Compute the H2 and Hinf norm of the error
    %     H2 = norm(sys-sysr);
    %     Hinf = norm(sys-sysr,'inf');
    %     
    %     ax(6) =subplot(2,3,6);
    %         Opts.lineH.h2 = plot(H2,'d','Color',TUM_Rot);hold on
    %         Opts.lineH.hinf = plot(Hinf,'p','Color',TUM_Rot);
    %     legend('H_2','H_{inf}');
    end
else
        %%  Bode diagram of sys and sysr
    [mag, phase, omega]     = get_bode(sys);
    [magr, phaser, omegar]  = get_bode(sysr);

            % amplitude
            set(Opts.lineH.bmsys,'XData',omega,'YData',mag{1});
            set(Opts.lineH.bmsysr,'XData',omegar,'YData',magr{1});
      
            % phase
            set(Opts.lineH.bpsys,'XData',omega,'YData',phase{1});
            set(Opts.lineH.bpsysr,'XData',omegar,'YData',phaser{1});

    clear mag phase omega

    if strcmp(Opts.mode, 'full')
        %%  Bode diagram of error system

        [mag, phase, omega]     = bode(sys-sysr);

                % amplitude
                set(Opts.lineH.bmsyse,'XData',omega,'YData',mag{1});

                % phase
                set(Opts.lineH.bpsyse,'XData',omega,'YData',phase{1});

        clear mag phase omega
        %% Compare the poles

            sys_poles = get_poles(sys);
            sysr_poles = get_poles(sysr);

            set(Opts.lineH.psys,'XData',real(sys_poles),'YData',imag(sys_poles));
            set(Opts.lineH.psysr,'XData',real(sysr_poles),'YData',imag(sysr_poles));

    %     %% Compute the H2 and Hinf norm of the error
    %         H2 = norm(sys-sysr);
    %         Hinf = norm(sys-sysr,'inf');
    %         xdat = [get(Opts.lineH.h2,'XData'),legth(get(Opts.lineH.h2,'XData'))+1];
    %     
    %         set(Opts.lineH.h2,'XData',xdat,...
    %             'YData',...
    %             [get(Opts.lineH.h2,'YData'),H2]);
    %         set(Opts.lineH.hinf,'XData',xdat,...
    %             'YData',...
    %             [get(Opts.lineH.hinf,'YData'),Hinf]);
    end
end
%   Generate output
lineH = Opts.lineH;
    
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