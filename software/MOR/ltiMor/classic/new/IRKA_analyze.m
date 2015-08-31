function [sysr, conv, V, W, s0] = IRKA_analyze(sys, s0, maxiter, epsilon,opts) 
% IRKA Algorithm for Model Order Reduction by Krylov Subspace Methods
% ------------------------------------------------------------------
% [sysr, V, W, s0] = IRKA(sys, [s0;q], maxiter, epsilon)
% Inputs:       * sys: an sss-object containing the LTI system
%               * s0: initial expansion point
%               * maxiter: maximum number of iterations
%               * epsilon: maximum step size for convergence
% Outputs:      * sysr: reduced system
%               * V, W: Projection matrices spanning Krylov subspaces
%               * s0: vector containing optimal expansion points
% ------------------------------------------------------------------
% For further explanation of the notation of s0, please refer to
% <a href="matlab:help RK">help RK</a>.
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de)
% Last Change:  23 Jan 2012
% ------------------------------------------------------------------

if ~exist('maxiter', 'var') || isempty(maxiter)
    maxiter=50;
end
if ~exist('epsilon', 'var') || isempty(epsilon)
    epsilon=1e-3;
elseif epsilon<=0 || ~isreal(epsilon)
    error('epsilon must be a real positive number.');
end

s0 = s0_vect(s0);
tic

sysr=sss([],[],[]);

s0_traj = zeros(maxiter+2, length(s0));
s0_traj(1,:) = s0;

    % new code: store convergence criterion
    es0_traj = zeros(maxiter, 1);


k=0;

% -------new analysis code-----------
if ~exist('opts','var') || isempty(opts)
    analysis_type = '';
    % analysis_type = 'red';
%     analysis_type = 'complete';
else
    analysis_type = opts{1};
    tpause = opts{2};
    
    [lh,th,ax] = initialize_plots(k,analysis_type,sys,s0);
    h = [];
    custompause(tpause);
end

%------------------------------
conv = 1; %initialize the convergence flag to 1
while true
    k=k+1;
    
    sysr_old = sysr;
    [sysr, V, W] = RK(sys, s0, s0, @(x,y) (x'*y));

    s0_old=s0;
    s0 = -eig(sysr)';
    curr_poles = -s0;

    s0(isnan(s0)) = 0;
    s0 = cplxpair(s0);

%     % mirror shifts with negative real part
    s0 = s0.*sign(real(s0)); %[Ale] how does it affect the convergence?

    s0_traj(k+1,:) = s0;
    
    %% --------new analysis code----------------
    if ~isempty(analysis_type)
        [h,lh] = plot_current(k,analysis_type,h,lh,th,ax,s0,s0_old,curr_poles,sys,sysr);
        custompause(tpause);
    end
%-----------------------------------   


es0 = stopping_crit(s0,s0_old); %compute the measure of convergence

    %   new code
    es0_traj(k) = es0;


    e=inf;
    if all(real(s0)>0)
        warning off
        if isreal(norm(sysr))
            e=norm(sysr-sysr_old)/norm(sysr);
        end
        warning on
    end
    
    if sys.n>500
        disp(['IRKA step ' num2str(k,'%03u') 9 ': Convergence '  num2str(es0, '%3.1e'), ', ' num2str(e, '%3.1e')]);
    end
    
    if es0 < epsilon || k >= maxiter %  || e < epsilon
        s0 = s0_old; % function return value
        s0_traj = s0_traj(1:(k+1),:);
        break
    end
end

assignin('base', 'IRKA_s0', s0_traj);

%new code
assignin('base', 'IRKA_es0', es0_traj);

if sys.n<=500
    disp(['IRKA step ' num2str(k,'%03u') 9 ': Convergence '  num2str(es0, '%3.1e'), ', ' num2str(e, '%3.1e')]);
end

disp(['IRKA exit. Elapsed time: ' num2str(toc) 's']);

if k>=maxiter
    warning('IRKA:no_converged', ['IRKA has not converged after ' num2str(k) ' steps.']);
    conv = 0; %overvrite the convergence flag to 0
    [s0_c, s0_nc] = test_convergence(s0_traj,epsilon);
    assignin('base', 'IRKA_s0_c', s0_c);
    assignin('base', 'IRKA_s0_nc', s0_nc);
    return
end

end

function s0=s0_vect(s0)
    % change two-row notation to vector notation
    if size(s0,1)==2
        temp=zeros(1,sum(s0(2,:)));
        for j=1:size(s0,2)
            k=sum(s0(2,1:(j-1))); k=(k+1):(k+s0(2,j));
            temp(k)=s0(1,j)*ones(1,s0(2,j));
        end
        s0=temp;
    end
    % sort expansion points
    s0 = cplxpair(s0);
    if size(s0,1)>size(s0,2)
        s0=transpose(s0);
    end
end

function es0 = stopping_crit(s0,s0_old)
    %     es0 = norm((s0-s0_old)./s0, 1)/sysr.n;
    es0 = norm((s0-s0_old), 1)/length(s0);
end

function [lh,th,ax] = initialize_plots(k,analysis_type,sys,s0)
    
    th = []; lh = []; ax = []; %initialization required(?)
    
   
    % plot initial shift configuration
    if strcmp(analysis_type,'complete');
        fh = figure;
            ax(1) = subplot(1,2,1);
                HFM_poles = eig(sys);
                plot(real(HFM_poles),imag(HFM_poles),'*k'); hold on
                plot([get(ax(1),'XLim'),0 0],[0 0 get(ax(1),'YLim')],'-k');
                xlabel('Re(lambda)'); ylabel('Im(lambda)');
                th(1) = title(['Poles of HFM and ROM (Iteration: ',num2str(k),')']);
            ax(2) = subplot(1,2,2);
                lh(1) = plot(real(s0),imag(s0),'*b'); hold on
                plot(-[get(ax(1),'XLim'),0 0],[0 0 get(ax(1),'YLim')],'-k');
                xlabel('Re(s0)'); ylabel('Im(s0)');
                th(2) = title(['IRKA shift selection (Iteration: ',num2str(k),')']);
    elseif strcmp(analysis_type,'red')
        fh = figure;
                lh(1) = plot(real(s0),imag(s0),'*b'); hold on
                xlabel('Re(s0)'); ylabel('Im(s0)');
                th = title(['IRKA shift selection (Iteration: ',num2str(k),')']);
    end  

end

function [h,lh] = plot_current(k,analysis_type,h,lh,th,ax,s0,s0_old,curr_poles,sys,sysr)

    if strcmp(analysis_type,'complete');
            if k == 1
                axes(ax(2));
                lh(2) = plot(real(s0_old),imag(s0_old),'or');

                axes(ax(1));
                lh(3) = plot(real(curr_poles),imag(curr_poles),'or');
                h = analyze_MOR(sys,sysr);
            else
                set(lh(3),'XData',real(curr_poles),'YData',imag(curr_poles));

                set(lh(2),'XData',real(s0_old), 'YData',imag(s0_old));

                analyze_MOR(sys,sysr,h);
            end
            set(lh(1),'XData',real(s0),'YData',imag(s0));
            axes(ax(2))
            legend([lh(2),lh(1)],'s_{0,k}','s_{0,k+1}','Location','EastOutside');

            set(th(1),'String',['Poles of HFM and ROM (Iteration: ',num2str(k),')']);
            set(th(2),'String',['IRKA shift selection (Iteration: ',num2str(k),')']);
        elseif strcmp(analysis_type,'red')
            if k == 1
                lh(2) = plot(real(s0_old),imag(s0_old),'or');
            else
                set(lh(2),'XData',real(s0_old), 'YData',imag(s0_old));
            end
            set(lh(1),'XData',real(s0),'YData',imag(s0));
            legend([lh(2),lh(1)],'s_{0,k}','s_{0,k+1}','Location','EastOutside');
            set(th,'String',['IRKA shift selection (Iteration: ',num2str(k),')']);
    end
end

function [s0_c, s0_nc] = test_convergence(s0,epsilon)
%
%   this function is used to determine which shifts have converged and
%   which have not
%
%   s0 is the (k+1 x n_rom) array of shifts over the iteration
s0_c = [];
s0_nc = [];

[kp,n] = size(s0);

for ii = 1:n
    if stopping_crit(s0(kp,ii),s0(kp-1,ii)) < 5e-2
        s0_c = [s0_c, s0(kp,ii)];
    else
        s0_nc = [s0_nc, s0(kp,ii)];
    end
end


end


function custompause(t)

%   pauses for user interaction for t = 'user' and pauses for some definite
%   time t if t is a double

if ischar(t) && strcmp(t,'user')
    pause()
else
    pause(t)
end
end