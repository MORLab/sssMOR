function sysr = CURE(sys,opts)
%CURE - CUmulative REduction framework
% ------------------------------------------------------------------
% sysr = CURE(sys, opts)
% Inputs:       * sys: an sss-object containing the LTI system
%               * opts: a structure containing following options
%                   * CURE.red:  reduction algorithm - {'RK','IRKA','SPARK',...}
%                   * CURE.nk:   reduced order at each iteration
%                   * CURE.fact: factorization mode - {'V','W'}
%                   * CURE.stop: stopping criterion - {'H2-error','HInf-error',...
%                               'nmax'}
%                   * CURE.stopval: value according to which the stopping
%                               criterion is evaluated (e.g. 'error'->tol, etc.)
%                   * test:     produce plots and verbatim to analyze
% Outputs:      * sysr: reduced system
% ------------------------------------------------------------------
% USAGE:  This function implements the CUmulative REduction framework
% (CURE) introduced by Panzer and Wolf (see [1,2]).
%
% Using the duality between Sylvester equation and Krylov subspaces, the 
% error is factorized at each step of CURE and only the high-dimensional
% factor is reduced in a subsequent step of the iteration.
%
% See also SPARK, RK, IRKA, PORK_V, PORK_W.
%
% ------------------------------------------------------------------
% REFERENCES:
% [1] Panzer (2014): Model Order Reduction by Krylov Subspace Methods
%     with Global Error Bounds and Automatic Choice of Parameters
% [2] Wolf (2014): H2 Pseudo-Optimal Moder Order Reduction
% ------------------------------------------------------------------
% This file is part of MORLab, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                      -> MORLab@tum.de <-
% ------------------------------------------------------------------
% Authors:      Heiko K.F. Panzer, Alessandro Castagnotto, 
%               Maria Cruz Varona
% Last Change:  26 April 2015
% Copyright 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

% Parse input and load default parameters
    % default values
    def.test = 0; %execute analysis code?
    def.warn = 0; %show warnings?
    def.verbose = 0; %show progress text?
    def.w = []; %frequencies for bode plot
    
    def.zeroThres = 1e-4; % define the threshold to replace "0" by
        def.CURE.red = 'SPARK'; %reduction algorithm
        def.CURE.nk = 2; % reduced order at each step
        def.CURE.stop = 'nmax'; %type of stopping criterion
        def.CURE.init = 0; %shift initialization type
        def.CURE.fact = 'V'; %error factorization
        def.CURE.SE_DAE = 0; %SE-DAE reduction
                def.SPARK.type = 'model'; %SPARK type
                def.MESPARK.ritz = 1;
                def.MESPARK.pertIter = 5; % # iteration at which perturbation begins
                def.MESPARK.maxIter = 20; %maximum number of model function updates
                
    opts = parseOpts(opts,def);
         

if opts.test
    fh = nicefigure('Reduction with CURE');
    if ~isfield(opts,'w'), opts.w = []; end
    [mag,phase,w] = bode(sys,opts.w); bodeOpts = {'Color',TUM_Blau,'LineWidth',2};
    redo_bodeplot(mag,phase,w,bodeOpts), hold on    
end

%   Initialize some variables
[~,m] = size(sys.b);  p = size(sys.c,1);
Er_tot = []; Ar_tot = []; Br_tot = []; Cr_tot = []; B_ = sys.b; C_ = sys.c;
BrL_tot = zeros(0,p); CrL_tot = zeros(p,0); 
BrR_tot = zeros(0,m); CrR_tot = zeros(m,0);

sysr = sss(Ar_tot,Br_tot,Cr_tot,zeros(p,m),Er_tot);

%   We reduce only the strictly proper part and add the feedthrough at the end   
Dr_tot = sys.d;

%   Start cumulative reduction
if opts.verbose, fprintf('\nBeginning CURE iteration...\n'); end
while ~stopcrit(sys,sysr,opts) && size(sysr.a,1)<=size(sys.a,1)
    
    %   Redefine the G_ system at each iteration
    sys = sss(sys.a,B_,C_,0,sys.e);
    
    %   Initializations
    s0 = initialize_shifts(sys,opts);
        
	% 1) Reduction
    switch opts.CURE.fact
        case 'V'
            % V-based decomposition, if A*V - E*V*S - B*Crt = 0
            switch opts.CURE.red
                case 'SPARK'               
                    [V,S_V,Crt] = SPARK(sys.a,sys.b,sys.c,sys.e,s0,opts); 
                    
                    [Ar,Br,Cr,Er] = PORK_V(V,S_V,Crt,C_);
                    
                case 'IRKA'
                    [sysr_temp,V,W] = IRKA(sys,s0);
                    
                    Ar = sysr_temp.a;
                    Br = sysr_temp.b;
                    Cr = sysr_temp.c;
                    Er = sysr_temp.e;
                    
                    Crt = [eye(m), zeros(m)]; 
            end
            n = size(V,2);
        case 'W'
        % W-based decomposition, if A'*W - E'*W*SW' - C'*Brt' = 0
            switch opts.CURE.red
                case 'SPARK'               
                    [W,S_W,Brt] = SPARK(sys.a',sys.c',sys.b',sys.e',s0,opts);
                    Brt = Brt';
                    S_W = S_W';
                    
                    [Ar,Br,Cr,Er] = PORK_W(W,S_W,Brt,B_);
            end
            n = size(W,2);
    end

    
	%Er = W.'*E*V;  Ar = W.'*A*V;  Br = W.'*B_;  Cr = C_*V;
    Er_tot = blkdiag(Er_tot, Er);
    Ar_tot = [Ar_tot, BrL_tot*Cr; Br*CrR_tot, Ar]; %#ok<*AGROW>
    Br_tot = [Br_tot; Br];  Cr_tot = [Cr_tot, Cr];
    if opts.CURE.fact=='V'
        B_ = B_ - sys.e*(V*(Er\Br));    % B_bot
        BrL_tot = [BrL_tot; zeros(n,p)];    BrR_tot = [BrR_tot; Br];
        CrL_tot = [CrL_tot, zeros(p,n)];    CrR_tot = [CrR_tot, Crt];
    elseif opts.CURE.fact=='W'
        C_ = C_ - Cr/Er*W.'*sys.e;		% C_bot
        BrL_tot = [BrL_tot; Brt];   BrR_tot = [BrR_tot; zeros(n,m)];
        CrL_tot = [CrL_tot, Cr];    CrR_tot = [CrR_tot, zeros(m,n)];
    end

    sysr    = sss(Ar_tot, Br_tot, Cr_tot, zeros(p,m), Er_tot);
    
    if opts.test
        sysr_bode = sysr; %sysr_bode = sss(sysr);
        figure(fh);
        bode(sysr_bode,w,'--','Color',TUM_Orange);
        subplot(2,1,1)
        title(sprintf('n_{red} = %i',size(sysr.a,1)));
    end
end

if opts.verbose,fprintf('Stopping criterion satisfied. Exiting CURE...\n\n');end
if opts.test
        sysr_bode = sysr;
        figure(fh);
        bode(sysr_bode,w,'-','Color',TUM_Gruen,'LineWidth',2);
        subplot(2,1,1)
        title(sprintf('n_{red} = %i',size(sysr.a,1)));
end

% sysbot	= dss(A, B_, C_, zeros(p,m), E);	% caution: large-scale!
% 
% % truncate non controllable/observable states in sysrL, sysrR
% i = find(any(Ar_tot(any(BrL_tot~=0,2),:),1));
% sysrL   = dss(Ar_tot(i,i), BrL_tot(i,:), CrL_tot(:,i), eye(p), Er_tot(i,i));
% i = find(any(Ar_tot(:,any(CrR_tot~=0,1)),2));
% sysrR   = dss(Ar_tot(i,i), BrR_tot(i,:), CrR_tot(:,i), eye(m), Er_tot(i,i));

%--------------------------AUXILIARY FUNCTIONS---------------------------
function stop = stopcrit(sys,sysr,opts)
%   computes the stopping criterion for CURE iteration
switch opts.CURE.stop
    case 'H2-error'
        if size(sys,1)>500
            warning('System size might be to large for stopping criterion');
        end
        stop = (norm(sys-sysr,2)<= opts.CURE.stopval);
    case 'nmax'
        stop = (size(sysr.a,1)>=opts.CURE.stopval);
    otherwise
        error('The stopping criterion chosen does not exist or is not yet implemented');
end
function s0 = initialize_shifts(sys,opts)
 if ~isscalar(opts.CURE.init)
     s0 = opts.CURE.init;
 else
     switch opts.CURE.init
         case 0
             %zero
             s0 = opts.zeroThres*ones(1,opts.CURE.nk);
         case 1
             %random
             s0 = (.5+rand)*[1,1] + 100i*[1, -1]*rand;
         case 2
             % Convex combination of largest and smallest imaginary part eigs
             try
                s1 = -eigs(sys.a,sys.e,opts.CURE.nk,'li')';
             catch 
                warning('Largest imaginary eigs not found! Using rand*1e2 +-rand*1e2i instead.')
                s1 = rand*1e2*ones(1,2)+1e2*1i*rand*[1,-1];
             end
             try
                s2 = -eigs(sys.a,sys.e,opts.CURE.nk,'si')';
             catch 
                 warning('Smallest imaginary eigs not found! Using 0 instead.')
                 s2 = 0;
             end
    %          % test convex combination
    %          fh = figure;
    %          plot(-real(eig(sys)),imag(eig(sys)),'o');hold on;
    %          for ii = 1:100
             if all(real(s1)>= real(s2))
                s0 = (real(s1)+rand*(real(s2)-real(s1))) +...
                    1i*(imag(s2)+.5*rand*(imag(s1)-imag(s2)));
             else
                 s0 = (real(s2)+rand*(real(s1)-real(s2))) +...
                    1i*(imag(s2)+.5*rand*(imag(s1)-imag(s2)));
             end
    %          plot(real(s0),imag(s0),'xr');
    %          end
    %          pause;
    %          close(fh);

             s0 = s0';
         case 3
             % On a line between smallest imaginary and largest imaginary 
             s1 = -eigs(sys.a,sys.e,opts.CURE.nk,'li');
             try
                s2 = -eigs(sys.a,sys.e,opts.CURE.nk,'si');
             catch err
                 warning('Smallest imaginary eigs not found! Using 0 instead.')
                 s2 = zeros(2,1);
             end
             % test convex combination
             fh = figure;
             plot(-real(eig(sys)),imag(eig(sys)),'o');hold on;
             for ii = 1:100
                 s0 = s2 + (s1-s2)*rand;
                 plot(real(s0),imag(s0),'xr');
             end
             pause;
             close(fh);

             s0 = s0';
         case 4
             % As case 3, but projected onto real axis
             s1 = -eigs(sys.a,sys.e,opts.CURE.nk,'li');
             try
                s2 = -eigs(sys.a,sys.e,opts.CURE.nk,'si');
             catch err
                 warning('Smallest imaginary eigs not found! Using 0 instead.')
                 s2 = zeros(2,1);
             end
    %          % test convex combination
    %          fh = figure;
    %          plot(-real(eig(sys)),imag(eig(sys)),'o');hold on;
    %          for ii = 1:100
                 s0 = real(s2 + (s1-s2)*rand);
    %              plot(real(s0),imag(s0),'xr');
    %          end
    %          pause;
    %          close(fh);

             s0 = s0';
         case 5
             % random real shifts
             s0 = [rand,5*rand];

         case 6 %used for transmission line
    %          p = [6e7,1e20];
    %          s_opt = p(1)+[1 -1]*sqrt(p(1)^2-p(2))

             s1 = 6e4*ones(1,2)+5e7*1i*[1,-1];
             s2 = 1e4*ones(1,2)+5e5*1i*[1,-1];

             s0 = ((real(s1)+rand*(real(s2)-real(s1))) +...
                    1i*(imag(s2)+rand*(imag(s1)-imag(s2))));
     end
    %   replace 0 with thresh, where threshold is a small number
    %   (sometimes the optimizer complaints about cost function @0)
    s0(s0==0)=opts.zeroThres;
 end
                    