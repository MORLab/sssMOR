function sysr = CURE(sys,Opts)
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

%% Parse input and load default parameters

    % default values
    Def.warn = 0; %show warnings?
    Def.verbose = 0; %show progress text?
    Def.w = []; %frequencies for bode plot
    
    Def.zeroThres = 1e-4; % define the threshold to replace "0" by
        Def.CURE.red = 'SPARK'; %reduction algorithm
        Def.CURE.nk = 2; % reduced order at each step
        Def.CURE.stop = 'nmax'; %type of stopping criterion
        Def.CURE.stopval = round(sqrt(size(sys.a,1))); %default reduced order
        Def.CURE.init = 0; %shift initialization type
        Def.CURE.fact = 'V'; %error factorization
        Def.CURE.SE_DAE = 0; %SE-DAE reduction
        Def.CURE.test = 0; %execute analysis code?
        Def.CURE.gif = 0; %produce a .gif of the CURE iteration
        
    % create the options structure
    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end
        
    
         
%%  Plot for testing
if Opts.CURE.test
    fhOriginalSystem = nicefigure('CURE - Reduction of the original model');
    fhSystemBeingReduced = fhOriginalSystem; %the two coincide for the moment
    [mag,phase,w] = bode(sys,Opts.w); bodeOpts = {'Color',TUM_Blau,'LineWidth',2};
    redo_bodeplot(mag,phase,w,bodeOpts), hold on
    magHandle = subplot(2,1,1); magLim = ylim;
    phHandle  = subplot(2,1,2); phLim  = ylim;
    if Opts.CURE.gif, writeGif('create'); end
end

%%   Initialize some variables
[~,m] = size(sys.b);  p = size(sys.c,1);
Er_tot = []; Ar_tot = []; Br_tot = []; Cr_tot = []; B_ = sys.b; C_ = sys.c;
BrL_tot = zeros(0,p); CrL_tot = zeros(p,0); 
BrR_tot = zeros(0,m); CrR_tot = zeros(m,0);

sysr = sss(Ar_tot,Br_tot,Cr_tot,zeros(p,m),Er_tot);

%   Computation for semiexplicit index 1 DAEs (SE DAEs)
if Opts.CURE.SE_DAE
    [DrImp,A22InvB22] = implicitFeedthrough(sys,Opts.CURE.SE_DAE);
    
    % if we inted to used SPARK, the DAE has to be modified if DrImp~=0
    if DrImp && strcmp(Opts.CURE.red,'SPARK')
        B_ = adaptDaeForSpark(sys,Opts.CURE.SE_DAE,A22InvB22);
        
        % if Opts.CURE.test, add a new plot for the modified system
        if Opts.CURE.test
            %switch the plot to be updated in the loop
            fhSystemBeingReduced = nicefigure('CURE - modified DAE (strictly proper)');
            sys = sss(sys.a,B_,C_,0,sys.e);
            [mag,phase,w] = bode(sys,Opts.w); bodeOpts = {'Color',TUM_Blau,'LineWidth',2};
            redo_bodeplot(mag,phase,w,bodeOpts)  
            magHandle = subplot(2,1,1); magLim = ylim;
            phHandle  = subplot(2,1,2); phLim  = ylim;
            if Opts.CURE.gif, writeGif('create'); end
        end 
    end
else
    DrImp = zeros(size(sys.d));
end 

%   We reduce only the strictly proper part and add the feedthrough at the end   
Dr_tot = sys.d + DrImp;

%%   Start cumulative reduction
if Opts.verbose, fprintf('\nBeginning CURE iteration...\n'); end
while ~stopCrit(sys,sysr,Opts) && size(sysr.a,1)<=size(sys.a,1)
    
    %   Redefine the G_ system at each iteration
    sys = sss(sys.a,B_,C_,0,sys.e);
    
    %   Initializations
    s0 = initializeShifts(sys,Opts);
        
	% 1) Reduction
    switch Opts.CURE.fact
        case 'V'
            % V-based decomposition, if A*V - E*V*S - B*Crt = 0
            switch Opts.CURE.red
                case 'SPARK'               
                    [V,S_V,Crt] = SPARK(sys.a,sys.b,sys.c,sys.e,s0,Opts); 
                    
                    [Ar,Br,Cr,Er] = PORK_V(V,S_V,Crt,C_);                   
                case 'IRKA'
                    [sysr_temp,V,W] = IRKA(sys,s0);
                    
                    Ar = sysr_temp.a;
                    Br = sysr_temp.b;
                    Cr = sysr_temp.c;
                    Er = sysr_temp.e;
                    
                    Crt = [eye(m), zeros(m)]; 
                case 'RK+PORK'
                    [sysRedRK, V, ~, ~, Crt] = RK(sys,s0); 
                    S_V = sysRedRK.E\(sysRedRK.A-sysRedRK.B*Crt);
                    
                    [Ar,Br,Cr,Er] = PORK_V(V,S_V,Crt,C_);
                    
                    %   Adapt Cr for SE DAEs
                    Cr = Cr - DrImp*Crt;
                    %   Adapt Cr_tot for SE DAEs
                    if ~isempty(Cr_tot)
                        Cr_tot(:,end-n+1:end) = Cr_tot(:,end-n+1:end) + ...
                            DrImp*CrR_tot(:,end-n+1:end);
                    end
            end
            n = size(V,2);
        case 'W'
        % W-based decomposition, if A'*W - E'*W*SW' - C'*Brt' = 0
            switch Opts.CURE.red
                case 'SPARK'               
                    [W,S_W,Brt] = SPARK(sys.a',sys.c',sys.b',sys.e',s0,Opts);
                    Brt = Brt';
                    S_W = S_W';
                    
                    [Ar,Br,Cr,Er] = PORK_W(W,S_W,Brt,B_);
                case 'RK+PORK'
                    [sysRedRK, ~, W, ~, ~, ~, Brt] = RK(sys,[],s0);
                    S_W = sysRedRK.E'\(sysRedRK.A'-Brt*sysRedRK.C);
                    
                    [Ar,Br,Cr,Er] = PORK_W(W,S_W,Brt,B_);  
                    
                    %   Adapt Br for SE-DAEs
                    Br = Br - Brt*DrImp;
                    %   Adapt Cr_tot for SE DAEs
                    if ~isempty(Br_tot)
                        Br_tot(end-n+1:end,:) = Br_tot(end-n+1:end,:) + ...
                            BrR_tot(end-n+1:end,:)*DrImp;
                    end
                    
            end
            n = size(W,2);
    end

    
	%Er = W.'*E*V;  Ar = W.'*A*V;  Br = W.'*B_;  Cr = C_*V;
    Er_tot = blkdiag(Er_tot, Er);
    Ar_tot = [Ar_tot, BrL_tot*Cr; Br*CrR_tot, Ar]; %#ok<*AGROW>
    Br_tot = [Br_tot; Br]; Cr_tot = [Cr_tot, Cr];
    
    if Opts.CURE.fact=='V'
        B_ = B_ - sys.e*(V*(Er\Br));    % B_bot
        BrL_tot = [BrL_tot; zeros(n,p)];    BrR_tot = [BrR_tot; Br];
        CrL_tot = [CrL_tot, zeros(p,n)];    CrR_tot = [CrR_tot, Crt];
    elseif Opts.CURE.fact=='W'
        C_ = C_ - Cr/Er*W.'*sys.e;		% C_bot
        BrL_tot = [BrL_tot; Brt];   BrR_tot = [BrR_tot; zeros(n,m)];
        CrL_tot = [CrL_tot, Cr];    CrR_tot = [CrR_tot, zeros(m,n)];
    end

    sysr    = sss(Ar_tot, Br_tot, Cr_tot, zeros(p,m), Er_tot);
    
    if Opts.CURE.test
        sysr_bode = sysr; %sysr_bode = sss(sysr);
        figure(fhSystemBeingReduced);
        bode(sysr_bode,w,'--','Color',TUM_Orange);
        subplot(2,1,1)
        title(sprintf('n_{red} = %i',size(sysr.a,1)));
        set(magHandle,'YLim', magLim); set(phHandle,'YLim', phLim);
        if Opts.CURE.gif, writeGif('append'); end
    end
end

%   Add the feedthrough term before returning the reduced system
sysr.D = Dr_tot;

%%  Finishing execution
if Opts.verbose,fprintf('Stopping criterion satisfied. Exiting CURE...\n\n');end
if Opts.CURE.test
        sysr_bode = sysr;
        figure(fhOriginalSystem);
        bode(sysr_bode,w,'-','Color',TUM_Gruen,'LineWidth',2);
        subplot(2,1,1)
        title(sprintf('n_{red} = %i',size(sysr.a,1)));
        
        if Opts.CURE.gif, writeGif('append'), end
end

%--------------------------AUXILIARY FUNCTIONS---------------------------
function stop = stopCrit(sys,sysr,opts)
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
function s0 = initializeShifts(sys,opts)
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
function [DrImp,A22InvB22] = implicitFeedthrough(sys,dynamicOrder)
%   compute the implicit feedthrough of a SE DAE

    B22 = sys.b(dynamicOrder+1:end);
    C22 = sys.c(dynamicOrder+1:end);
    if norm(B22)>0 && norm(C22)>0
        %this is not suffiecient for Dr~=0, but it is necessary
        
        %   Computing A22\B22 since it can be reused before SPARK
        A22InvB22 = sys.a(dynamicOrder+1:end,dynamicOrder+1:end)\B22;
        if norm(A22InvB22)>0
            DrImp = -C22*A22InvB22;
        else
            DrImp = zeros(p,m);
        end
    else
        A22InvB22 = zeros(size(B22));
        DrImp = zeros(p,m);
    end
function Bnew = adaptDaeForSpark(sys,dynamicOrder,A22InvB22)
% adapt B if the system is SE-DAE with Dr,imp ~=0

    [~,A12] = partition(sys.a,dynamicOrder);
    B11 = sys.b(1:dynamicOrder);
    
    Bnew = [ B11 - A12*A22InvB22;
          zeros(size(sys.a,1)-dynamicOrder,size(B11,2))];
      
%         Alternatively, the same can be done with C
%         C11 = sys.c(1:opts.CURE.SE_DAE);
%         C22 = sys.c(opts.CURE.SE_DAE+1:end);
%         sysCURE.C = [ C11- C22*(A22\A21),zeros(size(C22))];

function writeGif(gifMode)
    filename = 'CURE.gif';
    dt = 1.5;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    switch gifMode
        case 'create'
            imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',dt);
        case 'append'
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',dt);
        otherwise
            error('Invalid gifMode')
    end