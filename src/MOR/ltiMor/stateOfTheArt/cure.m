function sysr = cure(sys,Opts)
% CURE - CUmulative REduction framework
%
% Syntax:
%       sysr = CURE(sys)
%       sysr = CURE(sys,Opts)
%
% Description:
%       This function implements the CUmulative REduction framework
%       (CURE) introduced by Panzer and Wolf (see [1,2]).
%
%       Using the duality between Sylvester equation and Krylov subspaces, the 
%       error is factorized at each step of CURE and only the high-dimensional
%       factor is reduced in a subsequent step of the iteration.
%
%       Currently, this function supports following reduction strategies at
%       each step of CURE:
%       spark (def.), irka, rk+pork (pseudo-optimal reduction)
%
%       You can reduce index 1 DAEs in semiexplicit form by selecting the
%       appropriate option. See [3] for more details.
%
%       //Note: Currently CUREd SPARK works only for SISO systems.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: An sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -Opts: A structure containing following fields
%           -.cure.redfun:  reduction algorithm
%                           [{'spark'} / 'irka' / 'rk+pork']
%           -.cure.nk:      reduced order at each iteration 
%                           [{'2'} / positive integer]
%           -.cure.fact:    factorization mode 
%                           [{'V'} / 'W']
%           -.cure.init:    shift initialization mode 
%                           [{'sm'} / 'zero' / 'lm' / 'slm']
%           -.cure.stop:    stopping criterion
%                           [{'nmax'} / 'h2Error']
%           -.cure.stopval: value according to which the stopping criterion is evaluated
%                           [{'round(sqrt(sys.n))'} / positive integer]
%           -.cure.verbose: display text during cure 
%                           [{'0'} / '1']
%           -.cure.SE_DAE:  reduction of index 1 semiexplicit DAE 
%                           [{'0'} / '1']
%           -.cure.test:    execute analysis code 
%                           [{'0'} / '1']
%           -.cure.gif:     produce a .gif file of the CURE iteration
%                           [{'0'} / '1']
%           -.cure.maxIter: maximum number of CURE iterations
%                           [{'20'} / positive integer]
%           -.warn:         show warnings
%                           [{'0'} / '1']
%           -.w:            frequencies for analysis plots
%                           [{''} / '{wmin,wmax}' / vector of frequencies]
%           -.zeroThres:    value that can be used to replace 0 
%                           [{'1e-4'} / postivie float]
%
% Output Arguments:     
%       -sysr: Reduced system
%
% Examples:
%       By default, cure reduces a given model sys to a reduced order of
%       sqrt(sys.n) by steps of nk = 2 using mespark (model function based
%       spark)
%> sys = loadSss('building');
%> sysr = cure(sys);
%
%       The behavior of the function can be highly customized using the
%       option structure Opts
%
%> Opts.cure = struct('nk',4, 'redfun', 'irka', 'verbose', 1, 'stopval',12);
%> sysr = cure(sys,Opts)
% 
% See Also: 
%       spark, rk, irka, porkV, porkW, getSylvester
%
% References:
%       * *[1] Panzer (2014)*, Model Order Reduction by Krylov Subspace Methods
%              with Global Error Bounds and Automatic Choice of Parameters
%       * *[2] Wolf (2014)*, H2 Pseudo-Optimal Moder Order Reduction
%       * *[3] Castagnotto et al. (2015)*, Stability-preserving, adaptive
%              model order reduction of DAEs by Krylov subspace methods
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  13 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Parse input and load default parameters
    % default values
    Def.warn = 0;%show warnings?
    Def.cure.verbose = 0; %show progress text?
    Def.w = []; %frequencies for bode plot
    Def.zeroThres = 1e-4; % define the threshold to replace "0" by
        Def.cure.redfun = 'spark'; %reduction algorithm
        Def.cure.nk = 2; % reduced order at each step
        Def.cure.stop = 'nmax'; %type of stopping criterion
        Def.cure.stopval = round(sqrt(sys.n)); %default reduced order
            if ~isEven(Def.cure.stopval), Def.cure.stopval = Def.cure.stopval +1;end
        Def.cure.init = 'sm'; %shift initialization type
        Def.cure.fact = 'V'; %error factorization
        Def.cure.SE_DAE = 0; %SE-DAE reduction
        Def.cure.test = 0; %execute analysis code?
        Def.cure.gif = 0; %produce a .gif of the CURE iteration
        Def.cure.maxIter = 20; %maximum number of iterations
        
    % create the options structure
    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end              
    
    % make sure reduced order does not exceed original
    if Opts.cure.maxIter > sys.n/Opts.cure.nk
        Opts.cure.maxIter = floor(sys.n/Opts.cure.nk);
    end
%%  Plot for testing
if Opts.cure.test
    fhOriginalSystem = figure('Name','CURE - Reduction of the original model');
    fhSystemBeingReduced = fhOriginalSystem; %the two coincide for the moment
    [m,w] = freqresp(sys,Opts.w);
    bode(frd(m,w),'-b'), hold on
    axH = findall(gcf,'type','axes');
    magHandle = axH(3); magLim = get(magHandle,'YLim');
    phHandle = axH(2); phLim  = get(phHandle,'YLim');
    if Opts.cure.gif, writeGif('create'); end
end
%%   Initialize some variables
[~,m] = size(sys.b);  p = size(sys.c,1);
Er_tot = []; Ar_tot = []; Br_tot = []; Cr_tot = []; B_ = sys.b; C_ = sys.c;
BrL_tot = zeros(0,p); CrL_tot = zeros(p,0); 
BrR_tot = zeros(0,m); CrR_tot = zeros(m,0);

sysr = sss(Ar_tot,Br_tot,Cr_tot,zeros(p,m),Er_tot);
%%   Computation for semiexplicit index 1 DAEs (SE DAEs)
if Opts.cure.SE_DAE
    [DrImp,A22InvB22] = implicitFeedthrough(sys,Opts.cure.SE_DAE);
    
    % if we inted to used spark, the DAE has to be modified if DrImp~=0
    if DrImp && strcmp(Opts.cure.redfun,'spark')
        B_ = adaptDaeForSpark(sys,Opts.cure.SE_DAE,A22InvB22);
        
        % if Opts.cure.test, add a new plot for the modified system
        if Opts.cure.test
            %switch the plot to be updated in the loop
            fhSystemBeingReduced = figure('Name', 'CURE - modified DAE (strictly proper)');
            sys = sss(sys.a,B_,C_,0,sys.e);
            [m,w] = freqresp(sys,Opts.w);
            bode(frd(m,w),bodeOpts)  
            axH = findall(gcf,'type','axes');
            magHandle = axH(3); magLim = get(magHandle,'YLim');
            phHandle = axH(2); phLim  = get(phHandle,'YLim');
            if Opts.cure.gif, writeGif('create'); end
        end 
    elseif DrImp
        error(['Cumulative reduction of non strictly proper index 1 DAEs ' ...
                'in semiexplicit form is currently supported only for spark'])
    end
else
    DrImp = zeros(size(sys.d));
end 

%   We reduce only the strictly proper part and add the feedthrough at the end   
Dr_tot = sys.d + DrImp;

%%   Start cumulative reduction
if Opts.cure.verbose, fprintf('\nBeginning CURE iteration...\n'); end

iCure = 0; %iteration counter
while ~stopCrit(sys,sysr,Opts) && iCure < Opts.cure.maxIter
    iCure = iCure + 1;
    if Opts.cure.verbose, fprintf('\tCURE iteration %03i\n',iCure');end
    %   Redefine the G_ system at each iteration
    sys = sss(sys.a,B_,C_,0,sys.e);
    
    %   Initializations
    [s0,Opts] = initializeShifts(sys,Opts,iCure);        
	% 1) Reduction
    switch Opts.cure.fact
        case 'V'
            % V-based decomposition, if A*V - E*V*S - B*Rv = 0
            switch Opts.cure.redfun
                case 'spark'               
                    [V,Sv,Rv,~,usedOpts] = spark(sys,s0,Opts); 
                    
                    sysrTemp = porkV(V,Sv,Rv,C_); 
                    [Ar,Br,Cr,~,Er] = dssdata(sysrTemp);
                case 'irka'
                    [sysrTemp,V,W,~,~,~,~,~,Rv] = irka(sys,s0);
                                      
                    [Ar,Br,Cr,~,Er] = dssdata(sysrTemp);
                    
                case 'rk+pork'
                    [~, V, ~, ~, Sv, Rv] = rk(sys,s0);
                    sysrTemp = porkV(V,Sv,Rv,C_);
                    [Ar,Br,Cr,~,Er] = dssdata(sysrTemp);
                    
                    %   Adapt Cr for SE DAEs
                    Cr = Cr - DrImp*Rv;
                    %   Adapt Cr_tot for SE DAEs
                    if ~isempty(Cr_tot)
                        Cr_tot(:,end-n+1:end) = Cr_tot(:,end-n+1:end) + ...
                            DrImp*CrR_tot(:,end-n+1:end);
                    end
                otherwise 
                    error('The selected reduction scheme (Opts.cure.redfun) is not availabe in cure');
            end
            n = size(V,2);
        case 'W'
        % W-based decomposition, if A.'*W - E.'*W*Sw.' - C.'*Lw = 0
            switch Opts.cure.redfun
                case 'spark'               
                    [W,Sw,Lw,~,usedOpts] = spark(sys.',s0,Opts);
                    Sw = Sw.'; %make Sw from Sw^T
                    
                    sysrTemp = porkW(W,Sw,Lw,B_);
                    [Ar,Br,Cr,~,Er] = dssdata(sysrTemp);
                case 'irka'
                    [sysrTemp,V,W,~,~,~,~,~,~,~,~,Lw] = irka(sys,s0);
                    
                    [Ar,Br,Cr,~,Er] = dssdata(sysrTemp);
                    
                case 'rk+pork'
                    [~, ~, W, ~, ~, ~, ~, Sw, Lw] = rk(sys,[],s0);
                    
                    sysrTemp = porkW(W,Sw,Lw,B_); 
                    [Ar,Br,Cr,~,Er] = dssdata(sysrTemp);
                    
                    %   Adapt Br for SE-DAEs
                    Br = Br - Lw.'*DrImp;
                    %   Adapt Cr_tot for SE DAEs
                    if ~isempty(Br_tot)
                        Br_tot(end-n+1:end,:) = Br_tot(end-n+1:end,:) + ...
                            BrR_tot(end-n+1:end,:)*DrImp;
                    end  
            end
            n = size(W,2);
    end
    %%  Cumulate the matrices and define sysr
	%Er = W.'*E*V;  Ar = W.'*A*V;  Br = W.'*B_;  Cr = C_*V;
    Er_tot = blkdiag(Er_tot, Er);
    Ar_tot = [Ar_tot, BrL_tot*Cr; Br*CrR_tot, Ar]; %#ok<*AGROW>
    Br_tot = [Br_tot; Br]; Cr_tot = [Cr_tot, Cr];
    
    if Opts.cure.fact=='V'
        B_ = B_ - sys.e*(V*(Er\Br));    % B_bot
        BrL_tot = [BrL_tot; zeros(n,p)];    BrR_tot = [BrR_tot; Br];
        CrL_tot = [CrL_tot, zeros(p,n)];    CrR_tot = [CrR_tot, Cr];
    elseif Opts.cure.fact=='W'
        C_ = C_ - Cr/Er*W.'*sys.e;		% C_bot
        BrL_tot = [BrL_tot; Lw.'];   BrR_tot = [BrR_tot; zeros(n,m)];
        CrL_tot = [CrL_tot, Cr];    CrR_tot = [CrR_tot, zeros(m,n)];
    end
    
    %%  Storing additional parameters
    %Stroring additional information about the reduction in the object 
    %containing the reduced model:
    %   1. Define a new field for the Opts struct and write the information
    %      that should be stored to this field
    %   2. Adapt the method "checkParamsStruct" of the class "ssRed" in such a
    %      way that the new defined field passes the check
      
    if strcmp(Opts.cure.redfun,'irka') || strcmp(Opts.cure.redfun,'rk+pork')
        usedOpts = sysrTemp.reductionParameters{end,1}.params;
        usedOpts.cure = Opts.cure;
    end
    usedOpts.originalOrder = sysr.n;
    if isa(sysr,'ssRed')
        sysr = ssRed(strcat('cure_',Opts.cure.redfun),usedOpts,Ar_tot, ...
                     Br_tot, Cr_tot, zeros(p,m), Er_tot,sysr.reductionParameters);
    else
        sysr = ssRed(strcat('cure_',Opts.cure.redfun),usedOpts,Ar_tot, ...
                     Br_tot, Cr_tot, zeros(p,m), Er_tot);
    end
    
    % display
    if Opts.cure.test
        sysr_bode = sysr; 
        figure(fhSystemBeingReduced);
        bode(sysr_bode,w,'--r');
        magHandle;
        title(sprintf('n_{red} = %i',size(sysr.a,1)));
        set(magHandle,'YLim', magLim); set(phHandle,'YLim', phLim);
        if Opts.cure.gif, writeGif('append'); end
    end
end

%%   Was maxIter achieved?
if iCure >= Opts.cure.maxIter
    warning('Iterations count reached maxIter. You might want to increase maxIter or the convergence criterion')
end
%%   Add the feedthrough term before returning the reduced system
sysr.D = Dr_tot;

%%  Finishing execution
if Opts.cure.verbose,fprintf('Stopping criterion satisfied. Exiting cure...\n\n');end
if Opts.cure.test
        sysr_bode = sysr;
        figure(fhOriginalSystem);
        bode(sysr_bode,w,'-g');
        subplot(2,1,1)
        title(sprintf('n_{red} = %i',size(sysr.a,1)));
        
        if Opts.cure.gif, writeGif('append'), end
end

%% --------------------------AUXILIARY FUNCTIONS---------------------------
function stop = stopCrit(sys,sysr,opts)
%   computes the stopping criterion for CURE iteration
switch opts.cure.stop
    case 'h2Error'
        if sys.isBig
            warning('System size might be to large for stopping criterion');
        end
        if isempty(sys.h2Norm), sys.h2Norm = norm(sys,2); end
        if sysr.n>0 %avoid computing when initializing
            stop = (norm(sys-sysr,2)/sys.h2Norm <= opts.cure.stopval);
        else
            stop = 0;
        end
    case 'nmax'
        stop = (sysr.n >=opts.cure.stopval);
    otherwise
        error('The stopping criterion chosen does not exist or is not yet implemented');
end
function [s0,Opts] = initializeShifts(sys,Opts,iCure)
 %%   parse
 if Opts.cure.init ==0, Opts.cure.init = 'zero'; end
 
 %%   initialization of the shifts
 if ~ischar(Opts.cure.init) %initial shifts were defined
     if length(Opts.cure.init) == Opts.cure.nk %correct amount for iteration
         s0 = Opts.cure.init;
     elseif length(Opts.cure.init)> Opts.cure.nk
         firstIndex = (iCure-1)*Opts.cure.nk+1;
            if firstIndex + Opts.cure.nk -1 > length(Opts.cure.init)
                % overlap with the previous set of shifst
                Opts.cure.init = [Opts.cure.init, Opts.cure.init];
            end
            % new set of shifts
            s0 = Opts.cure.init(firstIndex:firstIndex+Opts.cure.nk-1);
     else
         error('The initial vector of shifts s0 passed to CURE is invalid');
     end
     
 else %choose initialization option     
     %  choose the number of shifts to compute
     if strcmp(Opts.cure.stop,'nmax')
         ns0 = Opts.cure.stopval; %just as many as needed
     else %another stopping criterion was chosen
         ns0 = round(sqrt(sys.n));
         if ~isEven(ns0), ns0 = ns0+1; end
     end     
     %  compute the shifts
     try
         switch Opts.cure.init
             case 'zero' %zero initialization
                 Opts.cure.init = Opts.zeroThres*ones(1,ns0);
             case 'sm' %smallest magnitude eigenvalues
                 Opts.cure.init = -eigs(sys.a,sys.e,ns0,0, ...
                        struct('tol',1e-6,'v0',sum(sys.b,2)));
             case 'lm' %largest magnitude eigenvalues
                 Opts.cure.init = -eigs(sys.a,sys.e,ns0,'lm', ...
                        struct('tol',1e-6,'v0',sum(sys.b,2)));
             case 'slm' %(def.) smallest and largest eigs
                 %  decide how many 'lm' and 'sm' eigs to compute
                 if ns0 <=4
                     nSm = 2; nLm = 2;
                 else %ns0 > 4
                    if isEven(ns0)
                        nSm = ns0/2; nLm = nSm;
                        s0 = [];
                    else
                        nSm = (ns0-1)/2; nLm = nSm;
                        s0 = 0; %initialize the first shift at 0
                    end
                    if ~isEven(nSm), nSm = nSm +1; nLm = nLm -1;end
                 end
                 Opts.cure.init = [ -eigs(sys.a,sys.e,nSm,'sm', ...
                                      struct('tol',1e-6,'v0',sum(sys.b,2)));...
                                    -eigs(sys.a,sys.e,nLm,'lm', ...
                                      struct('tol',1e-6,'v0',sum(sys.b,2)))]';
         end
     catch err
         warning([getReport(err,'basic'),' Using 0.'])
         Opts.cure.init = Opts.zeroThres*ones(1,ns0);
     end
     s0 = Opts.cure.init(1:Opts.cure.nk);
 end 
    %   make sure the initial values for the shifts are complex conjugated
    if mod(nnz(imag(s0)),2)~=0 %if there are complex valued shifts...
        % find the s0 which has no compl.conj. partner
        % (assumes that only one shifts has no partner)
        s0(abs(imag(s0)-imag(sum(s0)))<1e-16) = 0;
    end
    %     cplxpair(s0); %only checking
    %   replace 0 with thresh, where threshold is a small number
    %   (sometimes the optimizer complaints about cost function @0)
    s0(s0==0)=Opts.zeroThres;  
function [DrImp,A22InvB22] = implicitFeedthrough(sys,dynamicOrder)
%   compute the implicit feedthrough of a SE DAE

    B22 = sys.b(dynamicOrder+1:end);
    C22 = sys.c(dynamicOrder+1:end);
    if norm(B22)>0 && norm(C22)>0
        %this is not suffiecient for Dr~=0, but it is necessary
        
        %   Computing A22\B22 since it can be reused before spark
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
% analogously, the same effect could be achieved by modifying C

    [~,A12] = partition(sys.a,dynamicOrder);
    B11 = sys.b(1:dynamicOrder);
    
    Bnew = [ B11 - A12*A22InvB22;
          zeros(size(sys.a,1)-dynamicOrder,size(B11,2))]; 
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
function isEven = isEven(a)
    isEven = round(a/2)== a/2;
function y = TUM_Blau()
    y = [0 101 189]/255;
function y = TUM_Gruen()
    y = [162 173 0]/255;
function y = TUM_Orange()
    y = [227 114 34]/255;