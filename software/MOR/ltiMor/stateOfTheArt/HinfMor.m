function [sysr, sysr0, Hinf] = HinfMor(sys, n, varargin) 
% HINFMOR - H-infinity reduction by tangential interpolation
% ------------------------------------------------------------------
%
% [sysr, Hinf] = HINFMOR(sys, varargin) 
% Inputs:       * sys: full oder model (sss)
%               * Opts (opt.) structure with execution parameters
% Outputs:      * sysr: reduced order model (sss)
%               * Hinf: H-infinity error (or indicator)
% ------------------------------------------------------------------
% USAGE:  
%
% See also IRKA
%
% ------------------------------------------------------------------
% REFERENCES:
% [1] Gugercin (2008), H2 model reduction for large-scale linear
%     dynamical systems
% [2] Flagg (2013), Interpolatory Hinf model reduction
% [3] Beattie (2014), Model reduction by rational interpolation
% ------------------------------------------------------------------
% This file is part of MORLab, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                   -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Last Change:  30 Oct 2015
% Copyright (c) ?
% ------------------------------------------------------------------

%%  Run IRKA
if sys.isSiso
    % initialize
    try s0 = -eigs(sys,n,'sm').'; catch , s0 = zeros(1,n); end
    % run IRKA
    [sysr0, ~, ~, ~, ~, Rt, Lt] = irka(sys,s0);
else %MIMO
    % initialize
    %   compute one step of tangential Krylov at 0 to get initial tangent 
    %   directions
    s0 = zeros(1,n); Rt = ones(sys.m,n); Lt = ones(sys.p,n);
    sysr = rk(sys,s0,s0,Rt,Lt);  [X,D,Y] = eig(sysr);
    Rt = full((Y.'*sysr.B).'); Lt = full(sysr.C*X); s0 = -diag(D).';
    %run IRKA
    [sysr0, ~, ~, ~, ~, Rt, Lt] = irka(sys,s0,Rt,Lt);
end

%%  Make Hinf correction
%
% steadyState: just take the current steady state error (can yield worse
%               results)
% steadyStateOpt: optimize over the difference in between the magnitude of
%                 the steady state error and Dr
% findGe0match: finds the feedthrough that best matches the steady-state
%               error amplitude
% normOpt:      Optimizes the actual inf norm of the error system!!
% steadyStateOpt+normOpt: combines the two

if nargin>2
    corrType = varargin{1};
else
    corrType = 'normOpt';
end

switch corrType
    case 'steadyState'
        warning(['This approach fails due to the same reasons stated for ',...
            'the forAll case']);
        %plus, Ge(0) changes depending on Dr
        
        Dr = freqresp(sys,0)-freqresp(sysr0,0);
        sysr = sss(sysr0.A+Lt.'*Dr*Rt, sysr0.B-Lt.'*Dr, ...
                    sysr0.C-Dr*Rt, Dr, sysr0.E);
    case 'steadyStateOpt'
        warning(['This approach fails due to the same reasons stated for ',...
            'the forAll case']);
        % however, is seems to work fine for build...
        
        G0 = freqresp(sys,0); %the only costly part
        Dr0 = G0-freqresp(sysr0,0);
        sysrfun = @(Dr) sss(sysr0.A+Lt.'*Dr*Rt, sysr0.B-Lt.'*Dr, ...
                                 sysr0.C-Dr*Rt, Dr, sysr0.E);
        cost = @(Dr) abs(...
                    abs(G0-...
                    freqresp(sysrfun(Dr),0)) - abs(Dr));
        DrOpt = fmincon(cost,Dr0)
        sysr = sss(sysr0.A+Lt.'*DrOpt*Rt, sysr0.B-Lt.'*DrOpt, ...
                   sysr0.C-DrOpt*Rt, Dr0, sysr0.E);
    case 'findGe0match'
        warning(['This approach fails since there does not seem to be a ',...
            'reduced order model that yields a completely flat magnitude',...
            'response in the error']);
        
        G0 = freqresp(sys,0); %the only costly part
        
        Dr0 = G0 - freqresp(sysr0,0) %get an initial feedthrough
        deltaDr = 100*abs(Dr0), nStep = 5000; DrSet = []; dSet = [];
        syse0 = sys-sysr0;
        Hinf0 = norm(syse0,Inf); d0 = abs(Dr0); dMin = d0; %initial error
        drawnow
        for Dr = Dr0-deltaDr:(2*deltaDr)/nStep:Dr0+deltaDr;
            sysr = sss(sysr0.A+Lt.'*Dr*Rt, sysr0.B-Lt.'*Dr, ...
                                 sysr0.C-Dr*Rt, Dr, sysr0.E);
            
            d = abs(abs(G0-freqresp(sysr,0)) - abs(Dr));
            if d < dMin
                dMin = d;
                dSet = [dSet, d];
                DrSet = [DrSet,Dr];
%                 syse = sys-sysr; sigma(syse,'Color',rand(1,3));
%                 drawnow, keyboard
            end
        end
        %   best feedthrough in term of minimizing the error between the
        %   response at 0 and Inf
        Dr = DrSet(end), dMin
        sysr = sss(sysr0.A+Lt.'*Dr*Rt, sysr0.B-Lt.'*Dr, ...
                                 sysr0.C-Dr*Rt, Dr, sysr0.E);
    case 'normOpt'
        
        sysrfun = @(Dr) sss(sysr0.A+Lt.'*Dr*Rt, sysr0.B-Lt.'*Dr, ...
                                 sysr0.C-Dr*Rt, Dr, sysr0.E);
        cost = @(Dr) norm(ss(sys-sysrfun(Dr)),Inf);
        warning('optimizing over the actual error norm');
        solver = 'fminsearch';
        switch solver
            case 'fminsearch'
            % initialization at steady-state error amplitude response
            % G0 = freqresp(sys,0); Dr0 = G0-freqresp(sysr0,0);
            Dr0 = 0;
            tic, [DrOpt, Hinf] = fminsearch(cost,Dr0), tOpt = toc
            
            case 'ga'
                options = gaoptimset('Display','iter','TimeLimit',5*60,...
                    'UseParallel',true, 'PopInitRange',[-1;1]);
                tic, [DrOpt, Hinf] = ga(cost,1,[],[],[],[],[],[],[],options); tOpt = toc
        end
        
        sysr = sss(sysr0.A+Lt.'*DrOpt*Rt, sysr0.B-Lt.'*DrOpt, ...
                   sysr0.C-DrOpt*Rt, DrOpt, sysr0.E);              
    case 'steadyStateOpt+normOpt'
        if sys.isSiso
            Lt = ones(1,sysr0.n); %new input/output vectors (from IRKA)
            Rt = ones(1,sysr0.n);
        else
            error('not for MIMO yet')
        end
        
        % initialization at steady-state error amplitude response
        G0 = freqresp(sys,0); Dr0 = G0-freqresp(sysr0,0);
        sysrfun = @(Dr) sss(sysr0.A+Lt.'*Dr*Rt, sysr0.B-Lt.'*Dr, ...
                                 sysr0.C-Dr*Rt, Dr, sysr0.E);
                             
        cost = @(Dr) abs(...
                    abs(G0-...
                    freqresp(sysrfun(Dr),0)) - abs(Dr));
        Dr0opt = fminsearch(cost,Dr0)
        
        cost = @(Dr) norm(sys-sysrfun(Dr),Inf);
        warning('optimizing over the actual error norm');
        tic, [DrOpt, Hinf] = fmincon(cost,Dr0opt), tOpt= toc
        sysr = sss(sysr0.A+Lt.'*DrOpt*Rt, sysr0.B-Lt.'*DrOpt, ...
                   sysr0.C-DrOpt*Rt, DrOpt, sysr0.E);
    case 'DRange'
       
        G0 = freqresp(sys,0); %the only costly part
        
        Dr0 = abs(G0 - freqresp(sysr0,0)); %get an initial feedthrough
        nStep = 50; DrSet = []; dSet = [];
        
        Dr = linspace(-Dr0/3,Dr0/3, nStep);
        normVec = zeros(1,length(Dr));
        minVal = inf; minDr = 0; %initializing
        normO = norm(ss(sys),inf); 
        for iDr = 1:length(Dr)
            sysr = sss(sysr0.A+Lt.'*Dr(iDr)*Rt, sysr0.B-Lt.'*Dr(iDr), ...
                   sysr0.C-Dr(iDr)*Rt, Dr(iDr), sysr0.E);
            normVec(iDr) = norm(ss(sys-sysr),inf)/normO;
            if normVec(iDr) < minVal
                minDr = Dr(iDr); minVal = normVec(iDr);
            end
        end
        figure; plot(Dr,normVec,'*-'); hold on; plot(minDr,minVal,'or');
        drawnow
    otherwise
        error('Specified Hinf optimization type not valid');
end
        

%% 
if nargout == 3
    Hinf = norm(sys-sysr,inf)/norm(sys,inf);
end

