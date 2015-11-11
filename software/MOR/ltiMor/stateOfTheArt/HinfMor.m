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
% Authors:      Alessandro Castagnotto
% Last Change:  10 Nov 2015
% Copyright (c) ?
% ------------------------------------------------------------------

%%  Run IRKA
if sys.isSiso
    % initialize
    try s0 = -eigs(sys,n,'sm').'; catch , s0 = zeros(1,n); end
    % run IRKA
    [sysr0, V, W, s0, ~, Rt, Lt] = irka(sys,s0);
else %MIMO
    % initialize
    %   compute one step of tangential Krylov at 0 to get initial tangent 
    %   directions
    s0 = -eigs(sys,n,'sm').'; Rt = ones(sys.m,n); Lt = ones(sys.p,n);
    sysr = rk(sys,s0,s0,Rt,Lt);  [X,D,Y] = eig(sysr);
    Rt = full((Y.'*sysr.B).'); Lt = full(sysr.C*X); s0 = -diag(D).';
    %run IRKA
    [sysr0, V, W, s0, ~, Rt, Lt] = irka(sys,s0,Rt,Lt);
end

%   Check moment matcing
rMatching = zeros(1,n); lMatching = zeros(1,n);
for iS = 1:length(s0)
    rMatching(iS) = norm((moments(sys,s0(iS),1)-moments(sysr0,s0(iS),1))*Rt(:,iS));
    lMatching(iS) = norm(Lt(:,iS).'*(moments(sys,s0(iS),1)-moments(sysr0,s0(iS),1)));
end
all([rMatching, lMatching]< 1e-5)

%   Transform (A- s0*E) to (s0*E- A)
sysr0.C = -sysr0.C; sysr0.B = -sysr0.B;
Rt = -Rt; Lt = -Lt;

%   Check moment matcing
rMatching = zeros(1,n); lMatching = zeros(1,n);
for iS = 1:length(s0)
    rMatching(iS) = norm((moments(sys,s0(iS),1)-moments(sysr0,s0(iS),1))*Rt(:,iS));
    lMatching(iS) = norm(Lt(:,iS).'*(moments(sys,s0(iS),1)-moments(sysr0,s0(iS),1)));
end
if ~all([rMatching, lMatching]< 1e-5)
    warning('Tangential moment matching failed?')
    keyboard
end

% Check that the generalized tangential directions are correct
R = getSylvester(sys,sysr0,-V); L = getSylvester(sys,sysr0,-W,'W'); 
if norm(Rt-R)> 1e-5 || norm(Lt-L)> 1e-5
    warning('Residuals could be wrong')
    keyboard
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
% steadyState+normOpt: Initializes optimization at Ge(0)
% DRange:       computes the actual error norm for a series of feedtrhoughs
%               and takes the one with minimum value

if nargin>2
    corrType = varargin{1};
else
    corrType = 'normOpt';
end

%   Parametrize reduced order models by Dr
sysrfun = @(Dr) sss(sysr0.A+Lt.'*Dr*Rt, sysr0.B+Lt.'*Dr, ...
                                 sysr0.C+Dr*Rt, Dr, sysr0.E);

switch corrType
    case 'steadyState'
        warning(['This approach fails due to the same reasons stated for ',...
            'the forAll case']);
        %plus, Ge(0) changes depending on Dr
        
        Dr = freqresp(sys,0)-freqresp(sysr0,0);
        sysr = sysrfun(Dr);
    case 'steadyStateOpt'
        warning(['This approach fails due to the same reasons stated for ',...
            'the forAll case']);
        % however, is seems to work fine for build...
        
        G0 = freqresp(sys,0); %the only costly part
        Dr0 = G0-freqresp(sysr0,0);
        
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
            sysr = sysrfun(Dr);
            
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
        sysr = sysrfun(Dr);
    case 'normOpt'
        cost = @(Dr) norm(sys-sysrfun(Dr),Inf);
        warning('optimizing over the actual error norm');
        solver = 'fminsearch';
        switch solver
            case 'fminsearch'
            %zero initialization
            Dr0 = zeros(sys.p,sys.m);
            warning off
            tic, [DrOpt, Hinf] = fminsearch(cost,Dr0), tOpt = toc
            warning on
            
            case 'ga'
                options = gaoptimset('Display','iter','TimeLimit',5*60,...
                    'UseParallel',true, 'PopInitRange',[-1;1]);
                tic, [DrOpt, Hinf] = ga(cost,1,[],[],[],[],[],[],[],options); tOpt = toc
        end
        
        sysr = sysrfun(DrOpt);  
        for iS = 1:length(s0)
            norm(moments(sys-sysr0,s0(iS),1)*Rt(:,iS)) %< 1e-5
            norm(Lt(:,iS).'*moments(sys-sysr0,s0(iS),1)) %< 1e-5
        end
        for iS = 1:length(s0)
            norm(moments(sys-sysr,s0(iS),1)*Rt(:,iS)) %< 1e-5
            norm(Lt(:,iS).'*moments(sys-sysr,s0(iS),1)) %< 1e-5
        end
    case 'steadyState+normOpt'
        % initialization at steady-state error amplitude response
        G0 = freqresp(sys,0); Dr0 = G0-freqresp(sysr0,0);
                             
        cost = @(Dr) norm(sys-sysrfun(Dr),Inf);
        solver = 'fminsearch';
        warning('optimizing over the actual error norm');
        warning off
        switch solver
            case 'fminsearch'
            %zero initialization
            tic, [DrOpt, Hinf] = fminsearch(cost,Dr0), tOpt = toc
            
            case 'ga'
                options = gaoptimset('Display','iter','TimeLimit',5*60,...
                    'UseParallel',true, 'PopInitRange',[-1;1]);
                tic, [DrOpt, Hinf] = ga(cost,1,[],[],[],[],[],[],[],options); tOpt = toc
        end
        warning on
        sysr = sysrfun(DrOpt);
        
        for iS = 1:length(s0)
            norm(moments(sys-sysr,s0(iS),1)*Rt(:,iS)) %< 1e-5
        end
    case 'DRange'
       
        G0 = freqresp(sys,0); 
        
        Dr0 = abs(G0 - freqresp(sysr0,0)); %get an initial feedthrough
        nStep = 100;
        
        Dr = linspace(-Dr0,Dr0, nStep); 
        normVec = zeros(1,length(Dr));
        normO = norm(ss(sys),inf); 
        %initializing
        normVec(1) = norm(ss(sys-sysr0),inf)/normO; 
        minVal = normVec(1); minDr = 0; Dr = [0,Dr];
        for iDr = 2:length(Dr)
            sysr = sysrfun(Dr(iDr));
            normVec(iDr) = norm(ss(sys-sysr),inf)/normO;
            if normVec(iDr) < minVal
                minDr = Dr(iDr); minVal = normVec(iDr);
            end
        end
        figure; plot(Dr,normVec,'*'); hold on; plot(minDr,minVal,'or');
        title(sprintf('minDr = %e',minDr));
        drawnow
    otherwise
        error('Specified Hinf optimization type not valid');
end
        
%% 
if nargout == 3
    Hinf = norm(sys-sysr,inf)/norm(sys,inf);
end

