function [s0_inp,s0_out,Rt,Lt] = getShifts(sys,sysr,nShifts,s0_inp,Rt,s0_out,Lt,basis1,basis2,Opts)
% GETSHIFTS - get shift update for iterative rtional Krylov subspace methods within the iteration process
%
% Syntax for 'eigs'-method:
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS([],sysr,nShifts,s0_inp,[],[],[],[],[])
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS([],sysr,nShifts,s0_inp,Rt,[],[],[],[])
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS([],sysr,nShifts,s0_inp,[],s0_out,[],[])
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS([],sysr,nShifts,[],[],s0_out,Lt,[],[])
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS([],sysr,nShifts,s0_inp,Rt,s0_out,Lt,[],[])
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS([],sysr,nShifts,s0_inp,...,Opts)
%
% Syntax for 'adaptive'-method:
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS(sys,sysr,[],s0_inp,[],[],[],V,[])
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS(sys,sysr,[],s0_inp,Rt,[],[],V,[])
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS(sys,sysr,[],s0_inp,[],s0_out,V,W)
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS(sys,sysr,[],[],[],s0_out,Lt,[],W)
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS(sys,sysr,[],s0_inp,Rt,s0_out,Lt,V,W)
%       [s0_inp,s0_out,Rt,Lt]       = GETSHIFTS(sys,sysr,[],s0_inp,...,Opts)
%       
%
% Description:
%
%       Generates a 1 x nShift row vector or nShift x 1 column vector of shift frequencies s0_inp 
%       according to the chosen strategy, which can be defined in Opts.strategy.
%       Besides, if necessary, tangential directions Rt and Lt as well as shift frequencies s0_out for an output
%       Krylov subspace are available. 
%       The GETSHIFTS function is mainly important for the Cumulative Rational Krylov Subspace
%       Method (CRKSM) of sssMOR to obtain e. g. the solution of a Lyapunov
%       equation or a reduced Model.
%
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:                   An sss-object containing the original LTI system (necessary for 'adaptive' method)
%       -sysr:                  sss-object containing an already reduced LTI system (necessary for both methods)
%       -nShifts                Number of new shifts (necessary for 'eigs' method)
%       -s0_inp:                previously used shift frequencies (necessary for both methods)
%
%		*Optional Input Arguments:*
%       -s0_out:                expansion points for output Krylov subspace
%       -Rt/Lt:                 right/left tangential directions (MIMO case)
%       -V/W:                   orthogonal Krylov input/output subspace bases
%       -Opts:                  a structure containing following options
%           -.strategy:         strategy for calculating new shifts [{'adaptive'} / 'eigs']
%                -'adaptive':       use the adaptive  shift generation method after Druskin (one new shift) 
%                -'eigs':           use a defined number (nShifts) of Ritz values of the current reduced order model 
%           -.multDir:          use the multiple tangential directions method after Druskin [{0} / 1]
%
%
% Output Arguments:
%       -s0_inp:                new, enlarged input shift vector 
%       -s0_out:                new, enlarged output shift vector, if it is
%                               not specified, GETSHIFTS delivers an empty array
%       -Rt:                    new, enlarged right tangential direction matrix, if it is
%                               not specified, GETSHIFTS delivers an empty array
%       -Lt:                    new, enlarged left tangential direction matrix, if it is
%                               not specified, GETSHIFTS delivers an empty array
%
%
% See Also: 
%       initializeShifts, crksm, mess_para, mess_get_ritz_vals 
%
% References:
%       * *[1] Druskin, Simoncini (2011)*, Adaptive Rational Krylov Subspaces
%       for large-scale dynamical systems
%       * *[2] Druskin, Simoncini, Zaslavsky (2014)*, Adaptive Tangential
%       Interpolation in Rational Krylov Subspaces for MIMO Dynamical Systems
%       * *[3] Kürschner (2016)*, Efficient Low-Rank Solution of Large-Scale Matrix Equations
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Paul Heidenreich, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Nov 2017
% Copyright (c) 2016-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Create Def-struct containing default values/options
Def.getShiftsStrategy    = 'adaptive';       % strategy for shift generation: [{'adaptive'} / 'eigs' / 'projection']
Def.multDir              =  0;               % choose strategy for tangential dierections [{0} / 1] 

% additional function intern defaults, only interessting for usage within CRKSM-function
Def.B_                   = [];               % recycle Bbot if it is already available
Def.Er_inv_Ar            = [];               % recycle Er^-1*Ar if it is already available
Def.C_                   = [];               % recycle Cbot if it is already available
Def.Er_invT_ArT          = [];               % recycle Er^-T*ArT if it is already available

%% Parsing of Inputs/Preprocessing
% create the options structure and check tangential directions
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end
clear Def

%% Compute New Shift/Shifts
% get new shifts: 'eigs' / 'adaptive'
switch Opts.getShiftsStrategy
    case 'eigs'   
        % getShifts-function uses initializeShifts-function, overwrite strategies
        Opts.initShiftsStrategy = Opts.getShiftsStrategy;
        % call INITIALIZESHIFTS to perform eigs
        if isempty(Rt) && isempty(s0_out) 
            snewOut = [];   Rtnew = [];     Ltnew = [];
            snewInp = initializeShifts(sysr,nShifts,1,Opts);
            snewInp = snewInp'; % row-vector
        elseif isempty(s0_out)
            snewOut = [];   Ltnew = [];
            [snewInp,Rtnew] = initializeShifts(sysr,nShifts,1,Opts);
            snewInp = snewInp';
        elseif isempty(s0_inp)
            snewInp = [];   Rtnew = [];    
            [~,~,snewOut,Ltnew] = initializeShifts(sysr,nShifts,1,Opts);
            snewOut = snewOut';
        elseif isempty(Rt) && ~isempty(s0_out) 
            Rtnew = [];     Ltnew = [];
            [snewInp,~,snewOut,~] = initializeShifts(sysr,nShifts,1,Opts);
            snewInp = snewInp';
            snewOut = snewOut';
        else
            [snewInp,Rtnew,snewOut,Ltnew] = initializeShifts(sysr,nShifts,1,Opts);
            snewInp = snewInp';
            snewOut = snewOut';
        end
                      
    case 'adaptive'
        % shift generation after adaptive Druskin method
        if isempty(s0_out)
            snewOut = [];        Ltnew = [];
            [snewInp,Rtnew] = newParaInp(sys,sysr,basis1,s0_inp,Rt,Opts);
            if isreal(snewInp) && ~isreal(Rtnew)
                snewInp = [snewInp snewInp];
            elseif ~isreal(snewInp) && isreal(Rtnew)
                Rtnew = [Rtnew Rtnew];
            end
        elseif isempty(s0_inp)
            snewInp = [];        Rtnew = [];
            [snewOut,Ltnew] = newParaOut(sys,sysr,basis1,s0_out,[],Lt,Opts);
            if isreal(snewOut) && ~isreal(Rtnew)
                snewOut = [snewOut snewOut];
            elseif ~isreal(snewOut) && isreal(Ltnew)
                Ltnew = [Ltnew Ltnew];
            end
        else
            [snewInp,Rtnew] = newParaInp(sys,sysr,basis1,s0_inp,Rt,Opts);
            [snewOut,Ltnew] = newParaOut(sys,sysr,basis2,s0_out,[],Lt,Opts);
            % check cocnsistence of vector/matrix length
            if size(snewInp,2) ~= size(snewOut,2)
                if size(snewInp,2) > size(snewOut,2)
                    snewOut = [snewOut snewOut];
                else
                    snewInp = [snewInp snewInp];
                end
            end
            if size(Rtnew,2) ~= size(Ltnew,2)
                if size(Rtnew,2) > size(Ltnew,2)
                    Ltnew = [Ltnew Ltnew];
                else
                    Rtnew = [Rtnew Rtnew];
                end
            end
        end    
end % end switch
% convert to single for numerical accuracy
snewInp = single(snewInp);      Rtnew = single(Rtnew);
snewOut = single(snewOut);      Ltnew = single(Ltnew); 

% check complex values
soldInp = snewInp;
soldOut = snewOut;
snewInp = cplxpair(snewInp);
snewOut = cplxpair(snewOut);   
[~,cplxSortingInp] = ismember(snewInp,soldInp);       [~,cplxSortingOut] = ismember(snewOut,soldOut); 
if ~isempty(Rtnew),     Rtnew = Rtnew(:,cplxSortingInp);        end                       
if ~isempty(Ltnew),     Ltnew = Ltnew(:,cplxSortingOut);        end

% check if s0_inp is column or row vector or matrix, build enlarged output parametres
if isrow(s0_inp) || isrow(s0_out)
    if ~isempty(s0_inp),    s0_inp = double([s0_inp snewInp]);   end
    if ~isempty(s0_out),    s0_out = double([s0_out snewOut]);   end
else 
    if ~isempty(s0_inp),    s0_inp = double([s0_inp reshape(snewInp,[1,size(snewInp,2)])]);                end    
    if ~isempty(s0_out),    s0_out = double([s0_out reshape(snewOut,[1,size(snewOut,1)])]);                end
end
Rt = double([Rt Rtnew]);        Lt = double([Lt Ltnew]);  
end % end of getShifts

%% ***************************** AUXILIARY ********************************

function [snewInp,Rtnew] = newParaInp(sys,sysr,V,s0_inp,Rt,Opts)  
    % compute Ritz-Values of reduced system 
    ritzVal = eig(sysr); 
    ritzVal = single(ritzVal);
    
    % delete already used Ritz values in s0_inp 
    [~,idx] = ismember(single(ritzVal),-single(unique(s0_inp)));
    idx = find(idx);
    ritzVal(idx,:) = [];
    ritzVal = sort(ritzVal);
    
    %Spiegeln,falls instabile eigs
    if all(real(ritzVal)) < 0
        idxUnstable = real(ritzVal)<0; 
        ritzVal(idxUnstable) = -ritzVal(idxUnstable);
    end
    
    % check if there are NaN or Inf entries in ritzVal
    if sum(isnan(ritzVal)) ~= 0 || sum(isinf(ritzVal)) ~= 0
        warning('The reduced system seems to be unstable because there are Nan and/or Inf Ritz Values. An error may occur! Try a onesided projection only with V or W basis to avoid instability');
    end
    ritzVal = double(ritzVal);
    
    % build convex hull
    if ~isreal(ritzVal)
        specSet = sort([s0_inp'; -ritzVal]);
        chull = convhull(real(specSet),imag(specSet));     % bulid convex hull of spectral set
        specSet = specSet(single(chull)); 
    else
        specSet = -ritzVal;
    end

    % qr-decompositons, residual matrices
    if isempty(Opts.B_)
        Er_inv_Br = solveLse(sysr.E,sysr.B);
        Opts.B_ = sys.B-(sys.E*V)*Er_inv_Br;
        Opts.Er_inv_Ar = solveLse(sysr.E,sysr.A,Opts);   
    end
    res0 = norm(Opts.B_);

    % solve max-problem for new shift, initialize and set variables
    resNorm = zeros(size(specSet,1),1);
    resNorm_last = 0;   
    for ii = 1:1:size(specSet,1)
        Y = solveLse((sysr.A-specSet(ii,1)*sysr.E),sysr.B);
        resNorm = norm(((sys.A*V-sys.E*V*Opts.Er_inv_Ar)*Y-Opts.B_))/res0; 
        if isempty(resNorm_last) || resNorm > resNorm_last
            resNorm_last = resNorm;
            snewInp = specSet(ii,1);
            if ~isreal(snewInp)
                snewInp = [snewInp conj(snewInp)];   
            end
            Ypic = Y;
        end
    end
    
    % compute (multiple) tangential directions
    if ~isempty(Rt)
        res = (sys.A*V-sys.E*V*Opts.Er_inv_Ar)*Ypic-Opts.B_;
        [~,S,rSingVec] = svd(res); 
        for ii = 1:1:size(S,2)
            if Opts.multDir == false
                resNorm(ii,1) = norm(res*rSingVec(:,ii));
                [~,index] = max(resNorm);
                Rtnew = rSingVec(:,index);
                if ~isreal(Rtnew)
                    Rtnew = [Rtnew conj(Rtnew)];
                end
            else              
                if S(ii,ii) > 0.1*S(1,1)
                    rt_new = rSingVec(:,ii);
                    if ~isreal(rt_new)
                        rt_new = [rt_new conj(rt_new)];
                    end
                    Rtnew = [Rt rt_new];
                    snewInp = [s0_inp S(ii,ii)];
                end    
            end
        end % end of for-loop
    else
        Rtnew = [];
    end
end
 
function [snewOut,Ltnew] = newParaOut(sys,sysr,W,s0_out,snewOut,Lt,Opts)
    % compute Ritz-Values of reduced system 
    ritzVal = eig(sysr);
    ritzVal = single(ritzVal);

    % delete already used Ritz values in s0_inp 
    [~,idx] = ismember(single(ritzVal),-single(unique(s0_out)));
    idx = find(idx);
    ritzVal(idx,:) = [];
    ritzVal = sort(ritzVal);

    %Spiegeln,falls instabile eigs
    if all(real(ritzVal)) < 0
        idxUnstable = real(ritzVal)<0; 
        ritzVal(idxUnstable) = -ritzVal(idxUnstable);
    end
    
    % check if there are NaN or Inf entries in ritzVal
    if sum(isnan(ritzVal)) ~= 0 || sum(isinf(ritzVal)) ~= 0
        warning('The reduced system seems to be unstable because there are Nan and/or Inf Ritz Values. An error may occur! Try a onesided projection only with V or W basis to avoid instability');
    end
    ritzVal = double(ritzVal);

    % build convex hull
    if ~isreal(ritzVal)
        specSet = sort([s0_out'; -ritzVal]);
        chull = convhull(real(specSet),imag(specSet));     % bulid convex hull of spectral set
        specSet = specSet(single(chull)); 
    else
        specSet = -ritzVal;
    end

    % qr-decompositons, residual matrices
    if isempty(Opts.C_)
        Er_invT_CrT = solveLse(sysr.E',sysr.C');
        Opts.C_ = sys.C'-(sys.E'*W)*Er_invT_CrT;
        Opts.Er_invT_ArT = solveLse(sysr.E',sysr.A',Opts);     
    end
    res0 = norm(Opts.C_);

    % solve max-problem for new shift, initialize and set variables
    resNorm = zeros(size(specSet,1),1);
    resNorm_last = 0;   
    for ii = 1:1:size(specSet,1)
        Y = solveLse((sysr.A-specSet(ii,1)*sysr.E)',sysr.C');
        resNorm = norm(((sys.A'*W-sys.E'*W*Opts.Er_invT_ArT)*Y-Opts.C_))/res0;         
        if isempty(resNorm_last) || resNorm > resNorm_last
            resNorm_last = resNorm;
            snewOut = specSet(ii,1);
            if ~isreal(snewOut)
                snewOut = [snewOut conj(snewOut)];  
            end
            Ypic = Y;
        end
    end
    
    % compute (multiple) tangential directions
    if ~isempty(Lt)
        res = (sys.A'*W-sys.E'*W*Opts.Er_invT_ArT)*Ypic-Opts.C_;
        [~,S,rSingVec] = svd(res); 
        for ii = 1:1:size(S,2)
            if Opts.multDir == false
                resNorm(ii,1) = norm(res*rSingVec(:,ii));
                [~,index] = max(resNorm);
                Ltnew = rSingVec(:,index);
                if ~isreal(Ltnew)
                    Ltnew = [Ltnew conj(Ltnew)];
                end
            else              
                if S(ii,ii) > 0.1*S(1,1)
                    lt_new = rSingVec(:,ii);
                    if ~isreal(lt_new)
                        lt_new = [lt_new conj(lt_new)];
                    end
                    Ltnew = [Rt lt_new];
                    snewOut = [s0_inp S(ii,ii)];
                end    
            end
        end % end of for-loop
    else
        Ltnew = [];
    end
end
