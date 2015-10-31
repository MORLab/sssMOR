function [sysr, Hinf] = HinfMor(sys, n, varargin) 
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
    sysr0 = irka(sys,s0);
else %MIMO
    % initialize
    %   compute one step of tangential Krylov at 0 to get initial tangent 
    %   directions
    s0 = zeros(1,n); Rt = ones(sys.m,n); Lt = ones(sys.p,n);
    sysr = rk(sys,s0,s0,Rt,Lt);  [X,D,Y] = eig(sysr);
    Rt = full((Y.'*sysr.B).'); Lt = full(sysr.C*X); s0 = -diag(D).';
    %run IRKA
    sysr0 = irka(sys,s0,Rt,Lt);
end

%%  Make Hinf correction
%
% steadyState: just take the current steady state error (can yield worse
%               results)
% steadyStateOpt: optimize over the difference in between the magnitude of
%                 the steady state error and Dr
corrType = 'steadyStateOpt';

switch corrType
    case 'steadyState'
        Lt = ones(1,sysr0.n);
        Rt = ones(1,sysr0.n);
        Dr = freqresp(sys,0)-freqresp(sysr0,0);
        sysr = sss(sysr0.A+Lt.'*Dr*Rt, sysr0.B-Lt.'*Dr, ...
                    sysr0.C-Dr*Rt, Dr, sysr0.E);
    case 'steadyStateOpt'
        Lt = ones(1,sysr0.n);
        Rt = ones(1,sysr0.n);
        G0 = freqresp(sys,0); %the only costly part
        Dr0 = G0-freqresp(sysr0,0);
        cost = @(Dr) abs(...
                    abs(G0-...
                    freqresp(sss(sysr0.A+Lt.'*Dr*Rt, sysr0.B-Lt.'*Dr, ...
                                 sysr0.C-Dr*Rt, Dr, sysr0.E),0)) - Dr);
        [DrOpt,optCost] = fmincon(cost,Dr0)
        sysr = sss(sysr0.A+Lt.'*DrOpt*Rt, sysr0.B-Lt.'*DrOpt, ...
                   sysr0.C-DrOpt*Rt, Dr0, sysr0.E);
end
        

%%  developing code
%analyze results
% figure; bode(sys); hold on; bode(sysr,'r--')
isstable(sysr0)
isstable(sysr)

syse0 = sys-sysr0;
Hinf0 = norm(syse0,Inf)/norm(sys,Inf)

syse = sys-sysr;
Hinf = norm(syse,Inf)/norm(sys,Inf)

% figure; bode(sys-sysr0,'-.b'); hold on; bode(sys-sysr,'-r');
errorPlot(sys,{sysr0,sysr});

