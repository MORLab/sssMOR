function sysIMr = tbrDAEimproper(sys,nu,Bim,Cim)
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
%           -.cure.checkEVB:check if [EV,B_] has full column rank (or dual)
%                           [{true},false]
%           -.cure.sEVBTol: rank tolerance for [EV,B_] matrix (or dual)
%                           [{1e-16}/ positive float]
%           -.warn:         show warnings
%                           [{'0'} / '1']
%           -.w:            frequencies for analysis plots
%                           [{''} / '{wmin,wmax}' / vector of frequencies]
%           -.zeroThers:    value that can be used to replace 0 
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
    Def.svalTol = 1e-16;
   
    % create the options structure
    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end    
    
%%  Computation of cholesky factors of the improper Gramians

OmegaIMC = zeros(sys.n,sys.m*nu);
OmegaIMO = zeros(sys.n,sys.p*nu);


[L,U,a,o,S]=lu(sys.A,'vector');

OmegaIMC(o,1:sys.m) = U\(L\(S(:,a)\Bim));
CimT = Cim.';
OmegaIMO(:,1:sys.p) = (S(:,a)).'\(L.'\(U.'\(CimT(o,:))));

for idx = 1:nu-1
    OmegaIMC(o,idx*sys.m+1:(idx+1)*sys.m) = U\(L\(S(:,a)\(sys.E* OmegaIMC(o,(idx-1)*sys.m+1:idx*sys.m))));
    OmegaIMO(:,idx*sys.p+1:(idx+1)*sys.p) = (S(:,a)).'\(L.'\(U.'\(OmegaIMO(o,(idx-1)*sys.p+1:idx*sys.p))));
end

%% SVD
[Zl,Lambda,Zr] = svd(OmegaIMO.'*sys.A*OmegaIMC);

%%  Truncate zero sv
idx = find(diag(Lambda)<Opts.svalTol);

Zl(:,idx:end) = [];
Lambda(idx:end,idx:end) = [];
Zr(:,idx:end) = [];

LambdaInvSqrt = Lambda^(-0.5);
Wim = OmegaIMO*Zl*LambdaInvSqrt;
Vim = OmegaIMC*Zr*LambdaInvSqrt;

%%  Minimal realization
sysIMr = sss(eye(size(Vim,2)),Wim.'*sys.B,sys.C*Vim, sys.D,Wim.'*sys.E*Vim);
