function [sysr,sysSPr,sysIMr,varargout] = DAEmor(sys,varargin)
% Model Order Reduction of Structured DAEs
%
% Syntax:
%       sysr = DAEmor(sys,DAEtype)
%       sysr = DAEmor(sys,DAEtype,Opts)
%       sysr = DAEmor(sys,SpecProjFinC,SpecProjFinB,nu)
%       sysr = DAEmor(sys,SpecProjFinC,SpecProjFinB,nu,Opts)
%
% Description:
%       This function implements the reduction of structured DAEs by
%       1.) reduction of the strictly proper subsystem
%       2.) finding a minimal realization of the improper subsystem by BT
%
%       So far for step 1.) only CUREd SPARK and IRKA are implemented.
%       Other reduction methods may be added by expanding the function with
%       the respective Opts entry.
%
%       For some typical types of DAEs, the spectral projections can be
%       constructed by the function itself ('ind1se','ind2stokes'). In
%       all other cases, the left and right projectors onto the invariant
%       subspace of finite eigenvalues needs to be added.
%       
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:           a sss-object containing the LTI-DAE-system
%       *Optional Input Arguments:*
%       -DAEtype        a string containing either 'ind1se' or 'ind2stokes'
%       -SpecProjFinC:  function handle for projecting the output matrix C
%                       onto the deflating subspace of (lambda*E-A) corr.
%                       to the finite eigenvalues
%                       (for explicit knowledge of spectral projectors
%                       use: "C*SpecFinRight")
%       -SpecProjFinB:  function handle for projecting the input matrix B
%                       onto the deflating subspace of (lambda*E-A) corr.
%                       to the finite eigenvalues
%                       (for explicit knowledge of spectral projectors
%                       use: "SpecFinLeft*B")
%       -nu:            Index of the DAE-system
%       -Opts: A structure containing reduction and execution options
%               -.redfun:  reduction algorithm for the strictly proper part
%                          [{'cure'} / 'irka']
%               -.projSide: choose if C (def.) or B is projected to obtain a strictly proper realization
%                          [{'C'} / 'B'] 
%               -.(refer to *cure* or *irka* for other options)
%           
%
% Output Arguments:     
%       -sysr:      reduced (overall) system
%       -sysSPr:    reduced strictly proper subsystem
%       -sysIMr:    reduced improper subsystem
%       -varargout: additional outputs depending on the reduction method
%
% See Also: 
%       cure, spark, irka, projInd1se, projInd2stokes
%
% References:
%       * *[1] Stykel (2004)*, Gramian-Based Model Reduction for Descriptor Systems
%       * *[2] Gugercin, Stykel, Wyatt (2013)*, Model Reduction of Descriptor Systems by 
%               Interpolatory Projection Methods
%       * *[3] Seiwald, Castagnotto, Stykel, Lohmann (tbd)*, H2 Pseudo-Optimal 
%               Reduction of Structured DAEs by Rational Interpolation
%
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
% Authors:      Alessandro Castagnotto, Philipp Seiwald
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  03 Mar 2017
% Copyright (c) 2016,2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Parse input and load default parameters

% Opts parameters
% =================
Def.redfun      = 'cure';
Def.projSide    = 'C';

Def.irka.s0     = zeros(1,2);

if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
else
    Opts = struct();
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% Other inputs
% =============
if isempty(varargin)
    error('sssMOR:DAEmor:nargin','Not enough input arguments');
else
    if ischar(varargin{1})
        DAEtype = varargin{1};
    elseif length(varargin) == 3
        SpecProjFinC = varargin{1};
        SpecProjFinB = varargin{2};
        nu           = varargin{3};
    else
        error('sssMOR:DAEmor:usage','The number of input arguments is not compatible.')
    end
end

%% Load spectral projectors if not passed
if exist('DAEtype','var')
    switch DAEtype
        case 'ind1se'
            nu = 1;
            [SpecProjFinC, SpecProjFinB] = projInd1se(sys);
        case 'ind2stokes'
            nu = 2;
            [SpecProjFinC, SpecProjFinB] = projInd2stokes(sys);
        otherwise
            errror('sssMOR:DAE:wrongDAEtype','The DAEtype specified is either wrong or not supported')
    end
end
%% Obtain a realization for the strictly proper subsystem

sysSP = sys; %preallocate

switch Opts.projSide
    case 'B'
        %partitioning by projection of input matrix B
        sysSP.B = SpecProjFinB(sysSP.B);
        Bim = sys.B - sysSP.B;
        Cim = sys.C - SpecProjFinC(sys.C);
    case 'C'
        %partitioning by projection of input matrix B
        sysSP.C = SpecProjFinC(sysSP.C);
        Bim     = sys.B - SpecProjFinB(sys.B);
        Cim     = sys.C - sysSP.C;
end

%% Processing of strictly proper subsystem
switch Opts.redfun
    case 'cure'
        % run CUREd SPARK
        [sysSPr,sysSPrAll] = cure(sysSP,Opts);
    case 'irka'
        sysSPr = irka(sysSP,Opts.irka.s0);
end

%% Processing of improper subsystem
% find minimal realization of polynomial part by balanced truncation
sysIMr = tbrDAEimproper(sys,nu,Bim,Cim,Opts);

%% Obtain the reduced model from summing the subsystems
sysr            = sysSPr + sysIMr;
varargout{1}    = cellfun(@(sysSPr) sysSPr + sysIMr, sysSPrAll,'UniformOutput',false);


function sysIMr = tbrDAEimproper(sys,nu,Bim,Cim,Opts)
% Compute Minimal Realization of Improper Subsystem by Balanced Truncation
%
% Syntax:
%       sysIMr = tbrDAEimproper(sys,nu,Bim,Cim)
%       sysIMr = tbrDAEimproper(sys,nu,Bim,Cim,Opts)
%
% Description:
%       This function implements Lyapunov balanced truncation for the
%       improper subsystem of a DAE. The minimal realization is found by
%       truncation of all ZERO improper HSVs.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:   a sss-object containing the original LTI-DAE-system
%       -nu:    Index of the DAE-system
%       -Bim:   input matrix B projected onto the deflating subspace of 
%               (lambda*E-A) corr. to the finite eig. 
%       -Cim:   output matrix C projected onto the deflating subspace of 
%               (lambda*E-A) corr. to the finite eig. 
%       *Optional Input Arguments:*
%       -Opts: A structure containing following fields
%           -.svalTol:  tolerance for zero improper HSV (for truncation)
%                       [{1e-16}/ positive float]
%
% Output Arguments:     
%       -sysIMr: minimal realization of the improper subsystem
%
% Examples:
%
% See Also: 
%
% References:
%
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
    
%% Computation of cholesky factors of the improper Gramians
% preallocate memory
OmegaIMC = zeros(sys.n,sys.m*nu);
OmegaIMO = zeros(sys.n,sys.p*nu);

% alternative NOT WORKING code
% ----------------------------
% [L,U,a,o,S]=lu(sys.A,'vector');
% OmegaIMC(o,1:sys.m) = U\(L\(S(:,a)\Bim));
% CimT = Cim.';
% OmegaIMO(:,1:sys.p) = (S(:,a)).'\(L.'\(U.'\(CimT(o,:))));
% for idx = 1:nu+4
%     OmegaIMC(o,idx*sys.m+1:(idx+1)*sys.m) = U\(L\(S(:,a)\(sys.E* OmegaIMC(o,(idx-1)*sys.m+1:idx*sys.m))));
%     OmegaIMO(:,idx*sys.p+1:(idx+1)*sys.p) = (S(:,a)).'\(L.'\(U.'\(OmegaIMO(o,(idx-1)*sys.p+1:idx*sys.p))));
% end

% compute first block
OmegaIMC(:,1:sys.m) = sys.A\Bim;
OmegaIMO(:,1:sys.p) = (Cim/sys.A).';

% compute the following blocks
for idx = 1:nu-1
    OmegaIMC(:,idx*sys.m+1:(idx+1)*sys.m) = sys.A\(sys.E*OmegaIMC(:,(idx-1)*sys.m+1:idx*sys.m));
    OmegaIMO(:,idx*sys.p+1:(idx+1)*sys.p) = (((OmegaIMO(:,(idx-1)*sys.p+1:idx*sys.p)).'*sys.E)/sys.A).';
end

%% perform (thin) SVD
[Zl,Lambda,Zr] = svd(OmegaIMO.'*sys.A*OmegaIMC,'econ');

%% truncate zero improper HSV
% get improper HSVs (sorted after SVD)
imHSV = diag(Lambda);

% determine first zero improper HSV
firstzero_imHSV = 0;
for idx = 1:length(imHSV)
    if imHSV(idx) < Opts.svalTol    % compare with user defined tolerance
        firstzero_imHSV = idx;
        break
    end
end

% check, if there are any zero improper HSVs
if firstzero_imHSV > 0
    % check, if all improper hankel singular values are zero
    if(firstzero_imHSV == 1)
        % there is no improper subsystem -> return empty system object
        sysIMr = sss([],[],[],[],[]);
        return
    else
        % otherwise truncate 
        Zl = Zl(:,1:firstzero_imHSV-1);
        Lambda = Lambda(1:firstzero_imHSV-1,1:firstzero_imHSV-1);
        Zr = Zr(:,1:firstzero_imHSV-1);
    end
else
    % all improper HSVs are nonzero -> don't truncate anything
end

%% prepare projection
LambdaInvSqrt = Lambda^(-0.5);
Wim = OmegaIMO*Zl*LambdaInvSqrt;
Vim = OmegaIMC*Zr*LambdaInvSqrt;

%% return minimal realization of improper subsystem
sysIMr = sss(speye(size(Vim,2)),Wim.'*sys.B,sys.C*Vim,sys.D,Wim.'*sys.E*Vim);

