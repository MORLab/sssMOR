function sysr = projectiveMor(sys,V,varargin)
% PROJECTIVEMOR - Reduces dynamic model by projection
% 
% Syntax:
%       sysr				= PROJECTIVEMOR(sys,V)
%		sysr				= PROJECTIVEMOR(sys,V,W)
%		sysr				= PROJECTIVEMOR(sys,...,Opts)
% 
% Description:
%       This function generates a reduced order model sysr by applying
%       (Petrov)-Galerkin projection to the model in sys using user-defined
%       projection matrices V (and W). Sys must be an sss-object.
%
%       If only one projection matrix V is specified, then the reduced
%       order model is obtained by Galerkin projection (W=V), otherwise by
%       Petrov-Galerkin projection.
% 
%       By default, the reduced order model is obtained by Petrov-Galerkin
%       projection of the type
%
%       $$ A_r = W^TAV,\; E_r = W^TEV,\; B_r = W^TB,\; C_r = CV,\; D_r = D $$
%
%       However, if the optional parameter Opts.trans is set to 'H', then
%       $$ W^T $$ is replaced by $$ W^H $$, i.e. the conjugate transpose of W.
%
% Input Arguments:
%		*Required Input Arguments:*
%		-sys: 		n-dimensional model of dynamic system (sss)
%		-V: 		projection Matrix (n-by-r, where r is the reduced
%                   order), enough to perfom Galerkin projection
%
%		*Optional Input Arguments:*
%		-W:   		additional projectio Matrix (same size as V) to perform
%                   Petrov-Galerkin projection
%		-Opts:      A structure containing following fields
%			-.trans:  	choose how W should be transposed;
% 						[ {T} / H  ]
%                       //Note: Make sure the choice is consistent to the rest of your code.
%
% Output Arguments:
%       -sysr:      Reduced order model 
%
% Examples:
%		This code loads a benchmark model included in the toolbox and
%		performs Galerkin projection
%
%> sys = loadSss('building');
%> s0 = 1; V = (sys.A-s0*sys.E)\sys.B;
%> sysr = projectiveMor(sys,V)
%
%       To perform Petrov-Galerkin projection using *conjugate transposition*
%       use
%
%> s0 = 1i; 
%> V = (sys.A-s0*sys.E)\sys.B; 
%> W = (sys.A-s0*sys.E)'\sys.C'; 
%> sysr = projectiveMor(sys,V,W,struct('trans','T'))
%
%
% See Also: 
%		rk, moments, tbr
%
% References:
%		* *[1] Grimme (1997)*, Krylov projection methods for model reduction
%		* *[2] Antoulas (2005)*, Approximation of large-scale dynamical systems
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  06 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Parsing

% Opts structure
Def.trans = 'T'; 

if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
    
    Opts = parseOpts(Opts,Def);
else
    Opts = Def;
end       

% Galerkin projection?
if isempty(varargin) || isempty(varargin{1}) %no or empty W passed
    W = V;
else
    W = varargin{1};
end

%%  Storing additional parameters
%Stroring additional information about thr reduction in the object 
%containing the reduced model:
%   1. Define a new field for the Opts struct and write the information
%      that should be stored to this field
%   2. Adapt the methods "checkParamsStruct" and "parseParamsStruct" of 
%      class "ssRed" in such a way that the new defined field passes the 
%      check
Opts.originalOrder = sys.n;

%%  Projection
switch Opts.trans
    case 'T'
        if isa(sys,'ssRed')
            sysr = ssRed('projectiveMor',Opts,W.'*sys.A*V, W.'*sys.B, sys.C*V, sys.D, W.'*sys.E*V,sys.reductionParameters);
        else
            sysr = ssRed('projectiveMor',Opts,W.'*sys.A*V, W.'*sys.B, sys.C*V, sys.D, W.'*sys.E*V);
        end
    case 'H'
        if isa(sys,'ssRed')
            sysr = ssRed('projectiveMor',Opts,W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V,sys.reductionParameters);
        else
            sysr = ssRed('projectiveMor',Opts,W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);
        end
end