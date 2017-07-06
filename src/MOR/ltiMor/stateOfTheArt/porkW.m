function varargout = porkW(W,Sw,Lw,B)
% PORKW - Pseudo-Optimal Rational Krylov (Output)
%
% Syntax: 
%       [Ar,Br,Cr,Er] = porkW(W,Sw,Lw,B)
%       sysrPO = porkW(W,Sw,Lw,B)
%
% Description:
%       This function implements the pseudo-optimal rational Krylov
%       algorithm introduced by Wolf and Panzer [1,2].
%
%       Given a projection matrix W spanning an output Krylov subspace and
%       the respective matrices from the Sylvester equation
%
%       $A^T W - E^T W S_w^T - C^T L_w = 0$
%
%       This function computes the reduced order matrices corresponding to
%       the H2-pseudo-optimal reduced order model, i.e. a model
%       interpolating the original according to (W,Sw,Lw) and having
%       eigenvalues as mirror images of the shifts.
%
%       If only one output is specified, this function returns an ssRed
%       object. Otherwise, the reduced system matrices are returned.
%
% Input Arguments:
%       *Required Input  Arguments:* 
%       -W,Sw,Lw:        solution of  A.'*W - E.'*W*Sw.' - C.'*Lw = 0
%       -B:              input matrix of the original model
%
% Output Arguments: 
%       - sysrPO:         Pseudo-optimal reduced order model 
%       - Ar,Br,Cr,Er:    ROM matrices
%
% Examples:
%       Following code computes an H2-pseudo-optimal reduced order model
%       with an output Krylov subspace
%> sys = loadSss('building');
%> s0 = -eigs(sys,4,'sm').';
%> [sysr, ~, W] = rk(sys,[],s0);
%> [Lw, ~, Sw] = getSylvester(sys, sysr, W, 'W');
%> sysrPO = porkW(W,Sw,Lw,sys.B)
% 
% See Also: 
%       porkV, spark, rk, getSylvester
%
% References:
%       * *[1] Wolf (2014)*, H2 Pseudo-Optimal Moder Order Reduction
%       * *[2] Panzer (2014)*, Model Order Reduction by Krylov Subspace Methods
%              with Global Error Bounds and Automatic Choice of Parameters
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
% Authors:      Thomas Wolf, Heiko Panzer
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  13 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Computation
Pr_c = lyapchol(-Sw.', Lw.');
Cr = -Lw/Pr_c/Pr_c.';
Ar = Sw.'+Lw.'*Cr;
Br = W.'*B;
Er = eye(size(Ar));

%% Preparing output
if nargout == 1
    %%  Storing additional parameters
    %Stroring additional information about thr reduction in the object 
    %containing the reduced model:
    %   1. Define a new field for the Opts struct and write the information
    %      that should be stored to this field
    %   2. Adapt the method "parseParamsStruct" of the class "ssRed" in such a
    %      way that the new defined field passes the check
    
    Opts.originalOrder = size(W,1);
    varargout{1} = ssRed(Ar,Br,Cr,zeros(size(Cr,1),size(Br,2)),Er,'porkW',Opts);
else
    varargout{1} = Ar;
    varargout{2} = Br;
    varargout{3} = Cr;
    varargout{4} = Er;
end
