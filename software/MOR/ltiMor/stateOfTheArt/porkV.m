function varargout = porkV(V,S_V,Crt,C)
% PORKV - Pseudo-Optimal Rational Krylov (Input)
%
% Syntax: 
%       [Ar,Br,Cr,Er] = porkV(V,S_V,Crt,C)
%
% Description:
%       This function implements the pseudo-optimal rational Krylov
%       algorithm introduced by Wolf and Panzer [1,2].
%
%       Given a projection matrix V spanning an input Krylov subspace and
%       the respective matrices from the Sylvester equation
%
%           $A V - E V S_V - B Crt = 0$
%
%       this function computes the reduced order matrices corresponding to
%       the H2-pseudo-optimal reduced order model, i.e. a model
%       interpolating the original according to (V,S_V,Crt) and having
%       eigenvalues as mirror images of the shifts.
%
%       If only one output is specified, this function returns an sss
%       object. Otherwise, the reduced system matrices are returned.
% 
% Input Arguments:
%       *Required Input Arguments:*
%       -V,S_V,Crt:      solution of  A*V - E*V*S_V - B*Crt = 0
%       -C:              output matrix of original model
%
% Output Arguments: 
%       - sysrPO:         Pseudo-optimal reduced order model 
%       - Ar,Br,Cr,Er:    ROM matrices
% 
% Examples:
%       Following code computes an H2-pseudo-optimal reduced order model
%       with an input Krylov subspace
%> sys = loadSss('building');
%> s0 = -eigs(sys,4,'sm').';
%> [sysr, V] = rk(sys,s0);
%> [Crt, ~, S] = getSylvester(sys, sysr, V);
%> sysrPO = porkV(V,S,Crt,sys.C)
% 
% See Also: 
%       porkW, spark, rk, getSylvester
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
% Last Change:  08 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


%%  Computation
Qr_c = lyapchol(-S_V', Crt');
Br = -Qr_c\(Qr_c'\Crt');
Ar = S_V+Br*Crt;
Cr = C*V;
Er = eye(size(Ar));

%% Preparing output
if nargout == 1
    varargout{1} = sss(Ar,Br,Cr,zeros(size(Cr,1),size(Br,2)),Er);
else
    varargout{1} = Ar;
    varargout{2} = Br;
    varargout{3} = Cr;
    varargout{4} = Er;
end
    
