function n = l2norm(sys)
% L2NORM - L2-norm of a dynamical system (ssRed)
% 
% Syntax:
%		n                   = L2NORM(sys)
% 
% Description:
%       This function computes the L2 norm of a dynamical system defined as
%       $|G|_{\mathcal{L}_2} := \sqrt(\frac{1}{2\pi}\int_{-\infty}^\infty)\text{tr}\left(G^H(j\omega)G(j\omega)\right)d\omega$
%
%       For unstable models |sys|, the $\mathcal{H}_2$ norm is undefined.
%       However, if the model has no poles on the imaginary axis, then the
%       L2 norm can be computed. For stable models |sys|, the
%       $\mathcal{H}_2$ and $\mathcal{L}_2$ norms coincide.
%
%       The algorithm is based on [2] and computes $\mathcal{H}_2$ norms of
%       the stable and (time reverted) antistable subsystems.
%
% Input Arguments:
%		*Required Input Arguments:*
%		-sys:   Dynamical system model (ssRed)
%
%
% Output Arguments:
%       -n:     L2 norm
%
% See Also: 
%		norm, ssRed.stabsep
%
% References:
%		* *[1] Antoulas (2005)*, Approximation of large-scale dynamical systems
%		* *[2] Magruder, Beattie, Gugercin (2010)*, Rational Krylov methods
%               for optimal L2 model reduction
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
% Last Change:  02 Aug 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

[sysS,sysA] = stabsep(sys);

if isempty(sysA)
    %system is stable! -> compute H2 norm
    n = norm(sys,2);
else
    %revert time for untistable model
    sysA.(sysA.a_) = -sysA.(sysA.a_); 
    sysA.(sysA.b_) = -sysA.(sysA.b_);
    n = sqrt(norm(sysS,2)^2+norm(sysA,2)^2);
end