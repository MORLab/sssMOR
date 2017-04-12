function sysD = shiftD(sys,D,L,R)
% SHIFTD - Add a feedthrough D while preserving tangential interpolation
% 
% Syntax:
%		sysD				= SHIFTD(sys,D,L,R)
% 
% Description:
%       This function adds a feedthrough term |D| to |sys| such that
%           sysD.D = sys.D + D
%
%       All other system matrices are modified accordingly such that Gd(s)
%       and G(s) share the same bi-tangential data specified by the Sylvester
%       equations 
%           AV  - EVS  - BR  = 0
%           A'W - EWS' - C'L = 0
%
%       For more details, refer to the references.
%           
% Input Arguments:
%		*Required Input Arguments:*
%		-sys: 			state-space model
%		-D: 			additional feed-through
%       -R:             matrix of right tangential directions
%       -L:             matrix of left tangential directions
%						
% Output Arguments:
%       -sysd: modified model with additonal feedthrough and
%                       preserved tangential data 
%
% See Also: 
%		hInfMor
%
% References:
%		* *[1] Mayo, Antoulas (2007)*, A framework for the solution of the generalized realization problem
%		* *[2] Beattie, Gugercin (2009)*, Interpolatory projection methods for structure-preserving model reduction
%		* *[3] Castagnotto et al. (2017)*, Interpolatory methods for Hinf model reduction of multi-input/multi-output systems
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
% Last Change:  06 Apr 2017
% Copyright (c) 2016,2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


switch class(sys)
    case 'sss'
        sysD = sss(sys.A+L.'*D*R, sys.B+L.'*D, sys.C+D*R, D, sys.E);
    case 'ss'
        sysD = dss(sys.A+L.'*D*R, sys.B+L.'*D, sys.C+D*R, D, sys.E);
    case 'sso'
        sysD = sso(sys.M,sys.D,sys.K-L.'*D*R, sys.B-L.'*D, sys.Cp-D*R, sys.Cv, sys.Df+D);
    otherwise
        error('sssMOR:shiftD:notSupportedClass','The sys object provided is not supported')
end
