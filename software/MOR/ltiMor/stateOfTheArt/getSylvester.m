function [R, S, B_] = getSylvester(sys,sysr,V,type)
% GETSYLVESTER - get matrices of Sylvester equation
%
% Syntax: 
%       [R, S, B_] = GETSYLVESTER(sys, sysr, V, type)
%
% Description:
%       Given the full order model sys and a reduced order model sysr
%       obtained through a projection with V being an input Krylov
%       subspace, this function reconstructs the matrices of the Sylvester
%       equation for the Krylov subspace
%           A V - E V S - B R = 0           (1)
%           A V - E V Er^-1 Ar - B_ R = 0   (2)
%
%       If the type is set to 'W', then the matrices are given for the
%       output Krylov Sylvester equation
%           A.' W - E.' W S.' - C.' R = 0           (3)
%           A.' W - E.' V Er^-T Ar.' - C_.' R = 0   (4)
%
% 
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:      full order model
%       -sysr:     reduced order model
%       -V:        input Krylov subspace
%       *Optional Input Arguments:*
%       -type:     {'V' (def), 'W'} specifies if V spans an input (def) or 
%                  output Krylov subspace
%
% Output Arguments: 
%       -R,S:      matrices of Sylvester equation (1) or (3)
%       -B_:       matrix of Sylvester equation (2) or (4)
% 
% Examples:
%       TODO
% 
% See Also: 
%       porkV, porkW, cure, spark
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  08 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  parsing
if nargin == 3
    type = 'V';
end
if strcmp(type,'W')
    % use the dual system
    sys = sys.'; sysr = sysr.';
end

%%  computations
B_ = sys.B - sys.E*V*(sysr.E\sysr.B);
R = (B_.'*B_)\(B_.'*(sys.A*V - sys.E*V*(sysr.E\sysr.A)));
if nargout > 1
    S = sysr.E\(sysr.A - sysr.B*R);
    if strcmp(type,'W')
        S = S.';
    end
end

%% control the accuracy by computing the residual
res = zeros(1,3);
res(1) = norm(sysr.A - sysr.E*S - sysr.B*R);
res(2) = norm(sys.A*V - sys.E*V*S - sys.B*R);
res(3) = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A)-B_*R);

if any( res > 1e-6 )
    warning('careful, the problem might be ill conditioned and the results of getSylvester inaccurate');
end