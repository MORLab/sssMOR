function m = moments(sys, s0, n)
% MOMENTS - Returns the moments or Markov parameters of an LTI system
%
% Syntax:
%       m = moments(sys, s0, n)
%
%
% Inputs:
%      -sys:    an sss-object containing the LTI system
%      -s0:     expansion point (inf -> Markov parameters)
%      -n:      number of moments to be computed
%
%
% Output: 
%      -m:     vector of moments / Markov parameters
%
%
% Examples:
%       No examples
% 
% 
%------------------------------------------------------------------
%   This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State Space, Model Order 
%   Reduction and System Analysis Toolbox developed at the Chair of 
%   Automatic Control, Technische Universitaet Muenchen. For updates 
%   and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
%   More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  03 Feb 2011
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

m=zeros(1,n);
if isinf(s0)
    m(1) = sys.D;
    [L,U,k,l,S]=lu(sys.E, 'vector');
    B=sys.B;
    for i=2:n
        b=S\B; b=b(k,:);
        x=L\b;
        x(l,:)=U\x;
        m(i) = sys.C*x - sys.D;
        B=sys.A*x;
    end
else
    [L,U,k,l,S]=lu(sys.A - s0*sys.E, 'vector');
    B=sys.B;
    for i=1:n
        b=S\B; b=b(k,:);
        x=L\b;
        x(l,:)=U\x;
        m(i) = sys.C*x - sys.D;
        B=sys.E*x;
    end
end
