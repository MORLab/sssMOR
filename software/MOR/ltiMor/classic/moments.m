function m = moments(sys, s0, n)
% Returns the moments or Markov parameters of an LTI system
% ------------------------------------------------------------------
% m = moments(sys, s0, n)
% Inputs:       * sys: an sss-object containing the LTI system
%               * s0: expansion point (inf -> Markov parameters)
%               * n: number of moments to be computed
% Output:       * m: vector of moments / Markov parameters
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de)
% Last Change:  03 Feb 2011
% ------------------------------------------------------------------

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
