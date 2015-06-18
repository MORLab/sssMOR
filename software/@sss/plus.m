function sum = plus(sys1, sys2)
% Computes sum of two LTI systems: u-->(sys1-sys2)-->y
% ------------------------------------------------------------------
% sum = sys1+sys2; sum = plus(sys1, sys2);
% Inputs:       * sys1, sys2: summand sss-object
% Outputs:      * sss-object representing sys1+sys2
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de)
% Last Change:  02 Mar 2011
% ------------------------------------------------------------------

if sys1.n == 0
    sum = sss(sys2.A, sys2.B, sys2.C, sys2.D, sys2.E);
    return
end
if sys2.n == 0
    sum = sss(sys1.A, sys1.B, sys1.C, sys1.D, sys1.E);
    return
end
if sys1.p ~= sys2.p
    error('sys1 and sys2 must have same number of inputs.')
end
if sys1.m ~= sys2.m
    error('sys1 and sys2 must have same number of outputs.')
end

sum = sss([sys1.A sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.A], ...
          [sys1.B; sys2.B], ...
          [sys1.C, sys2.C], ...
          sys1.D + sys2.D, ...
          [sys1.E sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.E]);
