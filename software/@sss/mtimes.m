function prod = mtimes(sys2, sys1)
% Computes product of two LTI systems: u --> sys2 --> sys1 --> y
% ------------------------------------------------------------------
% prod = sys1*sys2; prod = mtimes(sys1, sys2);
% Inputs:       sys1, sys2: factor sss-objects
% Outputs:      * sss-object representing sys1-sys2
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de)
% Last Change:  28 Nov 2013 Yvonne Stuerz
% ------------------------------------------------------------------

if (sys1.n==0 && isempty(sys1.D)) || (sys2.n==0 && isempty(sys2.D))
    prod = sss([], [], []);
    return
end
if sys2.p ~= sys1.m
    error('Number of outputs of sys1 does not match number of inputs of sys2.')
end

prod = sss([sys1.A sparse(sys1.n,sys2.n); sys2.B*sys1.C sys2.A], ...
          [sys1.B; sys2.B*sys1.D], ...
          [sys2.D*sys1.C, sys2.C], ...
          sys2.D*sys1.D, ...
          [sys1.E sparse(sys1.n,sys2.n); sparse(sys2.n,sys1.n) sys2.E]);
