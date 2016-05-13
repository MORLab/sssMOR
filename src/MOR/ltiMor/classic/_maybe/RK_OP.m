function alpha_opt = RK_OP(sys)
% Determination of optimal expansion point for Laguerre series
% ------------------------------------------------------------------
% aopt = RK_OP(sys)
% Input:        * sys: an sss-object containing the LTI system
% Output:       * alpha_opt: Optimal Laguerre parameter
% ------------------------------------------------------------------
% For details and extensibilities see the references:
% [1] R. Eid: Time domain Model Reduction by Moment Matching, Ph.D thesis,
%   Institute of Automatic Control, Technische Universitaet Muenchen, 2009.
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Rudy Eid
% Last Change:  03 Feb 2011
% ------------------------------------------------------------------

A = full(sys.E\sys.A); B = full(sys.E\sys.B); C = full(sys.C);
B=B(:,1); C=C(1,:);

try
    P=lyap(A,B*B');
    Y=lyap(A,P);
    alpha_opt=sqrt((C*A*Y*A'*C')/(C*Y*C'));
    if alpha_opt<0 || imag(alpha_opt)~=0
        disp('Warning. RK_OP yields negative or complex result.')
        alpha_opt = abs(alpha_opt);
    end
catch ex
    disp(['Error during calculation of RK_OP: ' ex.message]);
    alpha_opt = NaN;
end
