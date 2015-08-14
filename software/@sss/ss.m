function sys_ss = ss(sys_sss)
% Converts sparse LTI system (sss) to Matlab\control\ss
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% sys_ss = ss(sys_sss);
% Input:        * sys_sss: sparse state space (sss)-object
% Output:       * sys_ss:  ss- or dss-object, respectively
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer
% Last Change:  02 Feb 2011
% ------------------------------------------------------------------

if sys_sss.isdescriptor
    sys_ss=dss(full(sys_sss.A),full(sys_sss.B),full(sys_sss.C),full(sys_sss.D),full(sys_sss.E));
else
    sys_ss=ss(full(sys_sss.A),full(sys_sss.B),full(sys_sss.C),full(sys_sss.D));
end