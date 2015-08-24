function sys = fem2ss(M, D, K, B, Cx, Cv, F)
% Converts from 2nd order system to state space representation
% ------------------------------------------------------------------
% sys = fem2ss(M, D, K, B, Cx, Cv, F)
% Inputs:       * M, D, K, B, Cx, Cv: Matrices of Second Order System
%    [optional] * F: regular matrix establishing equality of x' (see [1])
% Outputs:      * sys: Sparse State Space Model
% ------------------------------------------------------------------
% The matrices of the second order system are defined as:
%    M x" + D x' + K x = B u
%                    y = Cx x + Cv x'
%
% D and either Cx or Cv may be empty; they are then set to zero.
%
% [1] B. Salimbahrami: "Structure Preserving Order Reduction of Large
%     Scale Second Order Models", PhD Thesis at TUM, p. 36f
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Rudy Eid
% Last Change:  23 Feb 2011
% ------------------------------------------------------------------

if size(M,1)~=size(M,2)
    error('M must be symmetric.')
elseif any(size(M)-size(K))
    error('M and K must have same size.');
end

if isempty(D) || ~any(any(D))
    D = sparse(size(M,1), size(M,1));
elseif any(size(D)-size(K))
    error('D must have same size as M and K.');
end

if ~exist('F', 'var') || isempty(F)
    % if no F is given, use F:=K to obtain E>0 and A+A'<=0
    F=K;
elseif any(size(F)-size(K))
    error('F must have same size as M and K.');
end

n = size(M,1);
O = sparse(n,n);

E = [F O; O M];
A = [O F; -K -D];

if isempty(Cx) && isempty(Cv)
    % both output vectors are empty
    error('All outputs are zero.');
elseif isempty(Cx)
    % only Cv is given
    Cx = sparse(size(Cv,1),n);
elseif ~exist('Cv', 'var') || isempty(Cv)
    % only Cx is given
    Cv = sparse(size(Cx,1),n);
elseif any(size(Cx)-size(Cv))
    error('Cx and Cv must have same size.');
elseif size(Cx,2)~=size(M,2)
    error('Cx, Cv must have same column dimension as M, D, K')
end

if isempty(B)
    error('All inputs zero.');
elseif size(B,1)~=size(M,1)
    error('B must have same row dimension as M, D, K')
end

B = [sparse(n,size(B,2)); B];
C = [Cx, Cv];

sys = sss(A, B, C, 0, E);
