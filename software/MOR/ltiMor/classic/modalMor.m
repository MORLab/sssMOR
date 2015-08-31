function [sysr, V, W] = modalMOR(sys, q, varargin)
% Modal order reduction of LTI SISO systems
% ------------------------------------------------------------------
% [sysr, V, W] = modalMOR(sys, q,handles)
% Inputs:       * sys: an sss-object containing the LTI system
%               * q: order of reduced system
%    [optional] * options: options to eigs command:
%                  - 'SM' for eigenvalues of Smallest Magnitude
%                  - 'LM' for eigenvalues of Largest Magnitude
%                  - scalar number for eigenvalues next to it
% Outputs:      * sysr: reduced system
%               * V, W: projection matrices
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer, Rudy Eid
% Last Change:  11 Feb 2011
% ------------------------------------------------------------------

if sys.is_dae
    E=sys.E;
else
    E=[];
end
if nargin==2
    options = {'SM'};
else
    options = varargin;
end

[V, p] = eigs(sys.A, E, q, options{:});
p=diag(p);
W=zeros(size(V));
for i=1:q
    warning off
    [wi, pi]=eigs(sys.A',sys.E', 1, p(i)); %#ok<NASGU>
    warning on
    W(:,i) = wi;
end

% ***
V=[real(V(:,1:2:q)) imag(V(:,1:2:q))];
V=orth(V);
W=[real(W(:,1:2:q)) imag(W(:,1:2:q))];
W=orth(W);
sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);
   