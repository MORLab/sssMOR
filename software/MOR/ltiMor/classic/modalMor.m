function [sysr, V, W] = modalMor(sys, q, Opts)
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

%   Default execution parameters
Def.type = 'SM'; 

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

if sys.isdescriptor
    E=sys.E;
else
    E=[];
end

[V, p] = eigs(sys.A, E, q, Opts.type);
p=diag(p);
W=zeros(size(V));
r=zeros(size(p));
for i=1:q
    warning off
    sigma = p(i) - 1e-6; %eigs fails if sigma is equal to an eigenvalue
    [wi, pi]=eigs(sys.A',sys.E', 1, sigma); %#ok<NASGU>
    warning on
    W(:,i) = wi;
    r(i)=pi;
end

%sort eigenvalues and eigenvectors in descending/ascending order 
%(complex conjugated pairs are sorted successively)
V=V';
tbl=table(p,V);
tbl=sortrows(tbl);
V=tbl.V';
p=tbl.p;

W=W';
tbl=table(r,W);
tbl=sortrows(tbl);
W=tbl.W';
r=tbl.r;

%split complex conjugated columns into real and imaginary
k=find(imag(p));
if ~isempty(k)
     if mod(length(k),2)==0
          for i=1:2:length(k)
               temp=V(:,i);
               V(:,i)=real(temp);
               V(:,i+1)=imag(temp);
          end
     else
          warning('Reduced system contains complex elements. Please try the reduction order q+1.');
     end
end

k=find(imag(r));
if ~isempty(k)
     if mod(length(k),2)==0
          for i=1:2:length(k)
               temp=W(:,i);
               W(:,i)=real(temp);
               W(:,i+1)=imag(temp);
          end
     else
          warning('Reduced system contains complex elements. Please try the reduction order q+1.');
     end
end

V=orth(V);
W=orth(W);

sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);
   