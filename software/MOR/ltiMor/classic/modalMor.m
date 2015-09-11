function [sysr, V, W] = modalMor(sys, q, Opts)
% MODALMOR - Modal model order reduction of LTI SISO systems
% ------------------------------------------------------------------
% [sysr, V, W] = MODALMOR(sys, q, Opts)
% Inputs:       * sys: an sss-object containing the LTI system
%               * q: order of reduced system
%               * opts: a structure containing following options
%                   * type: options to eigs command - {'SM','LM',...}
%                       - 'SM' for eigenvalues of Smallest Magnitude
%                       - 'LM' for eigenvalues of Largest Magnitude
%                       - scalar number for eigenvalues next to it
% Outputs:      * sysr: reduced system
%               * V, W: projection matrices
% ------------------------------------------------------------------
% USAGE:  This function computes the reduced order system sysr and the 
% projection matrices V and W by the modal reduction technique [1].
%
% Only a few eigenvalues and left and right eigenvectors of the pair (A,E) 
% are computed with the eigs command. The right eigenvectors build the
% columns of V, while the left eigenvectors build the columns of W.
%
% See also TBR, RK, IRKA
%
% ------------------------------------------------------------------
% REFERENCES:
% [1] Antoulas (2005), Approximation of large-scale dynamical systems
% [2] Lehoucq and Sorensen (1996), Deflation Techniques for an Implicitly 
% Re-Started Arnoldi Iteration.
% [3] Sorensen (1992), Implicit Application of Polynomial Filters in a 
% k-Step Arnoldi Method.
% ------------------------------------------------------------------
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                   -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer, Rudy Eid
% Last Change:  11 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

% Default execution parameters
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
   