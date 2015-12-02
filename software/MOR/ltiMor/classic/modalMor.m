function [sysr, V, W] = modalMor(sys, q, Opts)
% MODALMOR - Modal model order reduction of LTI systems
%
% Syntax:
%       sysr			= MODALMOR(sys, q)
%       sysr			= MODALMOR(sys, q, Opts)
%       [sysr, V, W]	= MODALMOR(sys,... )
%
% Description:
%       This function computes the reduced order system sysr and the 
%       projection matrices V and W by the modal reduction technique [1,4].
%         
%       Only a few eigenvalues and left and right eigenvectors of the pair (A,E) 
%       are computed with the eigs command using the (default) option 'SM' (smallest
%       magnitude). The right eigenvectors build the columns of V, while the 
%       left eigenvectors build the columns of W.
%
%       Depending on the options, the vectors composing the projection
%       matrices are not the left and right eigenvectors, but linear
%       combinations to keep the reduced system matrices real and/or
%       orthogonalize the projection matrices.
%
% Input Arguments:
%		*Required Input Arguments:*
%		-sys:			an sss-object containing the LTI system
%		-q:				order of reduced system
%
%		*Optional Input Arguments:*
%		-Opts:			a structure containing following options
%			-.type:		option to eigs command;
%						[{'SM'} / 'LM' / 'SA' / 'LA' / 'SR' / 'LR' / real or complex scalar]
%			-.orth:		orhtogonalization;
%						[{'0'} / 'qr']
%			-.real:		real reduced system;
% 						[{'0'} / 'real']
%
% Output Arguments:
%       -sysr:          reduced system
%       -V, W:          projection matrices
%
% Examples:
%		This code loads the MIMO benchmark model 'CDplayer' and produces a
%		reduced model that preserves the 10 eigenvalues with smallest
%		magnitude
%
%> sys = loadSss('CDplayer');
%> sysr = modalMor(sys,10);
%> norm(cplxpair(eig(sysr))-cplxpair(eigs(sys,10,'sm')))
%
% See Also: 
%		tbr, rk, irka, eigs
%
% References:
%		* *[1] Antoulas (2005)*, Approximation of large-scale dynamical systems
%		* *[2] Lehoucq and Sorensen (1996)*, Deflation Techniques for an 
%              Implicitly Re-Started Arnoldi Iteration.
%		* *[3] Sorensen (1992)*, Implicit Application of Polynomial Filters 
%              in a k-Step Arnoldi Method.
%		* *[4] Foellinger (2013)*, Regelungstechnik (pp. 305-319)
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Sylvia Cremer, Rudy Eid, Alessandro
%               Castagnotto, Lisa Jeschek
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Dec 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% Default execution parameters
Def.type = 'SM'; 
Def.orth = '0'; %orthogonalization ('0','qr')
Def.real = 'real'; %real reduced system ('0', 'real')

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%calculate right eigenvalues
[V, rlambda] = eigs(sys.A, sys.E, q, Opts.type);
rlambda=diag(rlambda);

%sort eigenvalues and right eigenvectors 
V=V';
tbl=table(rlambda,V);
tbl=sortrows(tbl);
V=tbl.V';
rlambda=tbl.rlambda;

%calculate left eigenvalues

W = zeros(size(V)); llambda = zeros(size(rlambda));
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:eigs:SigmaNearExactEig');
opts.tol=1e3*eps/norm(full(sys.A));
for iEig = 1:length(rlambda)
    [W(:,iEig),llambda(iEig)] = eigs((sys.A-(rlambda(iEig))*sys.E).',1,1e3*eps, opts);
    %Ritz residual: r=Ax-lambda*x
    if (sys.A.'*W(:,iEig)-rlambda(iEig)*sys.E.'*W(:,iEig))/norm(W(:,iEig)) >= norm(full(sys.A))*1e3*eps
        error('Eigenvectors belong to different eigenvalues. Please try again.');
    end
end
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:eigs:SigmaNearExactEig');

%split complex conjugated columns into real and imaginary
if strcmp(Opts.real,'real');
    fac = 2/sqrt(2); %factor needed to guarantee preservation of eig
    k=find(imag(rlambda));
    if ~isempty(k)
         if mod(length(k),2)==0
              for i=1:2:length(k)
                   temp=V(:,i);
                   V(:,i)=fac*real(temp);
                   V(:,i+1)=fac*imag(temp);
                   temp=W(:,i);
                   W(:,i)=fac*real(temp);
                   W(:,i+1)=fac*imag(temp);
              end
         else
              warning('Reduced system contains complex elements. Please try reduction order q+1.');
         end
    end
end

%orthogonalize
if strcmp(Opts.orth,'qr');
    V=qr(V,0);
    W=qr(W,0);
end

sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);   