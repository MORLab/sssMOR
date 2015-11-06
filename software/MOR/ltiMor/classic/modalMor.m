function [sysr, V, W] = modalMor(sys, q, Opts)
% MODALMOR - Modal model order reduction of LTI SISO systems
%
% Syntax:
%       [sysr, V, W] = MODALMOR(sys, q, Opts)
%
% Description:
%       This function computes the reduced order system sysr and the 
%       projection matrices V and W by the modal reduction technique [1].
%         
%       Only a few eigenvalues and left and right eigenvectors of the pair (A,E) 
%       are computed with the eigs command. The right eigenvectors build the
%       columns of V, while the left eigenvectors build the columns of W.
%
% Input Arguments:
%		*Required Input Arguments:*
%		-sys:           an sss-object containing the LTI system
%		-q:             order of reduced system
%
%		*Optional Input Arguments:*
%		-Opts:	 		a structure containing following options
%			-.type:  	option to eigs command; [{'SM'} / 'LM' / 'SA' / 'LA' / 'SR' / 'LR' / real or complex scalar]
%			-.orth:  	orhtogonalization; [{'0'} / 'qr']
%			-.real:  	real reduced system; [{'0'} / 'real']
%
% Output Arguments:
%       -sysr:          reduced system
%       -V, W:          projection matrices
%
% Examples:
%		TODO
%
% See Also: 
%		tbr, rk, irka
%
% References:
%		* *[1] Antoulas (2005)*, Approximation of large-scale dynamical systems
%		* *[2] Lehoucq and Sorensen (1996)*, Deflation Techniques for an Implicitly Re-Started Arnoldi Iteration.
%		* *[3] Sorensen (1992)*, Implicit Application of Polynomial Filters in a k-Step Arnoldi Method.
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
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer, Rudy Eid
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  11 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% Default execution parameters
Def.type = 'SM'; 
Def.orth = '0'; %orthogonalization ('0','qr')
Def.real = '0'; %real reduced system ('0', 'real')

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

% note: eigenvector direction corresponds to the nullspace of the shifted 
%       pencil (sys.A-rlambda(iEig)*sys.E)
%       Since the sparse algorithm does not allow to input the exact value,
%       we will perturb it a little

% W = zeros(size(V)); warning('off','MATLAB:nearlySingularMatrix');
% tic
% for iEig = 1:length(rlambda)
%     [~,~,W(:,iEig)] = svds((sys.A-(rlambda(iEig)+1e2*eps)*sys.E).',1,0);
% end
% t = toc
% warning('on','MATLAB:nearlySingularMatrix');

% this following code should be a little more accurate (it looks for the
% right eigenvalue at 0) and seems to be a little faster

W = zeros(size(V)); llambda = zeros(size(rlambda));
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:eigs:SigmaNearExactEig');
for iEig = 1:length(rlambda)
    [W(:,iEig),llambda(iEig)] = eigs((sys.A-(rlambda(iEig))*sys.E).',1,1e3*eps);
end
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:eigs:SigmaNearExactEig');

%check if the same eigenvalues have been found
% in this case it is equivalent to verifying that llambda is close to 0
if norm(llambda)/length(llambda) > 1e3*eps
     error('Eigenvectors belong to different eigenvalues. Please try again.');
end


%   OLD CODE
% [W, llambda] = eigs(sys.A', sys.E', q, Opts.type);
% llambda=diag(llambda);

% %sort eigenvalues and left eigenvectors
% W=W';
% tbl=table(llambda,W);
% tbl=sortrows(tbl);
% W=tbl.W';
% llambda=tbl.llambda;
% 
% %check if the same eigenvalues have been found
% for i=1:q
%     if abs(real(llambda(i))-real(rlambda(i)))+abs(imag(llambda(i))-imag(rlambda(i))) > 10e-6
%         error('Eigenvectors belong to different eigenvalues. Please try again.');
%     end
% end

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