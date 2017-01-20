function [sysr, V, W, D] = modalMor(sys, q, Opts)
% MODALMOR - Modal truncation order reduction of LTI systems
%
% Syntax:
%       sysr			= MODALMOR(sys, q)
%       sysr			= MODALMOR(sys, q, Opts)
%       [sysr, V, W, D]	= MODALMOR(sys,... )
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
%       A dominance analysis of the eigenvalues can be performed.
%       In addition, the system can be reduced so that the most dominant 
%       eigenvalues of the chosen options are preserved.
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
%						[false / {'qr'}]
%			-.real:		real reduced system;
% 						[{'real'} / '0']
%           -.tol:      tolerance for the eigenspace computation
%                       [{1e-6} / positive float]
%           -.dominance: perform dominance analysis
%                       [{0} / 'analyze' / '2q' / '3q' / '4q']
%           -.subspaceW:choose how the left invariant subspace W is computed;
%                       [{'eigs'} / '1by1' ]
%           -.lse:      solve linear system of equations
%                       [{'sparse'} / 'full']
%
% Output Arguments:
%       -sysr:          reduced system
%       -V, W:          projection matrices
%       -D:             dominance values
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
%       * *[5] Litz and Roth (1981)*, State decomposition for singular perturbation
%              
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
% Authors:      Heiko Panzer, Sylvia Cremer, Rudy Eid, 
%               Alessandro Castagnotto, Lisa Jeschek
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  23 Nov 2016
% Copyright (c) 2015,2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% Default execution parameters
Def.type        = 'SM'; 
Def.orth        = 'qr'; %orthogonalization (false,true,'qr')
Def.real        = 'real'; %real reduced system (true, false)
Def.tol         = 1e-6; % tolerance for SM/LM eigenspace
Def.dominance   = 0; %dominance analysis ('0','analyze','2q','3q',..,'9q')
Def.lse         = 'sparse'; % solveLse ('sparse', 'full')
Def.subspaceW   = 'eigs'; % 'eigs', '1by1': choose how to compute W

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

if nnz(Opts.dominance) && ~strcmp(Opts.dominance,'analyze')
    q=str2double(Opts.dominance(1))*q;
    if q>size(sys.A,1)
        q=size(sys.A,1);
    end
end

%% Invariant subspace computation
%compute right and left eigenvectors
if q>=size(sys.A,1)-1 %eigs: q<n-1
    if ~sys.isSym
        [V,~,W]=eig(full(sys.A),full(sys.E));
        V=V(:,1:q);
        W=W(:,1:q);
    else
        [V,~]=eig(full(sys.A),full(sys.E));
        V=V(:,1:q);
        W=V;
    end
else
    [V, W]=eigenspace(q);
end


%% Post processing
%dominance analysis
if Opts.dominance
    [V, W, D] = dominanceAnalysis(q, V, W);
end

%split complex conjugated columns into real and imaginary
if strcmp(Opts.real,'real')
    [V, W] = makeReal(V, W);
end

%orthogonalize
if strcmp(Opts.orth,'qr')
    [V,~]=qr(V,0);
    [W,~]=qr(W,0);
end

sysr = projectiveMor(sys,V,W);
if ~ischar(Opts.type), Opts.type = num2str(Opts.type); end
if ~ischar(Opts.dominance), Opts.dominance = num2str(Opts.dominance); end

%%  Storing additional parameters
%Stroring additional information about thr reduction in the object 
%containing the reduced model:
%   1. Define a new field for the Opts struct and write the information
%      that should be stored to this field
%   2. Adapt the method "parseParamsStruct" of the class "ssRed" in such a
%      way that the new defined field passes the check
Opts.originalOrder = sys.n;

% Convert the reduced system to an ssRed-object
sysr = ssRed('modalMor',Opts,sysr);
sysr.Name = sprintf('%s_%i_modal_%s',sys.Name,sysr.n,Opts.type);

%% ------------------ AUXILIARY FUNCTIONS --------------------------
function [V, W, rlambda] = eigenspace(q)
    %calculate right eigenvectors
    [V, rlambda] = eigs(sys.A, sys.E, q, Opts.type);
    rlambda=diag(rlambda);

    %calculate left eigenvectors
    if sys.isSym
        W=V;
    else
        ceA = condest(sys.A);
        resTol = ceA*1e3*eps;  
        switch Opts.subspaceW
            case 'eigs'
                opts.v0 = real(V(:,1));
                [W,llambda] = eigs(sys.A.', sys.E.', q, Opts.type, opts);
                llambda=diag(llambda);
                % check if we converged to same eigenvalues
                if norm(sort(llambda)-sort(rlambda))/norm(rlambda) > resTol
                    warning('sssMOR:modalMor:eigsConvergence',...
                        'Eigs did not converge to same eigenvalues for V and W. Consider setting Opts.subspaceW = ''exact''.')
                end
            case '1by1'
                W = zeros(size(V));
                warning('off','MATLAB:nearlySingularMatrix');
                warning('off','MATLAB:eigs:SigmaNearExactEig');
                opts.tol=1e3*eps/ceA;
                for iEig = 1:length(rlambda)                    
                    %check if eigenvalue is complex conjugated
                    if iEig==1 || isreal(rlambda(iEig))
                        cplxconjpartner = false;
                    else %complex eigenvalue
                        if abs(rlambda(iEig-1)-conj(rlambda(iEig)))>1e-6
                            cplxconjpartner = true;
                        else
                            cplxconjpartner = false;
                        end
                    end
                        
                    if cplxconjpartner
                        %switch order of complex conjugated W -> W'*A*V diagonal
                        W(:,iEig)   = W(:,iEig-1);
                        W(:,iEig-1) = conj(W(:,iEig-1)); 
                    else
                        opts.v0 = V(:,iEig);
                        [W(:,iEig),~, flag] = eigs((sys.A-(rlambda(iEig))*sys.E).',1,1e3*eps, opts);
                        %Ritz residual: r=Ax-lambda*E*x
                        if flag~=0
                            warning('sssMOR:modalMor:eigsConvergence','Eigs did not converge for selected eigenvalue.');
                        elseif norm(sys.A.'*W(:,iEig)-rlambda(iEig)*sys.E.'*W(:,iEig))/norm(W(:,iEig)) >= resTol
                            error('sssMOR:modalMor:ritzResidual','Eigenvectors in W and V belong to different eigenvalues.');
                        end
                    end
                end
                warning('on','MATLAB:nearlySingularMatrix');
                warning('on','MATLAB:eigs:SigmaNearExactEig');
        end
    end
end
function [V, W] = makeReal(V, W)
    fac = 2/sqrt(2); %factor needed to guarantee preservation of eig
    k=find(any(imag(V)));
    if ~isempty(k)
         if mod(length(k),2)==0 %make sure they come in pairs
              for iK = 1:2:length(k)
                idx = k(iK);
                temp=V(:,idx);
                V(:,idx)=fac*real(temp);
                V(:,idx+1)=fac*imag(temp);
                temp=W(:,idx);
                W(:,idx)=fac*real(temp);
                W(:,idx+1)=fac*imag(temp);
              end
         else
              warning('Reduced system contains complex elements. Please try reduction order q+1.');
         end
    end
end
function [V, W, D] = dominanceAnalysis(q, V, W)  
    D=zeros(q,1);
    tempSys=sys.B*sys.C;
    for j=1:q
        D(j)=norm((W(:,j)'*tempSys*V(:,j))/(W(:,j)'*sys.A*V(:,j)));
    end
    if ~strcmp(Opts.dominance,'analyze')
        %take eigenvectors of most dominant eigenvalues
        T=speye(q);
        tbl=table(-D,T);
        tbl=sortrows(tbl);
        q=q/str2double(Opts.dominance(1));
        D=-tbl.Var1(1:q);
        V=V*tbl.T';
        W=W*tbl.T';
        V=V(:,1:q);
        W=W(:,1:q);
    end 
end
end