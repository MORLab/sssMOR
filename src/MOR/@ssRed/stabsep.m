function [sysS, sysAs, Vs, Ws, Vas, Was] = stabsep(sys)
% STABSEP - Stable-unstable decomposition
% 
% Syntax:
%		sysS                            = STABSEP(sys)
%		[sysS, sysAs]                   = STABSEP(sys)
%		[sysS, sysAs, Vs, Ws, Vas, Was]	= STABSEP(sys)
% 
% Description:
%       This function return the decomposition of a dynamic system
%       (|ssRed| object) into a stable and antistable part such that
%       |sys = sysS + sysAs|. Eigenvalues on the imaginary axis will be
%       attributed to the antistable part.
%
%       In addition, it returns the projection matrices |Vs|, |Ws| and 
%       |Vas|, |Was| used to obtain |sysS| and |sysAs| from a projection 
%       of |sys| (cp. |projectiveMor|).
%
% Input Arguments:
%		*Required Input Arguments:*
%		-sys:           ssRed object with original model
%
% Output Arguments:
%       -sysS:      Stable subsystem 
%       -sysAs:     Antistable subsystem
%       -Vs,Ws:     projection matrices wrt. stable subsystem
%       -Vas,Was:   projection matrices wrt. antistable subsystem
%
% Examples:
%		This code generate a synthetic unstable ssRed object
%		with eigenvalues left, right, and on the imaginary axis:
%
%> A    = diag([1 0 -1]); b = ones(3,1); c = ones(1,3); 
%> sys  = ssRed(A,b,c);
%> disp(sys)
%> figure; plot(complex(eig(sys)),'bx');hold on
%
%       Using |stabsep|, we can separate the asymptotically stable from the
%       antistable part:
%
%> [sysS, sysAs] = stabsep(sys);
%> disp(sysS), disp(sysAs)
%
%       We can compare the eigenvalues of the subsystem to the eigenvalues
%       of the original model:
%
%> plot(complex(eig(sysS)),'go'); plot(complex(eig(sysAs)),'ro');
%> legend('original','stable','antistable')
%
% See Also: 
%		stabsep, projectiveMor, modalMor
%
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto, Siyang Hu
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sssMOR">www.rt.mw.tum.de/?sssMOR</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Aug 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

    %% Eigendecomposition
    if ~sys.isSym
        [V,D,W] = eig(sys);
    else
        [V,D]   = eig(sys);
        W = V;
    end
    
    %% Find unstable eigenvalues
    iAs = find(real(diag(D))>=0);
       
    %% Perform decomposition
    if isempty(iAs)
        % System is stable
        sysS    = sys;
        sysAs   = ssRed([]);
        Vs      = eye(sys.n);
        Ws      = Vs;
        Vas     = [];
        Was     = [];
    else
        % System is unstable
        iS  = 1:sys.n; iS(iAs) = [];
        
        Vs  = V(:,iS);  Ws  = W(:,iS);        
        Vas = V(:,iAs); Was = W(:,iAs);
        
        %split complex conjugated columns into real and imaginary
        [Vs, Ws]   = makeReal(Vs,   Ws);
        [Vas, Was] = makeReal(Vas,  Was);
        
        %orthogonalize
        [Vs,~]  = qr(Vs, 0);    [Ws,~]  = qr(Ws, 0);
        [Vas,~] = qr(Vas,0);    [Was,~] = qr(Was,0);
        
        sysS    = projectiveMor(sys,Vs, Ws);
        sysAs   = projectiveMor(sys,Vas,Was);
        
        %Update information of ssRed object
        params.originalOrder    = sys.n;
        params.reducedOrder     = sysS.n;
        sysS = ssRed(sysS.(sys.a_),sysS.(sys.b_), ...
                             sysS.(sys.c_),sysS.(sys.d_), ...
                             sysS.(sys.e_),'stabsep', ...
                             params, sys);
                         
        params.reducedOrder     = sysAs.n;
        sysAs = ssRed(sysAs.(sys.a_),sysAs.(sys.b_), ...
                             sysAs.(sys.c_),sysAs.(sys.d_), ...
                             sysAs.(sys.e_),'stabsep', ...
                             params, sys);

        % make robust to computations with sys.e_
        if isempty(sysS.(sys.e_)), sysS.(sys.e_)    = eye(sysS.n); end
        if isempty(sysAs.(sys.e_)), sysAs.(sys.e_) 	= eye(sysAs.n); end
    end
end
    
%% ------------------ AUXILIARY FUNCTIONS --------------------------
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