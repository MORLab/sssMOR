function [varargout] = solveLse(varargin)
% ARNOLDI - Arnoldi algorithm for Krylov subspaces with multiple shifts
% 
% Syntax:
% E,A,V,W,s0,Rt,Lt,Opts
%       X                                = SOLVELSE(A,B)
%       [X,Y]                            = SOLVELSE(A,B,C)
%       X                                = SOLVELSE(A,B,E,s0)
%       X                                = SOLVELSE(A,B,E,s0,Rt)
%       [X,Sx,Rx]                        = SOLVELSE(A,B,E,s0)
%       [X,Sx,Rx]                        = SOLVELSE(A,B,E,s0,Rt)
%       [X,Y]                            = SOLVELSE(A,B,C,E,s0)
%       [X,Y,Sx,Rx,Sy,Ly]                = SOLVELSE(A,B,C,E,s0)
%       [X,Y]                            = SOLVELSE(A,B,C,E,s0,Rt,Lt)
%       [X,Y,Sx,Rx,Sy,Lx]                = SOLVELSE(A,B,C,E,s0,Rt,Lt)
%       [X,Y,Sx,Rx,Sy,Lx]                = SOLVELSE(A,B,C,E,s0,...,IP,Opts)
% 
% Description:
%       This function solves linear systems of equations X=(A-s0*E)\B 
%       corresponding to shifts s0. The order of the shifts is crucial for
%       reusing already computed factorizations, so it is recommended to 
%       sort the shifts in advance. If the matrix E is empty or not 
%       specified, X=A\B is computed. If s0 is Inf, then the Markov 
%       parameter X=(A-s0*E)*B is computed. If Opts.krylov is set, a Krylov
%       subspace [1-3] is created from the solutions.
%
%       In addition, if the output matrix C is passed, then SOLVELSE
%       computes the solutions X=(A-s0*E)\B and Y=C/(A-s0*E)'. If
%       Opts.krylov is specified, this means that input and output Krylov 
%       subspaces corresponding to the same shifts are computed. Optionally
%       this function computes the Sylvester matrices and orthogonalizes 
%       the Krylov subspaces [4,5]. The orthogonalization is conducted with 
%       respect to the inner product defined in IP. For more details, 
%       please refer to arnoldi.
% 
%       In case of MIMO models, matrices of tangential directions Rt 
%       (and Lt) have to be defined. They must have the same number of 
%       columns as the shifts s0, so that for each tangential direction it 
%       is clear to which shift it belongs.
%
%
% Input Arguments:
%       *Required Input Arguments:*
%       -A/B:  System matrices or right/left side matrices
%
%       *Optional Input Arguments:*
%       -s0:                Vector of complex conjuate expansion points
%       -E/C:               System matrices
%       -Rt,Lt:             Matrix of right/left tangential directions
%       -IP:                function handle for inner product
%       -Opts:              a structure containing following options
%           -.real:         keep the projection matrices real
%                           [{true} / false]
%           -.orth:         orthogonalization of new projection direction
%                           [{'2mgs'} / 0 / 'dgks' / 'mgs']
%           -.reorth:       reorthogonalization
%                           [{'gs'} / 0 / 'qr']
%           -.lse:          use LU or hessenberg decomposition
%                           [{'sparse'} / 'full' / 'hess' / 'iterative']
%           -.dgksTol:      tolerance for dgks orthogonalization
%                           [{1e-12} / positive float]
%           -.krylov:       standard or cascaded krylov basis
%                           [{0} / 'standardKrylov' / 'cascadedKrylov']
%           -.maxiterlse:   maximum number of iterations in iterSolve
%                           [{1e3} / positive integer]
%           -.tollse:       residual tolerance in iterSolve
%                           [{1e-6} / positive float]
%           -.solver:       preferred solver in iterSolve
%                           [{'cgs'} / 'bicgstab' / 'bicg']
%           -.verbose:      show warnings?
%                           [{1} / 0]
%           -.force:        force solve iteratively when not converging
%                           [{0} / 1]
%
% Output Arguments:
%       -X:        Lse solution corresp. to B/Orthonormal basis spanning the input Krylov subsp. 
%       -Sv:       Matrix of input Sylvester Eq. (1)
%       -Rv:       Right tangential directions of Sylvester Eq. (1), (mxq) matrix
%       -Y:        Lse solution corresp. to C/Orthonormal basis spanning the output Krylov subsp.
%       -Sw:       Matrix of output Sylvester Eq. (2)
%       -Lw:       Left tangential directions of Sylvester Eq. (2), (pxq) matrix
%
% See Also: 
%       arnoldi, rk, irka, projectiveMor
%
% References:
%       * *[1] Grimme (1997)*, Krylov projection methods for model reduction
%       * *[2] Antoulas (2005)*, Approximation of large-scale dynamical systems
%       * *[3] Antoulas (2010)*, Interpolatory model reduction of large-scale...
%       * *[4] Giraud (2005)*, The loss of orthogonality in the Gram-Schmidt... 
%       * *[5] Daniel (1976)*, Reorthogonalization and stable algorithms...
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
% Authors:      Heiko Panzer, Alessandro Castagnotto, Maria Cruz Varona,
%               Lisa Jeschek
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  19 Jul 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

Def.real = false; %keep the projection matrices real?
Def.orth = 0; %orthogonalization after every direction {0,'dgks','mgs','2mgs'}
Def.reorth = 0; %reorthogonaliation at the end {0, 'mgs', 'qr'}
Def.lse = 'sparse'; %use sparse or full LU or lse with Hessenberg decomposition {'sparse', 'full','hess','iterative'}
Def.dgksTol = 1e-12; %orthogonality tolerance: norm(V'*V-I,'fro')<tol
Def.krylov = 0; %standard or cascaded krylov basis (only for siso) {0,'cascade'}  

Def.solver = 'cgs'; %first iterative solver to try
Def.maxiterlse = 1000; %maximum number of iterations in iterative solver
Def.tollse = 1e-6; %residual tolerance in iterSolve 
Def.verbose = 1; %display warnings when iterative methods fail
Def.force = 0;  %not converging in iterSolve leads to error (0) or warning (1)

%% parsing of inputs
if isa(varargin{end},'struct')
    Opts=varargin{end};
    varargin=varargin(1:end-1);
end

if isa(varargin{end},'function_handle')
    IP=varargin{end};
    varargin=varargin(1:end-1);
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end  

if length(varargin)>1
    A=varargin{1};
    B=varargin{2};

    switch length(varargin)
        case 2
            hermite=false;
        case 3
            if size(varargin{3},1)==size(A,1) && size(varargin{3},2)==size(A,2)
                error('Please specify s0.');
            else
                C=varargin{3};
                hermite=true;
            end
        case 4
            if size(varargin{3},1)==size(A,1) && size(varargin{3},2)==size(A,2)
                E=varargin{3};
                s0=varargin{4};
                hermite=false;
            else
                C=varargin{3};
                E=varargin{4};
                s0=1;
                hermite=true;
            end
        case 5
            if size(varargin{3},1)==size(A,1) && size(varargin{3},2)==size(A,2)
                E=varargin{3};
                s0=varargin{4};
                Rt=varargin{5};
                hermite=false;
            else
                C=varargin{3};
                E=varargin{4};
                s0=varargin{5};
                hermite=true;
            end
        case 7
            C=varargin{3};
            E=varargin{4};
            s0=varargin{5};
            Rt=varargin{6};
            Lt=varargin{7};
            hermite=true;
        otherwise
            error('Wrong inputs');
    end
else
    error('More inputs required.');
end


% check E-matrix, tangential directions and IP
if ~exist('E','var') || isempty(E)
    withoutE=true;
    s0=0;
else
    withoutE=false;
end

if ~exist('IP', 'var') 
   IP=@(x,y) (x'*y); %seems to be better conditioned that E norm
end

if ~exist('Rt','var')
    if size(B,2)==1 %siso
        Rt=ones(1,length(s0));
        if hermite && size(C,1)==1
            Lt=ones(1,length(s0));
        end
    else
        error('Please specify tangential directions.');
    end
end

% If the 'full' option is selected for LU, convert E,A once to full
if withoutE
    if strcmp(Opts.lse,'full')
        A = full(A);
    elseif strcmp(Opts.lse,'hess')
        [P,A] = hess(full(A)); B = P*B; if hermite, C = C*P'; end
    else
        A=sparse(A);
    end
else
    if strcmp(Opts.lse,'full')
        E = full(E); A = full(A);
    elseif strcmp(Opts.lse,'hess')
        [A,E,Q,Z] = hess(full(A),full(E)); B = Q*B; if hermite, C = C*Z; end
    else
        E = sparse(E); A=sparse(A);
    end
end

% preallocate memory
q=length(s0)+nnz(imag(s0));
V=zeros(length(B),q);
Rv=zeros(size(B,2),q);
Sv=zeros(q);
if hermite 
    W = zeros(length(B),q); 
    Lw = zeros(size(C,1),q);
    Sw=zeros(q);
end

for jCol=1:length(s0)
    if hermite
        [V, SRsylv, Rsylv, W, SLsylv, Lsylv] = nextDirection(jCol, s0, V, W);
    else
        [V, SRsylv, Rsylv] = nextDirection(jCol, s0, V);
    end
    Sv(:,jCol) = SRsylv;
    Rv(:,jCol) = Rsylv*Rt(:,jCol);
    if hermite
        Sw(jCol,:) = SLsylv.';
        Lw(:,jCol) = Lsylv*Lt(:,jCol);
    end

    % split complex conjugate columns into real (->j) and imag (->j+length(s0c)/2
    if Opts.real
        if hermite
            [V, Sv, Rv, W, Sw, Lw] = realSubspace(jCol, q, s0, V, Sv, Rv, W, Sw, Lw);
        else
            [V, Sv, Rv] = realSubspace(jCol, q, s0, V, Sv, Rv);
        end
    end

    if Opts.orth
        if hermite
            [V, TRv, W, TLw] = gramSchmidt(jCol, V, W);
        else
            [V, TRv] = gramSchmidt(jCol, V);
        end
        Rv=Rv*TRv;
        Sv=TRv\Sv*TRv;
        if hermite
            Lw=Lw*TLw;
            Sw=TLw\Sw*TLw;
        end
    end
end

%orthogonalize columns from imaginary components
if Opts.orth
    for jCol=length(s0)+1:q
        if hermite
            [V, TRv, W, TLw] = gramSchmidt(jCol, V, W);
        else
            [V, TRv] = gramSchmidt(jCol, V);
        end
        Rv=Rv*TRv;
        Sv=TRv\Sv*TRv;
        if hermite
            Lw=Lw*TLw;
            Sw=TLw\Sw*TLw;
        end
    end
end

% reorthogonalization  
% Even modified Gram-Schmidt is not able to yield an orthonormal basis
% if the dimensions are high. Therefore, a reorthogonalization might be
% needed. On can choose to run modified GS again. From a theoretical 
% standpoint, this does not change the basis. However,
% numerically it is necessary to keep the numerics well behaved if the 
% reduced order is large
% The QR algorithm is much faster, however it does change the basis

switch Opts.reorth
    case 'mgs' %reorthogonalized GS
        Opts.orth='mgs'; %overwrite
        for jCol = 2:q        
            if hermite
                [V, TRv, W, TLw] = gramSchmidt(jCol, V, W);
            else
                [V, TRv] = gramSchmidt(jCol, V);
            end
            Rv=Rv*TRv;
            Sv=TRv\Sv*TRv;
            if hermite
                Lw=Lw*TLw;
                Sw=TLw\Sw*TLw;
            end
        end
    case 'qr'
       [V,Rq] = qr(V); %A=QR -> Q=A*inv(R) with transformation matrix inv(R)
       V=V(:,1:q);
       Rq=Rq(1:q,1:q);
       Rinv=Rq\eye(q);
       Rv=Rv*Rinv;
       Sv=Rq*Sv*Rinv;
       if hermite
           [W,Rq] = qr(W);
           W=W(:,1:q);
           Rq=Rq(1:q,1:q);
           Lw=Lw*Rinv;
           Sw=Rq*Sw*Rinv;
       end  
    case 0
    otherwise
        error('The orthogonalization chosen is incorrect or not implemented')
end

% output
if hermite
    varargout{1}=V;
    varargout{2}=W;
    varargout{3}=Sv;
    varargout{4}=Rv;
    varargout{5}=Sw;
    varargout{6}=Lw;
else
    varargout{1}=V;
    varargout{2}=Sv;
    varargout{3}=Rv;
end
    
   
    
    

% b) SECONDARY
function [V, TRv, W, TLw] = gramSchmidt(jCol, V, W)
%   Gram-Schmidt orthonormalization
%   Input:  jCol:  Column to be treated
%           V, W:  Krylov-Subspaces
%   Output: V, W:  orthonormal basis of Krylov-Subspaces
%           TRv, TLw: Transformation matrices

TRv=eye(size(V,2));
TLw=eye(size(V,2));
if jCol>1
    switch Opts.orth
        case 'dgks'
            % iterates standard gram-schmidt
            orthError=1;
            count=0;
            while(orthError>Opts.dgksTol)
                h=IP(V(:,1:jCol-1),V(:,jCol));
                V(:,jCol)=V(:,jCol)-V(:,1:jCol-1)*h;
                TRv(:,jCol)=TRv(:,jCol)-TRv(:,1:jCol-1)*h;
                if hermite
                    h=IP(W(:,1:jCol-1),W(:,jCol));
                    W(:,jCol)=W(:,jCol)-W(:,1:jCol-1)*h;
                    TLw(:,jCol)=TLw(:,jCol)-TLw(:,1:jCol-1)*h;
                end
                orthError=norm(IP([V(:,1:jCol-1),V(:,jCol)/sqrt(IP(V(:,jCol),V(:,jCol)))],...
                    [V(:,1:jCol-1),V(:,jCol)/sqrt(IP(V(:,jCol),V(:,jCol)))])-speye(jCol),'fro');
                if count>50 % if dgksTol is too small, Matlab can get caught in the while-loop
                    error('Orthogonalization of the Krylov basis failed due to the given accuracy.');
                end
                count=count+1;
            end
        case 'mgs'
            for iCol=1:jCol-1
              h=IP(V(:,jCol),V(:,iCol));
              V(:,jCol)=V(:,jCol)-V(:,iCol)*h;
              TRv(:,jCol)=TRv(:,jCol)-h*TRv(:,iCol);
              if hermite
                h=IP(W(:,jCol),W(:,iCol));
                W(:,jCol)=W(:,jCol)-W(:,iCol)*h;
                TLw(:,jCol)=TLw(:,jCol)-h*TLw(:,iCol);
              end 
            end
       case '2mgs'
            for k=0:1
                for iCol=1:jCol-1
                  h=IP(V(:,jCol),V(:,iCol));
                  V(:,jCol)=V(:,jCol)-V(:,iCol)*h;
                  TRv(:,jCol)=TRv(:,jCol)-h*TRv(:,iCol);
                  if hermite
                    h=IP(W(:,jCol),W(:,iCol));
                    W(:,jCol)=W(:,jCol)-W(:,iCol)*h;
                    TLw(:,jCol)=TLw(:,jCol)-h*TLw(:,iCol);
                  end 
                end
            end
        otherwise
            error('Opts.orth is invalid.');
    end  
end

% normalize new basis vector
h = sqrt(IP(V(:,jCol),V(:,jCol)));
V(:,jCol)=V(:,jCol)/h;
TRv(:,jCol) = TRv(:,jCol)/h;
if hermite
    h = sqrt(IP(W(:,jCol),W(:,jCol)));
    W(:,jCol)=W(:,jCol)/h;
    TLw(:,jCol) = TLw(:,jCol)/h;
end
end

function [V, Sv, Rv, W, Sw, Lw] = realSubspace(jCol, q, s0, V, Sv, Rv, W, Sw, Lw)
%   Split Krylov direction into real and imaginary to create a real 
%   Krylov subspace
%   Input:  jCol:  Column to be treated
%           q:     Reduction order
%           s0:    Vector containing the expansion points
%           V, W:  Krylov-Subspaces
%           Sv, Rsylv, Sw, Lsylv: Sylvester matrices
%   Output: V, W:  real basis of Krylov-Subspaces
%           Sv, Rv, Sw, Lw: real Sylvester matrices
nS0c=q-length(s0);
if ~isreal(s0(jCol))
    V(:,jCol+nS0c)=imag(V(:,jCol)); 
    V(:,jCol)=real(V(:,jCol));
    Rv(:,jCol+nS0c) = imag(Rv(:,jCol));
    Rv(:,jCol) = real(Rv(:,jCol));
    Sv(jCol, jCol+nS0c)=imag(Sv(jCol, jCol));
    Sv(jCol+nS0c, jCol)=-imag(Sv(jCol, jCol));
    Sv(jCol+nS0c, jCol+nS0c)=real(Sv(jCol, jCol));
    Sv(jCol, jCol)=real(Sv(jCol,jCol));
    if hermite, 
        W(:,jCol+nS0c)=imag(W(:,jCol));
        W(:,jCol)=real(W(:,jCol)); 
        Lw(:,jCol+nS0c) = imag(Lw(:,jCol));
        Lw(:,jCol) = real(Lw(:,jCol));
        Sw(jCol, jCol+nS0c)=imag(Sw(jCol, jCol));
        Sw(jCol+nS0c, jCol)=-imag(Sw(jCol, jCol));
        Sw(jCol+nS0c, jCol+nS0c)=real(Sw(jCol, jCol));
        Sw(jCol, jCol)=real(Sw(jCol,jCol));
    end
end
end

function [V, SRsylv, Rsylv, W, SLsylv, Lsylv] = nextDirection(jCol, s0, V, W)  
%   Get the next direction by solving the lse
%   Input:  jCol:  Column to be treated
%           s0:    Vector containing the expansion points
%           V, W:  solution of lse/Krylov subspaces
%   Output: V, W:  Updated lse solutions/Krylov subspace
%           SRsylv: update of column jCol of the Sylvester matrices
%                  Sv (e.g. SRsylv(:,jCol)=SRsylv)
%           Rsylv: update of column jCol of the Sylvester matrices 
%                  Rv (Rsylv either eye(size(B,2)) or 
%                  zeros(size(B,2)), e.g. Rsylv(:,jCol)=Rsylv*Rt(:,jCol)
%           SLsylv: update of column jCol of the Sylvester matrices
%                  Sw (e.g. SLsylv(:,jCol)=SLsylv)
%           Lsylv: update of column jCol of the Sylvester matrices 
%                  Lw (Lsylv either eye(size(C,1)) or 
%                  zeros(size(C,1)), e.g. Lsylv(:,jCol)=Lsylv*Lt(:,jCol)

SRsylv=zeros(size(V,2),1);
if hermite
    SLsylv=zeros(size(W,2),1);
end

% build Krylov subspace or just solve lse with current s0, Rt and B
switch Opts.krylov
    case 0
        tempV=B*Rt(:,jCol);
        newlu=1;
        newtan=1;
        SRsylv(jCol)=s0(jCol);
        Rsylv=eye(size(B,2));
        if hermite
            SLsylv(jCol)=s0(jCol);
            Lsylv=eye(size(C,1));
            tempW = C.'*Lt(:,jCol);
        end
        if hermite
            [V(:,jCol), W(:,jCol)] = lse(newlu, newtan, jCol, s0, tempV, tempW);
        else
            V(:,jCol) = lse(newlu, newtan, jCol, s0, tempV);
        end
    case 'standardKrylov'
        % new basis vector
        tempV=B*Rt(:,jCol); newlu=1; newtan=1;
        SRsylv(jCol)=s0(jCol);
        Rsylv=eye(size(B,2));
        if hermite
            SLsylv(jCol)=s0(jCol);
            Lsylv=eye(size(C,1));
            tempW = C.'*Lt(:,jCol);
        end
        if jCol>1
            if s0(jCol)==s0(jCol-1)
                newlu=0;
                if Rt(:,jCol) == Rt(:,jCol-1)
                    % Higher order moments, for the SISO and MIMO case
                    newtan = 0;
                    tempV = V(:,jCol-1); %overwrite
                    SRsylv(jCol-1)=1;
                    Rsylv=zeros(size(B,2));
                    if hermite
                        SLsylv(jCol-1)=1;
                        Lsylv=zeros(size(C,1));
                        tempW = W(:,jCol-1); 
                    end
                else
                    newtan = 1;
                end
            end
        end
        if hermite
            [V(:,jCol), W(:,jCol)] = lse(newlu, newtan, jCol, s0, tempV, tempW);
        else
            V(:,jCol) = lse(newlu, newtan, jCol, s0, tempV);
        end
    case 'cascadedKrylov'
        if size(B,2)==1
            newlu=1; newtan=1;
            SRsylv(jCol)=s0(jCol);
            if hermite
                SLsylv(jCol)=s0(jCol);
            end
            if jCol==1
                tempV=B;
                Rsylv=1;
                if hermite 
                    tempW=C.';
                    Lsylv=1;
                end
            else
                if s0(jCol)==s0(jCol-1)
                    newlu=0;
                    tempV=V(:,jCol-1);
                    if hermite
                        tempW=W(:,jCol-1);
                    end
                else
                    tempV=E*V(:,jCol-1);
                    if hermite
                        tempW=E*W(:,jCol-1);
                    end
                end
                Rsylv=0;
                SRsylv(jCol-1)=1;
                if hermite
                    Lsylv=0;
                    SLsylv(jCol-1)=1;
                end
            end
            if hermite
                [V(:,jCol), W(:,jCol)] = lse(newlu, newtan, jCol, s0, tempV, tempW);
            else
                V(:,jCol) = lse(newlu, newtan, jCol, s0, tempV);
            end
        else
            error('A cascaded Krylov basis is only available for SISO systems.');
        end
    otherwise 
        error('Opts.krylov is invalid.');
end
end

function [tempV, tempW] = lse(newlu, newtan, jCol, s0, tempV, tempW)
%   Solve linear system of equations to obtain the new direction
%   Moment matching:  newlu=0: tempV=(A-s0_old*E)^-1*(E*tempV)
%                     newlu=1: tempV=(A-s0*E)^-1*tempV
%   Markov parameter: newlu=0: tempV=E^-1*(A*tempV)
%                     newlu=1: tempW=E^-1*tempV
%   Input:  newlu: new lu decomposition required
%           newtan: new tangential direction required
%           jCol: Column to be treated
%           s0: Vector containing the expansion points
%           tempV, tempW: previous direction
%   Output: tempV, tempW: new direction
persistent R S L U a o;

if isinf(s0(jCol)) %Realization problem (match Markov parameters)
    if newlu==0 || strcmp(Opts.krylov,'cascade')
        tempV=A*tempV;
        if hermite
            tempW=A*tempW;
        end
    end
    if newlu==1
        try
            % compute Cholesky factors of E
            U=[];
            [R,p,S] = chol(E);
            if p~=0 % different factorization?
                error('solveLse:cholp','chol: p~=0');
            end
        catch err
            if strcmp(err.identifier,'MATLAB:posdef') || strcmp(err.identifier,'solveLse:cholp')|| ~strcmp(Opts.lse,'sparse')
                % E is not pos. def -> use LU instead
                switch Opts.lse
                    case 'sparse'
                        [L,U,a,o,S]=lu(E,'vector');
                    case 'full'
                        [L,U]=lu(E);
                end
            else
                rethrow(err);
            end
        end
    end
    if ~isempty(U) || strcmp(Opts.lse,'hess')
        switch Opts.lse
            case 'sparse'
                tempV(o,:) = U\(L\(S(:,a)\tempV)); %LU x(o,:) = S(:,a)\b 
                if hermite
                    tempW(o,:) = U\(L\(S(:,a)\tempW));
                end
            case 'full'
                tempV = U\(L\tempV);
                if hermite
                    tempW = U\(L\tempW);
                end
            case 'hess'
                tempV = E\tempV;
                if hermite
                    tempW = E\tempW;
                end
        end
    else
        tempV = S*(R\(R.'\(S.'*tempV)));
        if hermite
            tempW = S*(R\(R.'\(S.'*tempW)));
        end
    end
else %Rational Krylov
    if ~strcmp(Opts.lse,'iterative') %direct methods
        if newlu==0
            if size(B,2)==1 %SISO
                tempV=E*tempV;
                if hermite, tempW = E.'*tempW; end
            elseif newtan==0
                % Tangential matching of higher order moments
                tempV=E*tempV;
                if hermite, tempW = E.'*tempW; end
            end
        end
        if newlu==1
            if withoutE
                switch Opts.lse
                    case 'sparse'
                        % vector LU for sparse matrices
                        [L,U,a,o,S]=lu(A,'vector');
                    case 'full'
                        [L,U] = lu(A);
                end
            else
                switch Opts.lse
                case 'sparse'
                    % vector LU for sparse matrices
                    [L,U,a,o,S]=lu(A-s0(jCol)*E,'vector');
                case 'full'
                    [L,U] = lu(A-s0(jCol)*E);
                end
            end
        end
        % Solve the linear system of equations
        switch Opts.lse
            case 'sparse'
                tempV(o,:) = U\(L\(S(:,a)\tempV)); %LU x(o,:) = S(:,a)\b 
                if hermite, tempW = (S(:,a)).'\(L.'\(U.'\(tempW(o,:)))); end %U'L'S(:,a) x = c'(o,:) 
            case 'full'
                tempV = U\(L\tempV);
                if hermite, tempW = (L.'\(U.'\(tempW))); end 
            case 'hess'
                tempV = (A-s0(jCol)*E)\tempV;
                if hermite, tempW = (A-s0(jCol)*E).'\tempW; end 
        end
    else %iterative methods
        if withoutE
            [tempV,flag,method] = iterSolve(A,tempV,newlu,hermite);
            if Opts.verbose, disp(method); end
            if hermite
                [tempW,flag,method] = iterSolve((A).',tempW,0,hermite);
                if Opts.verbose, disp('entering hermite section'); end
            end
        else                
            if newlu
                [tempV,flag,method] = iterSolve(A-s0(jCol)*E,tempV,newlu,hermite);
                if Opts.verbose, disp(method); end
                if hermite
                    [tempW,flag,method] = iterSolve((A-s0(jCol)*E).',tempW,0,hermite);
                    if Opts.verbose, disp('entering hermite section'); end
                end
            else
                [tempV,flag,method] = iterSolve(A-s0(jCol)*E,E*tempV,newlu,hermite);
                if Opts.verbose, disp(method); end
                if hermite
                    [tempW,flag,method] = iterSolve((A-s0(jCol)*E).',E.'*tempW,0,hermite);
                    if Opts.verbose, disp('entering hermite section'); end
                end
            end
        end
    end
end
end

function [ x, varargout ] = iterSolve( A, b, newlu, hermite)
%ITERSOLVE function to solve the linear system A*x=b iteratively
%   Detailed explanation goes here

persistent first L U sym pd flag method failed nolu solver;   

if isempty(first);
    first = 1;
end

if first == 1           % set persistent variables on first call of function

    failed = {};

    if hermite==1 && newlu==0;             % LU factorization for transposed system: A~L*U => A'~U'*L' 
        if ~isempty(L) && ~isempty(U)      % ilu, not ichol is already computed
            temp = L;
            L = U.';
            U = temp.';
            clear temp;               
        end

        if sym==1 && pd==1
            A = -A;
            b = -b;
        end            

    end


end

if newlu==1 && first==1            % new L and U required (different s0)

    nolu = 0;        
    L = [];
    U = [];

    if isempty(solver)
        solver = Opts.solver;
    end

    % check for symmetry and definiteness

    if norm(A-A','fro')/norm(A,'fro') < 1e-10
        sym = 1;
    else
        sym = 0;
    end

    if ispd(-A)
        A = -A;
        b = -b;
        pd = 1;
    else
        pd = 0;
    end 
end

first = 0;                   % first round over

solver = selectSolver(solver, failed);

%     analyse = ['maxiterlse: ',num2str(Opts.maxiterlse),', tol: ',num2str(Opts.tollse)];
%     disp(analyse);

if sym && pd      % pcg
    if nolu~=1
        if isempty(L) && isempty(U)
            try
                L = ichol(A);
                [x,flag] = pcg(A,b,Def.tollse,Opts.maxiterlse,L,L.',b);
                method = 'pcg with ichol';
            catch 
                try
                    [L,U] = ilu(A);
                    [x,flag] = pcg(A,b,Opts.tollse,Opts.maxiterlse,L,U,b);
                    method = 'pcg with ilu';
                catch
                    [x,flag] = pcg(A,b,Opts.tollse,Opts.maxiterlse,[],[],b);
                    method = 'cg (no preconditioner)';
                    nolu = 1;
                end
            end
        elseif isempty(U)   % ichol already available
            [x,flag] = pcg(A,b,Def.tollse,Opts.maxiterlse,L,L.',b);
            method = 'pcg with ichol';
        else                % ilu already available
            [x,flag] = pcg(A,b,Opts.tollse,Opts.maxiterlse,L,U,b);
            method = 'pcg with ilu';
        end
    else
        [x,flag] = pcg(A,b,Opts.tollse,Opts.maxiterlse,[],[],b);
        method = 'cg (no preconditioner)';
    end
else
    switch solver
        case 'cgs'
            if nolu ~= 1                          % LU factorization not yet failed
                if isempty(L) && isempty(U)        % LU factorization not yet computed
                    try
                        [L,U] = ilu(A);
                        method = 'cgs with ilu';
                    catch
                        method = 'cgs (no preconditioner)';
                        nolu = 1;   
                    end                          
                else
                    method = 'cgs with ilu';                       
                end
            else
                method = 'cgs (no preconditioner)';
            end

            [x,flag] = cgs(A,b,Opts.tollse,Opts.maxiterlse,L,U,b);

            if flag ~= 0
                failed{length(failed)+1} = 'cgs';

                if Opts.verbose == 1              % print failure
                    msg = ['cgs did not converge.'];
                    warning(msg);
                end
            end

        case 'bicgstab'
            if nolu ~= 1                          % LU factorization not yet failed
                if isempty(L) && isempty(U)        % LU factorization not yet computed
                    try
                        [L,U] = ilu(A);
                        method = 'bicgstab with ilu';
                    catch
                        method = 'bicgstab (no preconditioner)';
                        nolu = 1;   
                    end                          
                else
                    method = 'bicgstab with ilu';                       
                end
            else
                method = 'bicgstab (no preconditioner)';
            end

            [x,flag] = bicgstab(A,b,Opts.tollse,ceil(Opts.maxiterlse/2),L,U,b);

            if flag ~= 0
                failed{length(failed)+1} = 'bicgstab';
                if Opts.verbose == 1              % print failure
                    msg = ['bicgstab did not converge.'];
                    warning(msg);
                end                    
            end

        case 'bicg'
            if nolu ~= 1                           % LU factorization not yet failed
                if isempty(L) && isempty(U)        % LU factorization not yet computed
                    try
                        [L,U] = ilu(A);
                        method = 'bicg with ilu';
                    catch
                        method = 'bicg (no preconditioner)';
                        nolu = 1;   
                    end                          
                else
                    method = 'bicg with ilu';                       
                end
            else
                method = 'bicg (no preconditioner)';
            end

            [x,flag] = bicg(A,b,Opts.tollse,ceil(Opts.maxiterlse/2),L,U,b);

            if flag ~= 0
                failed{length(failed)+1} = 'bicg';
                if Opts.verbose == 1              % print failure
                    msg = ['bicg did not converge.'];
                    warning(msg);
                end                    

            end
    end

    if flag ~= 0 && length(failed) < 3
        iterSolve(A,b,newlu,hermite);
    end

end

first = [];

if(nargout == 2)
    varargout{1} = flag;
end 

if(nargout == 3)
    varargout{1} = flag;
    varargout{2} = method;
end

if flag ~= 0 
    msg = ['Could not achieve desired tolerance (',num2str(Opts.tollse),') within ',num2str(Opts.maxiterlse),' iterations. Results may be inaccurate!'];
    if Opts.force == 1;
        warning(msg);
    else
        error(msg);
    end
    solver = Opts.solver;
end

end

function [solver] = selectSolver( solver, failed )
   solvers = {'cgs','bicgstab','bicg'}; 
   switch length(failed)
       case 0
       case 1
           index = strcmp(failed, solvers);
           solver = solvers{find(index==0,1)};

       case 2
           index1 = strcmp(failed{1}, solvers);
           index2 = strcmp(failed{2}, solvers);
           solver = solvers{index1==index2};
   end

end

end