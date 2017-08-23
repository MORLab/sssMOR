function [V, Sv, Rv, W, Sw, Lw, nLU] = arnoldi(E,A,B,varargin)
% ARNOLDI - Arnoldi algorithm for Krylov subspaces with multiple shifts
% 
% Syntax:
%       V                                = ARNOLDI(E,A,B,s0)
%       [V,Sv,Rv]                        = ARNOLDI(E,A,B,s0)
%       [V,Sv,Rv]                        = ARNOLDI(E,A,B,s0,IP)
%       [V,Sv,Rv]                        = ARNOLDI(E,A,B,s0,Rt)
%       [V,Sv,Rv]                        = ARNOLDI(E,A,B,s0,Rt,IP)
%       [V,Sv,Rv,W,Sw,Lw]                = ARNOLDI(E,A,B,C,s0)
%       [V,Sv,Rv,W,Sw,Lw]                = ARNOLDI(E,A,B,C,s0,IP)
%       [V,Sv,Rv,W,Sw,Lw]                = ARNOLDI(E,A,B,C,s0,Rt,Lt)
%       [V,Sv,Rv,W,Sw,Lw]                = ARNOLDI(E,A,B,C,s0,Rt,Lt,IP)
%       [V,Sv,Rv,W,Sw,Lw,nLU]         	 = ARNOLDI(E,A,B,...)
%       [V,...]                          = ARNOLDI(E,A,B,...,Opts)
% 
% Description:
%       This function is used to compute the matrix V spanning the 
%       rational input Krylov subspace corresponding to E, A, b and s0 [1-3].
%
%       The input Krylov subpspace of order q correpsonding to a single 
%       complex expansion point s_0 is defined as
%
%       $$ Im(V) = span\left\{ (A-s_0E)^{-1} b_t,\; \dots,\, \left[(A-s_0E)^{-1}E\right]^{q-1}(A-s_0E)^{-1}b_t\right\}. $$
%
%       In this case, $$ b_t $$ is either:
%
%       * the input vector of a SISO model,
%       * the input matrix of a MIMO model (block Krylov)
%       * the input matrix multiplied by a tangential direction (tangential Krylov)
%
%       s0 must be a vector of complex frequencies closed under conjugation. 
%       In case of MIMO models, if matrices of tangential directions Rt 
%       (and Lt) are defined, they must have the same number of columns as 
%       the shifts, so that for each tangential direction it is clear to 
%       which shift it belongs. If no tangential directions are specified,
%       then block Krylov subspaces are computed.
%
%       //Note: For MIMO models, block Krylov subpspaces 
%       with multiplicities in the shifts are not supported so far.
%
%       If in addition, the output matrix C is passed, then ARNOLDI
%       computes input and output Krylov subspaces corresponding to the
%       same expansion points. The resulting matrices V, W can be used for
%       Hermite interpolation.
%
%       In this case, the output Krylov subspace is defined as
%
%       $$ Im(W) = span\left\{ (A-s_0E)^{-T} c_t^T,\; \dots,\, \left[(A-s_0E)^{-T}E^T\right]^{q-1}(A-s_0E)^{-T}c_t^T\right\}. $$
%
%       The columns of V build an orthonormal basis of the input Krylov 
%       subspace. The orthogonalization is conducted using a 
%       reorthogonalized modified Gram-Schmidt procedure [4] or dgks
%       orthogonalization [5] with respect to the inner product defined in 
%       IP (optional). If no inner product is specified, then the euclidian
%       product corresponding to I is chosen by default:
%
%                       IP=@(x,y) (x.'*y)
%
%       If specified, this function computes the Sylvester matrices
%       corresponding to the Krylov subspaces. The matrices Sv and 
%       Rsylv satisfy the input Sylvester equation given as
%
%       $$ A V - E V S_v - B R_v = 0 \quad          (1)$$
%
%       and the output Sylvester matrices Sw and Lsylv are
%       accordingly defined by
%
%       $$ A^T W - E^T W S_w^T - C^T L_w = 0 \quad          (2)$$
%
%       Note that this function does not solve the Sylvester equations, 
%       but constructs the Sylvester matrices together with the Krylov 
%       subspaces.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -E/A/B/C:  System matrices
%       -s0:       Vector of complex conjuate expansion points
%
%       *Optional Input Arguments:*
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
%                           [{'sparse'} / 'full' / 'hess' /'iterative']
%           -.dgksTol:      tolerance for dgks orthogonalization
%                           [{1e-12} / positive float]
%           -.krylov:       standard or cascaded krylov basis
%                           [{'standardKrylov'} / 'cascadedKrylov']
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
%       -V:        Orthonormal basis spanning the input Krylov subsp. 
%       -Sv:       Matrix of input Sylvester Eq. (1)
%       -Rv:       Right tangential directions of Sylvester Eq. (1), (mxq) matrix
%       -W:        Orthonormal basis spanning the output Krylov subsp.
%       -Sw:       Matrix of output Sylvester Eq. (2)
%       -Lw:       Left tangential directions of Sylvester Eq. (2), (pxq) matrix
%       -nLU:      Number of LU decompositions required
%
% See Also: 
%       rk, irka, projectiveMor
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
% Last Change:  09 Aug 2017
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Define execution parameters
if ~isempty(varargin) && isstruct(varargin{end});
    %Options defined
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end

Def.real = true; %keep the projection matrices real?
Def.orth = '2mgs'; %orthogonalization after every direction {0,'dgks','mgs','2mgs'}
Def.reorth = 0; %reorthogonaliation at the end {0, 'mgs', 'qr'}
Def.lse = 'sparse'; %use sparse or full LU or lse with Hessenberg decomposition {'sparse', 'full','hess'}
Def.dgksTol = 1e-12; %orthogonality tolerance: norm(V'*V-I,'fro')<tol
Def.krylov = 'standardKrylov'; %standard or cascaded krylov basis (only for siso) {'standardKrylov','cascadedKrylov'}
        
% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end              
 
if islogical(Opts.orth) && Opts.orth
    Opts.orth = Def.orth; %set default orthogonalization method
end
%%  Parse input
if length(varargin) == 1
    % usage: ARNOLDI(E,A,B,s0)
    s0 = varargin{1};
    hermite = 0; % same shifts for input and output Krylov?
elseif length(varargin) > 1
    %   Do the classification depending on the properties of the objects
    %   ARNOLDI(E,A,B,s0,...) or ARNOLDI(E,A,B,C,...)
    if size(varargin{1},2) == size(A,1)
        % usage: ARNOLDI(E,A,B,C,s0,...)
        hermite = 1;
        C = varargin{1};
        s0 = varargin{2};
        if length(varargin) == 3
            % usage: ARNOLDI(E,A,B,C,s0,IP)
            IP = varargin{3};
        elseif length(varargin) == 4
            % usage: ARNOLDI(E,A,B,C,s0,Rt,Lt)
            Rt = varargin{3};
            Lt = varargin{4};
        elseif length(varargin) == 5
            % usage: ARNOLDI(E,A,B,C,s0,Rt,Lt,IP)
            Rt = varargin{3};
            Lt = varargin{4};
            IP = varargin{5};
        end
    else
        % usage: ARNOLDI(E,A,B,s0,...)
        hermite = 0;
        s0 = varargin{1};
        if length(varargin) == 2
            if size(varargin{2},2) == size(s0,2)
                % usage: ARNOLDI(E,A,B,s0,Rt)
                Rt = varargin{2};
            else   
                % usage: ARNOLDI(E,A,B,s0,IP)
                IP = varargin{2};
            end
        else
            % usage: ARNOLDI(E,A,b,s0,Rt,IP)
            Rt = varargin{2};
            IP = varargin{3};
        end
    end
end

if size(s0,1)>1
    error('s0 must be a vector containing the expansion points.')
end

if exist('Rt','var') && ~isempty(Rt)
    if length(s0) ~= size(Rt,2),
        error('Rt must have the same columns as s0')
    end
    %   The reduced order is equivalent to the number of shifts
    q = length(s0);
else
    %   Block Krylov subspaces will be performed
    q = length(s0)*size(B,2);
end

if exist('Lt','var') && ~isempty(Lt)
    if length(s0) ~= size(Lt,2),
        error('Lt must have the same columns as s0')
    end
end

% IP
if ~exist('IP', 'var') 
   IP=@(x,y) (x'*y); %seems to be better conditioned that E norm
end

%% ---------------------------- CODE -------------------------------
% Real reduced system
if Opts.real
    s0 = updateS0(s0);
    nLU = length(unique(s0)); %get number of LU decompositions required
else
    nLU = length(unique(updateS0(s0)));
end

% Tangential directions
if ~exist('Rt', 'var') || isempty(Rt)%   Compute block Krylov subspaces
    s0 = tangentialDirection(s0);
end

% preallocate memory
V   = zeros(length(B),q);
Rv  = zeros(size(B,2),q);
Sv  = zeros(q);
if hermite 
    W   = zeros(length(B),q); 
    Lw  = zeros(size(C,1),q);
    Sw  = zeros(q);
else % initialize outputs
    W   = [];
    Lw  = [];
    Sw  = [];
end

% Compute hess
if strcmp(Opts.lse,'hess')
    [A,E,Q,Z] = hess(full(A),full(E)); B = Q*B; 
    if hermite
        C = C*Z; 
    end
    Opts.lse='gauss';
end

% Compute the Krylov subspaces
for jCol=1:length(s0)
    if hermite
        [V, W, SRsylv, Rsylv, SLsylv, Lsylv]    = solveLse(jCol, V, W, A, B, C, E, s0, Rt, Lt, Opts);
    else
        [V, SRsylv, Rsylv]                      = solveLse(jCol, V, A, B, E, s0, Rt, Opts);
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
            [V, Sv, Rv, W, Sw, Lw]  = realSubspace(jCol, q, s0, V, Sv, Rv, W, Sw, Lw);
        else
            [V, Sv, Rv]             = realSubspace(jCol, q, s0, V, Sv, Rv);
        end
    end

    if Opts.orth
        if hermite
            [V, TRv, W, TLw]    = gramSchmidt(jCol, V, W);
        else
            [V, TRv]            = gramSchmidt(jCol, V);
        end
        Rv  = Rv*TRv;
        Sv  = TRv\Sv*TRv;
        if hermite
            Lw  = Lw*TLw;
            Sw  = TLw\Sw*TLw;
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

% transform hess back
if exist('Z','var') && ~isempty(Z)
    V=Z*V;
    if hermite
        W=Q.'*W;
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

function [s0] = tangentialDirection(s0)
%   Update s0 and define Rt for calculation of tangential directions
%   Input:  s0: Vector containing the expansion points
%   Output: s0: Vector containing the updated expansion points
if size(B,2) == 1; %SISO -> tangential directions are scalars
    Rt = ones(1,length(s0));
    % these siso "tangential directions" are not used for the
    % computatoin of the Krylov subspaces but just for the computation
    % of the transformed tangential directions 
else %MIMO -> fill up s0 and define tangential blocks

    % tangential matching of higher order moments not implemented so
    % far! Therefore, if two shifts are the same, give an error

    if any(diff(sort(s0))==0)
        s0=sort(s0);
        s0(end+1)=NaN;
        [us0,ia,~]=unique(s0);
        ns0=diff(ia);
        us0=us0(1:end-1);

        Rt=[];
        s0=[];
        for i=1:length(us0)
            tempRt=[];
            for j=1:size(B,2)
               tempRt=blkdiag(tempRt,ones(1,ns0(i)));
               s0=[s0,us0(i)*ones(1,ns0(i))];
            end
            Rt=[Rt,tempRt];
        end
    else
        s0old = s0; nS0=length(s0); s0 = [];
        for iShift = 1:nS0
            s0 = [s0, s0old(iShift)*ones(1,size(B,2))];
        end
        Rt = repmat(speye(size(B,2)),1,nS0);
    end
end
if hermite
    if size(B,2) ~=size(C,1)
        error('Block Krylov for m~=p is not supported in arnoldi');
    else
        Lt = Rt;
    end
end
end

function [s0] = updateS0(s0)
%   Remove one element of complex expansion points
%   Input:  s0: Vector containing the expansion points
%   Output: s0: Sorted vector containing only real or one element of
%               complex conjugated expansion points
%           nS0c: Number of complex conjugated expansion points

% remove one element of complex pairs (must be closed under conjugation)
k=find(imag(s0));
if ~isempty(k)
    % make sure shift are sorted and come in complex conjugate pairs
    try 
        s0cUnsrt = s0(k);
        s0c = cplxpair(s0cUnsrt);
        % get permutation indices, since cplxpair does not do it for you
        [~,cplxSorting] = ismember(s0c,s0cUnsrt); %B(idx) = A
    catch 
        error(['Shifts must come in complex conjugated pairs and be sorted',...
            ' before being passed to arnoldi.'])
    end

    % take only one shift per complex conjugate pair
    s0(k) = []; 
    s0 = [s0 s0c(1:2:end)];

    % take only one residue vector for each complex conjugate pair
    if exist('Rt','var') && ~isempty(Rt)
        RtcUnsrt = Rt(:,k); 
        Rtc = RtcUnsrt(:,cplxSorting);
        Rt(:,k) = []; 
        Rt = [Rt,Rtc(:,1:2:end)]; 
        if exist('Lt','var') && ~isempty(Lt)
            LtcUnsrt = Lt(:,k);
            Ltc = LtcUnsrt(:,cplxSorting);
            Lt(:,k) = [];
            Lt = [Lt,Ltc(:,1:2:end)];
        end
    end
end
end

%b) SECONDARY
function [V, TRv, W, TLw]   = gramSchmidt(jCol, V, W)
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

end
