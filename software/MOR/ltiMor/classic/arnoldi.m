function [V, SRsylv, CRsylv, W, SLsylv, CLsylv] = arnoldi(E,A,B,varargin)
% ARNOLDI - Arnoldi algorithm for Krylov subspaces with multiple shifts
% 
% Syntax:
%       V                                = ARNOLDI(E,A,B,s0)
%       [V,SRsylv,Rsylv]                 = ARNOLDI(E,A,B,s0)
%       [V,SRsylv,Rsylv]                 = ARNOLDI(E,A,B,s0,IP)
%       [V,SRsylv,Rsylv]                 = ARNOLDI(E,A,B,s0,Rt)
%       [V,SRsylv,Rsylv]                 = ARNOLDI(E,A,B,s0,Rt,IP)
%       [V,SRsylv,Rsylv,W,SLsyslv,Lsylv] = ARNOLDI(E,A,B,C,s0)
%       [V,SRsylv,Rsylv,W,SLsyslv,Lsylv] = ARNOLDI(E,A,B,C,s0,IP)
%       [V,SRsylv,Rsylv,W,SLsyslv,Lsylv] = ARNOLDI(E,A,B,C,s0,Rt,Lt)
%       [V,SRsylv,Rsylv,W,SLsyslv,Lsylv] = ARNOLDI(E,A,B,C,s0,Rt,Lt,IP)
%       [V,...]	= ARNOLDI(E,A,B,C,s0,...,Opts)
% 
% Description:
%       This function is used to compute the matrix V spanning the 
%       rational input Krylov subspace corresponding to E, A, b and s0 [1-3].
%
%       s0 must be a vector of complex frequencies closed under conjugation. 
%       In case of MIMO systems, if matrices of tangential directions Rt 
%       (and Lt) are defined, they must have the same number of columns as 
%       the shifts, so that for each tangential direction it is clear to 
%       which shift it belongs. If not tangential directions are specified,
%       then block Krylov subspaces are computed.
%
%       If in addition, the output matrix C is passed, then ARNOLDI
%       computes input and output Krylov subspaces corresponding to the
%       same expansion points. The resulting matrices V, W can be used for
%       Hermite interpolation.
%
%       //Note: for MIMO systems, block Krylov subpspaces 
%       with multiplicities in the shifts are not supported so far.
%
%       The columns of V build an orthonormal basis of the input Krylov 
%       subspace. The orthogonalization is conducted using a 
%       reorthogonalized modified Gram-Schmidt procedure [4] or dgks
%       orthogonalization [5] with respect to the inner product defined in 
%       IP (optional). If no inner product is specified, then the euclidian
%       product corresponding to I is chosen by default:
%
%                       IP=@(x,y) (x'*y)
%
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
%                           [{'real'} / '0']
%           -.orth:         orthogonalization of new projection direction
%                           [{'2mgs'} / 0 / 'dgks' / 'mgs']
%           -.reorth:       reorthogonalization
%                           [{'gs'} / 0 / 'qr']
%           -.lse:          use LU or hessenberg decomposition
%                           [{'sparse'} / 'full' / 'hess']
%           -.dgksTol:      tolerance for dgks orthogonalization
%                           [{1e-12} / positive float]
%           -.krylov:       standard or cascaded krylov basis
%                           [{0} / 'cascade]
%
% Output Arguments:
%       -V:        Orthonormal basis spanning the input Krylov subsp. 
%       -SRsylv:   Matrix of input Sylvester Eq.
%       -Rsylv:    Right tangential directions of Sylvester Eq.
%       -W:        Orthonormal basis spanning the output Krylov subsp.
%       -SLsylv:   Matrix of output Sylvester Eq.
%       -Lsylv:    Left tangential directions of Sylvester Eq.
%
% See Also: 
%       rk, irka
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
% Authors:      Heiko Panzer, Alessandro Castagnotto 
%               Lisa Jeschek
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  13 Jan 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Define execution parameters
if ~isempty(varargin) && isstruct(varargin{end});
    %Options defined
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end

Def.real = 'real'; %keep the projection matrices real?
Def.orth = '2mgs'; %orthogonalization after every direction {0,'dgks','mgs','2mgs'}
Def.reorth = 0; %reorthogonaliation at the end {0, 'mgs', 'qr'}
Def.lse = 'sparse'; %use sparse or full LU or lse with Hessenberg decomposition {'sparse', 'full','hess'}
Def.dgksTol = 1e-12; %orthogonality tolerance: norm(V'*V-I,'fro')<tol
Def.krylov = 0; %standard or cascaded krylov basis (only for siso) {0,'cascade'}
        
% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
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

% If the 'full' option is selected for LU, convert E,A once to full
if strcmp(Opts.lse,'full')
    E = full(E); A = full(A);
elseif strcmp(Opts.lse,'hess')
    [A,E,Q,Z] = hess(full(A),full(E)); B = Q*B; if hermite, C = C*Z; end
else
    E = sparse(E); A=sparse(A);
end

%% ---------------------------- CODE -------------------------------
% Real reduced system
if strcmp(Opts.real,'real')
    s0 = updateS0(s0);
end

% Tangential directions
if ~exist('Rt', 'var') || isempty(Rt)%   Compute block Krylov subspaces
    s0 = tangentialDirection(s0);
end

% Compute the Krylov subspaces
if hermite
    [V, SRsylv, CRsylv, W, SLsylv, CLsylv] = krylovSubspace(s0, q);
else
    [V, SRsylv, CRsylv] = krylovSubspace(s0, q);
end
    

%% ------------------ AUXILIARY FUNCTIONS --------------------------
% a) PRIMARY
    function [V, SRsylv, CRsylv, W, SLsylv, CLsylv] = krylovSubspace(s0, q)
    %   Calculate Krylov Subspace of s0
    %   Input:  s0:  Vector containing the expansion points
    %           q:   Original length of s0 with complex conjugated elements
    %   Output: V, W:  Krylov-Subspace of s0
    %           Rsylv, Lsylv: Sylvester matrix
        
    % preallocate memory
    V=zeros(length(B),q);
    CRsylv=zeros(size(B,2),q);
    SRsylv=zeros(q);
    if hermite 
        W = zeros(length(B),q); 
        CLsylv = zeros(size(C,1),q);
        SLsylv=zeros(q);
    end
    for jCol=1:length(s0)
        if hermite
            [V, Ssylv, Csylv, W] = krylovDirection(jCol, s0, V, W);
        else
            [V, Ssylv, Csylv] = krylovDirection(jCol, s0, V);
        end
        SRsylv(:,jCol) = Ssylv;
        CRsylv(:,jCol) = Csylv*Rt(:,jCol);
        if hermite
            SLsylv(:,jCol) = Ssylv;
            CLsylv(:,jCol) = Csylv*Lt(:,jCol);
        end

        % split complex conjugate columns into real (->j) and imag (->j+length(s0c)/2
        if strcmp(Opts.real,'real')
            if hermite
                [V, SRsylv, CRsylv,W, SLsylv, CLsylv] = realSubspace(jCol, q, s0, V, SRsylv, CRsylv, W, SLsylv, CLsylv);
            else
                [V, SRsylv, CRsylv] = realSubspace(jCol, q, s0, V, SRsylv, CRsylv);
            end
        end

        if Opts.orth
            if hermite
                [V, TRsylv, W, TLsylv] = gramSchmidt(jCol, V, W);
            else
                [V, TRsylv] = gramSchmidt(jCol, V);
            end
            CRsylv=CRsylv*TRsylv;
            SRsylv=TRsylv\SRsylv*TRsylv;
            if hermite
                CLsylv=CLsylv*TLsylv;
                SLsylv=TLsylv\SLsylv*TLsylv;
            end
        end
    end

    %orthogonalize columns from imaginary components
    if Opts.orth
        for jCol=length(s0)+1:q
            if hermite
                [V, TRsylv, W, TLsylv] = gramSchmidt(jCol, V, W);
            else
                [V, TRsylv] = gramSchmidt(jCol, V);
            end
            CRsylv=CRsylv*TRsylv;
            SRsylv=TRsylv\SRsylv*TRsylv;
            if hermite
                CLsylv=CLsylv*TLsylv;
                SLsylv=TLsylv\SLsylv*TLsylv;
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
                    [V, TRsylv, W, TLsylv] = gramSchmidt(jCol, V, W);
                else
                    [V, TRsylv] = gramSchmidt(jCol, V);
                end
                CRsylv=CRsylv*TRsylv;
                SRsylv=TRsylv\SRsylv*TRsylv;
                if hermite
                    CLsylv=CLsylv*TLsylv;
                    SLsylv=TLsylv\SLsylv*TLsylv;
                end
            end
        case 'qr'
           [V,R] = qr(V); %A=QR -> Q=A*inv(R) with transformation matrix inv(R)
           V=V(:,1:q);
           R=R(1:q,1:q);
           Rinv=R\eye(q);
           CRsylv=CRsylv*Rinv;
           SRsylv=R*SRsylv*Rinv;
           if hermite
               [W,R] = qr(W);
               W=W(:,1:q);
               R=R(1:q,1:q);
               CLsylv=CLsylv*Rinv;
               SLsylv=R*SLsylv*Rinv;
           end  
        case 0
        otherwise
            error('The orthogonalization chosen is incorrect or not implemented')
    end 
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
            error(['Multiplicities in the shifts detected. Tangential '...
                'matching of higher order moments with block'...
                'Krylov not implemented (yet)!']);
        end
        
        s0old = s0; nS0=length(s0); s0 = [];
        for iShift = 1:nS0
            s0 = [s0, s0old(iShift)*ones(1,size(B,2))];
        end
        Rt = repmat(speye(size(B,2)),1,nS0);
        
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

% b) SECONDARY
    function [V, TRsylv, W, TLsylv] = gramSchmidt(jCol, V, W)
    %   Gram-Schmidt orthonormalization
    %   Input:  jCol:  Column to be treated
    %           V, W:  Krylov-Subspaces
    %   Output: V, W:  orthonormal basis of Krylov-Subspaces
    %           TRsylv, TLsylv: Transformation matrices
    
    TRsylv=eye(size(V,2));
    TLsylv=eye(size(V,2));
    if jCol>1
        switch Opts.orth
            case 'dgks'
                % iterates standard gram-schmidt
                orthError=1;
                count=0;
                while(orthError>Opts.dgksTol)
                    h=IP(V(:,1:jCol-1),V(:,jCol));
                    V(:,jCol)=V(:,jCol)-V(:,1:jCol-1)*h;
                    TRsylv(:,jCol)=TRsylv(:,jCol)-TRsylv(:,1:jCol-1)*h;
                    if hermite
                        h=IP(W(:,1:jCol-1),W(:,jCol));
                        W(:,jCol)=W(:,jCol)-W(:,1:jCol-1)*h;
                        TLsylv(:,jCol)=TLsylv(:,jCol)-TLsylv(:,1:jCol-1)*h;
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
                  TRsylv(:,jCol)=TRsylv(:,jCol)-h*TRsylv(:,iCol);
                  if hermite
                    h=IP(W(:,jCol),W(:,iCol));
                    W(:,jCol)=W(:,jCol)-W(:,iCol)*h;
                    TLsylv(:,jCol)=TLsylv(:,jCol)-h*TLsylv(:,iCol);
                  end 
                end
           case '2mgs'
                for k=0:1
                    for iCol=1:jCol-1
                      h=IP(V(:,jCol),V(:,iCol));
                      V(:,jCol)=V(:,jCol)-V(:,iCol)*h;
                      TRsylv(:,jCol)=TRsylv(:,jCol)-h*TRsylv(:,iCol);
                      if hermite
                        h=IP(W(:,jCol),W(:,iCol));
                        W(:,jCol)=W(:,jCol)-W(:,iCol)*h;
                        TLsylv(:,jCol)=TLsylv(:,jCol)-h*TLsylv(:,iCol);
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
    TRsylv(:,jCol) = TRsylv(:,jCol)/h;
    if hermite
        h = sqrt(IP(W(:,jCol),W(:,jCol)));
        W(:,jCol)=W(:,jCol)/h;
        TLsylv(:,jCol) = TLsylv(:,jCol)/h;
    end
    end

    function [V, SRsylv, CRsylv, W, SLsylv, CLsylv] = realSubspace(jCol, q, s0, V, SRsylv, CRsylv, W, SLsylv, CLsylv)
    %   Split Krylov direction into real and imaginary to create a real 
    %   Krylov subspace
    %   Input:  jCol:  Column to be treated
    %           nS0c:  Number of complex conjugated expansion points
    %           s0:    Vector containing the expansion points
    %           V, W:  Krylov-Subspaces
    %           CRsylv, CLsylv: Sylvester matrices
    %   Output: V, W:  real basis of Krylov-Subspaces
    %           CRsylv, CLsylv: real Sylvester matrices
    nS0c=q-length(s0);
    if ~isreal(s0(jCol))
        V(:,jCol+nS0c)=imag(V(:,jCol)); 
        V(:,jCol)=real(V(:,jCol));
        CRsylv(:,jCol+nS0c) = imag(CRsylv(:,jCol));
        CRsylv(:,jCol) = real(CRsylv(:,jCol));
        SRsylv(jCol, jCol+nS0c)=imag(SRsylv(jCol, jCol));
        SRsylv(jCol+nS0c, jCol)=-imag(SRsylv(jCol, jCol));
        SRsylv(jCol+nS0c, jCol+nS0c)=real(SRsylv(jCol, jCol));
        SRsylv(jCol, jCol)=real(SRsylv(jCol,jCol));
        if hermite, 
            W(:,jCol+nS0c)=imag(W(:,jCol));
            W(:,jCol)=real(W(:,jCol)); 
            CLsylv(:,jCol+nS0c) = imag(CLsylv(:,jCol));
            CLsylv(:,jCol) = real(CLsylv(:,jCol));
            SLsylv(jCol, jCol+nS0c)=imag(SLsylv(jCol, jCol));
            SLsylv(jCol+nS0c, jCol)=-imag(SLsylv(jCol, jCol));
            SLsylv(jCol+nS0c, jCol+nS0c)=real(SLsylv(jCol, jCol));
            SLsylv(jCol, jCol)=real(SLsylv(jCol,jCol));
        end
    end
    end
    
    function [V, Ssylv, Csylv, W] = krylovDirection(jCol, s0, V, W)  
    %   Get new Krylov direction
    %   Input:  jCol:  Column to be treated
    %           s0:    Vector containing the expansion points
    %           V, W:  Krylov subspace
    %   Output: V, W:  Updated Krylov subspace
    %           Ssylv: update of column jCol of the Sylvester matrices
    %                  SRsylv and SLsylv (e.g. SRsylv(:,jCol)=Ssylv)
    %           Csylv: update of column jCol of the Sylvester matrices 
    %                  CRsylv and CLsylv (Csylv either eye(size(B,2)) or 
    %                  zeros(size(B,2)), e.g. CRsylv(:,jColt)=Csylv*Rt(:,jCol)
    
    Ssylv=zeros(size(V,2),1);
    switch Opts.krylov
        case 0
            % new basis vector
            tempV=B*Rt(:,jCol); newlu=1; newtan=1;
            Ssylv(jCol)=s0(jCol);
            Csylv=eye(size(B,2));
            if hermite
                tempW = C'*Lt(:,jCol); 
            end
            if jCol>1
                if s0(jCol)==s0(jCol-1)
                    newlu=0;
                    if Rt(:,jCol) == Rt(:,jCol-1)
                        % Higher order moments, for the SISO and MIMO case
                        newtan = 0;
                        tempV = V(:,jCol-1); %overwrite
                        Csylv=zeros(size(B,2));
                        Ssylv(jCol-1)=1;
                        if hermite
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
        case 'cascade'
            if size(B,2)==1
                newlu=1; newtan=1;
                Ssylv(jCol)=s0(jCol);
                if jCol==1
                    tempV=B;
                    Csylv=1;
                    if hermite
                        tempW=C';
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
                    Csylv=0;
                    Ssylv(jCol-1)=1;
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
    %   Solve linear system of equations to obtain the new Krylov direction
    %   Moment matching:  newlu=0: tempV=(A-s0_old*E)^-1*(E*tempV)
    %                     newlu=1: tempV=(A-s0*E)^-1*tempV
    %   Markov parameter: newlu=0: tempV=E^-1*(A*tempV)
    %                     newlu=1: tempW=E^-1*tempV
    %   Input:  newlu: new lu decomposition required
    %           newtan: new tangential direction required
    %           jCol: Column to be treated
    %           s0: Vector containing the expansion points
    %           tempV, tempW: previous Krylov direction
    %   Output: tempV, tempW: new Krylov direction
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
                [R,~,S] = chol(E);
%                 R = chol(sparse(E));
            catch err
                if (strcmp(err.identifier,'MATLAB:posdef')) || ~strcmp(Opts.lse,'sparse')
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
            tempV = S*(R\(R'\(S'*tempV)));
            if hermite
                tempW = S*(R\(R'\(S'*tempW)));
            end
        end
    else %Rational Krylov
        if newlu==0
            if size(B,2)==1 %SISO
                tempV=E*tempV;
                if hermite, tempW = E'*tempW; end
            elseif newtan==0
                % Tangential matching of higher order moments
                tempV=E*tempV;
                if hermite, tempW = E'*tempW; end
            end
        end
        if newlu==1
            switch Opts.lse
                case 'sparse'
                    % vector LU for sparse matrices
                    [L,U,a,o,S]=lu(A-s0(jCol)*E,'vector');
                case 'full'
                    [L,U] = lu(A-s0(jCol)*E);
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
    end
    end
end
