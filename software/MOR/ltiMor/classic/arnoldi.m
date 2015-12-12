function [V,Rsylv,W,Lsylv] = arnoldi(E,A,B,varargin)
% ARNOLDI - Arnoldi algorithm for Krylov subspaces with multiple shifts
% 
% Syntax:
%       V					= ARNOLDI(E,A,B,s0)
%       [V,Rsylv]			= ARNOLDI(E,A,B,s0)
%       [V,Rsylv]			= ARNOLDI(E,A,B,s0,IP)
%       [V,Rsylv]			= ARNOLDI(E,A,B,s0,Rt)
%       [V,Rsylv]			= ARNOLDI(E,A,B,s0,Rt,IP)
%       [V,Rsylv,W,Lsylv]	= ARNOLDI(E,A,B,C,s0)
%       [V,Rsylv,W,Lsylv]	= ARNOLDI(E,A,B,C,s0,IP)
%       [V,Rsylv,W,Lsylv]	= ARNOLDI(E,A,B,C,s0,Rt,Lt)
%       [V,Rsylv,W,Lsylv]	= ARNOLDI(E,A,B,C,s0,Rt,Lt,IP)
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
%       reorthogonalized modified Gram-Schmidt procedure [4] with respect 
%       to the inner product  defined in IP (optional). If no inner product 
%       is specified, then the euclidian product corresponding to I is 
%       chosen by default:
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
%       -Rt,Lt:    Matrix of right/left tangential directions
%       -IP:       function handle for inner product
%
% Output Arguments:
%       -V:        Orthonormal basis spanning the input Krylov subsp. 
%       -Rsylv:    Right tangential directions of Sylvester Eq.
%       -W:        Orthonormal basis spanning the output Krylov subsp.
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
% Last Change:  12 Dec 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Define execution parameters
if ~isempty(varargin) && isstruct(varargin{end});
    %Options defined
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end

Def.makeOrth = 1; %make orthogonal?
Def.makeReal = 1; %keep the projection matrices real?
Def.reorth = 'dgks'; %gram schmidt reorthogonalization {0,'dgks','gs','qr'}
Def.lse = 'sparse'; %use sparse or full LU or lse with Hessenberg decomposition {'sparse', 'full','hess'}
Def.dgksTol = 1e-12; %orthogonality tolerance: norm(V'*V-I,'fro')<tol
        
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
if Opts.makeReal
    s0 = updateS0(s0);
end

% Tangential directions
if ~exist('Rt', 'var') || isempty(Rt)%   Compute block Krylov subspaces
    s0 = tangentialDirection(s0);
end

% Compute the Krylov subspaces
if hermite
    [V, Rsylv, W, Lsylv] = krylovSubspace(s0, q);
else
    [V, Rsylv] = krylovSubspace(s0, q);
end
    

%% ------------------ AUXILIARY FUNCTIONS --------------------------
% a) PRIMARY
    function [V, Rsylv, W, Lsylv] = krylovSubspace(s0, q)
    %   Calculate Krylov Subspace of s0
    %   Input:  s0:  Vector containing the expansion points
    %           q:   Original length of s0 with complex conjugated elements
    %   Output: V, W:  Krylov-Subspace of s0
    %           Rsylv, Lsylv: Sylvester matrix
        
    % preallocate memory
    V=zeros(length(B),q);
    Rsylv=zeros(size(B,2),q);
    if hermite 
        W = zeros(length(B),q); 
        Lsylv = zeros(size(C,1),q);
    end
    for jCol=1:length(s0)
        if hermite
            [V, Rsylv, W, Lsylv] = krylovDirection(jCol, s0, V, Rsylv, W, Lsylv);
        else
            [V, Rsylv] = krylovDirection(jCol, s0, V, Rsylv);
        end

        % split complex conjugate columns into real (->j) and imag (->j+length(s0c)/2
        if Opts.makeReal
            if hermite
                [V, Rsylv,W, Lsylv] = realSubspace(jCol, q, s0, V, Rsylv, W, Lsylv);
            else
                [V, Rsylv] = realSubspace(jCol, q, s0, V, Rsylv);
            end
        end

        if Opts.makeOrth
            if hermite
                [V, Rsylv, W, Lsylv] = gramSchmidt(jCol, V, Rsylv, W, Lsylv);
            else
                [V, Rsylv] = gramSchmidt(jCol, V, Rsylv);
            end
        end
    end

    %orthogonalize columns from imaginary components
    if Opts.makeOrth
        for jCol=length(s0)+1:q
            if hermite
                [V, Rsylv, W, Lsylv] = gramSchmidt(jCol, V, Rsylv, W, Lsylv);
            else
                [V, Rsylv] = gramSchmidt(jCol, V, Rsylv);
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
    
    if Opts.reorth
       switch Opts.reorth
           case 'gs' %reorthogonalized GS
                for jCol = 2:q        
                    if hermite
                        [V, Rsylv, W, Lsylv] = gramSchmidt(jCol, V, Rsylv, W, Lsylv);
                    else
                        [V, Rsylv] = gramSchmidt(jCol, V, Rsylv);
                    end
                end
           case 'qr' 
               [V,~] = qr(V,0); if hermite, [W,~] = qr(W,0); end
           case 'dgks'
           otherwise
               error('The orthogonalization chosen is incorrect or not implemented')
       end
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
    function [V, Rsylv, W, Lsylv] = gramSchmidt(jCol, V, Rsylv, W, Lsylv)
    %   Gram-Schmidt orthonormalization
    %   Input:  jCol:  Column to be treated
    %           V, W:  Krylov-Subspaces
    %           Rsylv, Lsylv: Sylvester matrices
    %   Output: V, W:  orthonormal basis of Krylov-Subspaces
    %           Rsylv, Lsylv: orthonormal Sylvester matrices
    
    if strcmp(Opts.reorth, 'dgks') && jCol>1
        % iterates standard gram-schmidt
        orthError=1;
        count=0;
        while(orthError>Opts.dgksTol)
            h=IP(V(:,1:jCol-1),V(:,jCol));
            V(:,jCol)=V(:,jCol)-V(:,1:jCol-1)*h;
            Rsylv(:,jCol)=Rsylv(:,jCol)-Rsylv(:,1:jCol-1)*h;
            if hermite
                h=IP(W(:,1:jCol-1),W(:,jCol));
                W(:,jCol)=W(:,jCol)-W(:,1:jCol-1)*h;
                Lsylv(:,jCol)=Lsylv(:,jCol)-Lsylv(:,1:jCol-1)*h;
            end
            orthError=norm(IP([V(:,1:jCol-1),V(:,jCol)/sqrt(IP(V(:,jCol),V(:,jCol)))],...
                [V(:,1:jCol-1),V(:,jCol)/sqrt(IP(V(:,jCol),V(:,jCol)))])-speye(jCol),'fro');
            if count>50 % if dgksTol is too small, Matlab can get caught in the while-loop
                error('Orthogonalization of the Krylov basis failed due to the given accuracy.');
            end
            count=count+1;
        end
    else
        for iCol=1:jCol-1
          h=IP(V(:,jCol),V(:,iCol));
          V(:,jCol)=V(:,jCol)-V(:,iCol)*h;
          Rsylv(:,jCol)=Rsylv(:,jCol)-h*Rsylv(:,iCol);
          if hermite
            h=IP(W(:,jCol),W(:,iCol));
            W(:,jCol)=W(:,jCol)-W(:,iCol)*h;
            Lsylv(:,jCol)=Lsylv(:,jCol)-h*Lsylv(:,iCol);
          end 
        end
    end     

    % normalize new basis vector
    h = sqrt(IP(V(:,jCol),V(:,jCol)));
    V(:,jCol)=V(:,jCol)/h;
    Rsylv(:,jCol) = Rsylv(:,jCol)/h;
    if hermite
        h = sqrt(IP(W(:,jCol),W(:,jCol)));
        W(:,jCol)=W(:,jCol)/h;
        Lsylv(:,jCol) = Lsylv(:,jCol)/h;
    end
    end

    function [V, Rsylv,W, Lsylv] = realSubspace(jCol, q, s0, V, Rsylv, W, Lsylv)
    %   Split Krylov direction into real and imaginary to create a real 
    %   Krylov subspace
    %   Input:  jCol:  Column to be treated
    %           nS0c:  Number of complex conjugated expansion points
    %           s0:    Vector containing the expansion points
    %           V, W:  Krylov-Subspaces
    %           Rsylv, Lsylv: Sylvester matrices
    %   Output: V, W:  orthonormal basis of Krylov-Subspaces
    %           Rsylv, Lsylv: orthonormal Sylvester matrices
    nS0c=q-length(s0);
    if ~isreal(s0(jCol))
        V(:,jCol+nS0c)=imag(V(:,jCol)); 
        V(:,jCol)=real(V(:,jCol));
        Rsylv(:,jCol+nS0c) = imag(Rsylv(:,jCol));
        Rsylv(:,jCol) = real(Rsylv(:,jCol));
        if hermite, 
            W(:,jCol+nS0c)=imag(W(:,jCol));
            W(:,jCol)=real(W(:,jCol)); 
            Lsylv(:,jCol+nS0c) = imag(Lsylv(:,jCol));
            Lsylv(:,jCol) = real(Lsylv(:,jCol));
        end
    end
    end
    
    function [V, Rsylv, W, Lsylv] = krylovDirection(jCol, s0, V, Rsylv, W, Lsylv)  
    %   Get new Krylov direction
    %   Input:  jCol:  Column to be treated
    %           s0:    Vector containing the expansion points
    %           V, W:  Krylov subspace
    %           Rsylv, Lsylv: Sylvester matrices
    %   Output: V, W:  Updated Krylov subspace
    %           Rsylv, Lsylv: updated Sylvester matrices
        
    % new basis vector
    tempV=B*Rt(:,jCol); newlu=1; newtan=1;
    Rsylv(:,jCol) = Rt(:,jCol);
    if hermite, tempW = C'*Lt(:,jCol); Lsylv(:,jCol) = Lt(:,jCol); end
    if jCol>1
        if s0(jCol)==s0(jCol-1)
            newlu=0;
            if Rt(:,jCol) == Rt(:,jCol-1)
                % Higher order moments, for the SISO and MIMO case
                newtan = 0;
                tempV = V(:,jCol-1); %overwrite
                Rsylv(:,jCol)=zeros(size(B,2),1);
                if hermite
                    tempW = W(:,jCol-1); 
                    Lsylv(:,jCol)=zeros(size(C,1),1); 
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
    end

    function [tempV, tempW] = lse(newlu, newtan, jCol, s0, tempV, tempW)
    %   Solve linear system of equations to obtain the new Krylov direction
    %   Input:  newlu: new lu decomposition required
    %           newtan: new tangential direction required
    %           jCol: Column to be treated
    %           s0: Vector containing the expansion points
    %           tempV, tempW: previous Krylov direction
    %   Output: tempV, tempW: new Krylov direction
    persistent R S L U a o;
    
    if isinf(s0(jCol)) %Realization problem (match Markov parameters)
        if newlu==0
            tempV=A*tempV;
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
                case 'full'
                    tempV = U\(L\tempV);
                case 'hess'
                    tempV = E\tempV;
            end
        else
            tempV = S*(R\(R'\(S'*tempV)));
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