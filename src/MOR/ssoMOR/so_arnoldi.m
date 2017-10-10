function [V, R] = so_arnoldi(sys, s0, Rt, IP, Opts)
% ARNOLDI - Arnoldi algorithm for Krylov subspaces with multiple shifts
% 
% Syntax:
%       [V,R] = ARNOLDI(sys, s0, Rt, IP, Opts)
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  01 Apr 2017
% Copyright (c) 2016,2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Define execution parameters
Def.real    = true; %keep the projection matrices real?
Def.orth    = '2mgs'; %orthogonalization after every direction {0,'dgks','mgs','2mgs'}
Def.reorth  = 0; %reorthogonaliation at the end {0, 'mgs', 'qr'}
Def.lse     = 'sparse'; %use sparse or full LU or lse with Hessenberg decomposition {'sparse', 'full','hess'}
Def.dgksTol = 1e-12; %orthogonality tolerance: norm(V'*V-I,'fro')<tol
Def.krylov  = 'standardKrylov'; %standard or cascaded krylov basis (only for siso) {'standardKrylov','cascadedKrylov'}

if  exist('Opts','var') && ~isempty(Opts)
    Opts = parseOpts(Opts,Def);
else
    Opts = Def;
end          
 
if islogical(Opts.orth) && Opts.orth
    Opts.orth = Def.orth; %set default orthogonalization method
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
    q = length(s0)*size(sys.B,2);
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
end

% Tangential directions
if ~exist('Rt', 'var') || isempty(Rt)%   Compute block Krylov subspaces
    s0 = tangentialDirection(s0);
end

% preallocate memory
V   = zeros(length(sys.B),q);
R   = zeros(size(sys.B,2),q);

% Compute the Krylov subspaces
for jCol=1:length(s0)
    [V,~,Rsylv] = updateV(jCol, V, sys, s0, Rt, Opts);    
    R(:,jCol)   = Rsylv*Rt(:,jCol);

    % split complex conjugate columns into real (->j) and imag (->j+length(s0c)/2
    if Opts.real
        [V, R] = realSubspace(jCol, q, s0, V, R);
    end

    if Opts.orth
        [V, TRv]    = gramSchmidt(jCol, V);
        R           = R*TRv;
    end
end

%orthogonalize columns from imaginary components
if Opts.orth
    for jCol=length(s0)+1:q
        [V, TRv]        = gramSchmidt(jCol, V);
        R=R*TRv;
    end
end

% transform hess back
if exist('Z','var') && ~isempty(Z)
    V=Z*V;
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
            R=R*TRv;
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
       R=R*Rinv;
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

%a) PRIMARY
function [s0] = tangentialDirection(s0)
%   Update s0 and define Rt for calculation of tangential directions
%   Input:  s0: Vector containing the expansion points
%   Output: s0: Vector containing the updated expansion points
if size(sys.B,2) == 1; %SISO -> tangential directions are scalars
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
            for j=1:size(sys.B,2)
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
function [V,s0,Rsylv] = updateV(jCol, V, sys, s0, Rt,Opts)
    persistent L U a o S
    
    if jCol == 1
        [L,U,a,o,S] = lu(s0(jCol)^2*sys.M +s0(jCol)*sys.D + sys.K ,'vector');
        V(o,jCol) = U\(L\(S(:,a)\(sys.B*Rt(:,jCol))));     %LU x(o,:) = S(:,a)\b
        Rsylv=eye(size(sys.B,2));
    elseif jCol == 2
        if s0(jCol) == s0(jCol-1)
            V(o,jCol) = -U\(L\(S(:,a)\((2*s0(jCol)*sys.M+sys.D)*V(:,jCol-1))));
            Rsylv=zeros(size(sys.B,2));
        else
            [L,U,a,o,S] = lu(s0(jCol)^2*sys.M +s0(jCol)*sys.D + sys.K ,'vector');
            V(o,jCol) = U\(L\(S(:,a)\(sys.B*Rt(:,jCol))));     %LU x(o,:) = S(:,a)\b
            Rsylv=eye(size(sys.B,2));
        end
    else %jCol >= 3
        nRep = sum(s0(jCol)==s0(jCol-[1:2]));
        switch nRep
            case 0
                [L,U,a,o,S] = lu(s0(jCol)^2*sys.M +s0(jCol)*sys.D + sys.K ,'vector');
                V(o,jCol) = U\(L\(S(:,a)\(sys.B*Rt(:,jCol))));     %LU x(o,:) = S(:,a)\b
                Rsylv=eye(size(sys.B,2));
            case 1
                V(o,jCol) = -U\(L\(S(:,a)\((2*s0(jCol)*sys.M+D)*V(:,jCol-1))));
                Rsylv=zeros(size(sys.B,2));
            case 2
                V(o,jCol) = -U\(L\(S(:,a)\((2*s0(jCol)*sys.M+sys.D)*V(:,jCol-1)+sys.M*V(:,jCol-2))));
                Rsylv=zeros(size(sys.B,2));
        end
    end
end

%b) SECONDARY
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
end

function [V, Rv, W, Sw, Lw] = realSubspace(jCol, q, s0, V, Rv, W, Sw, Lw)
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
    V(:,jCol+nS0c)  = imag(V(:,jCol)); 
    V(:,jCol)       = real(V(:,jCol));
    Rv(:,jCol+nS0c) = imag(Rv(:,jCol));
    Rv(:,jCol)      = real(Rv(:,jCol));
end
end

end
