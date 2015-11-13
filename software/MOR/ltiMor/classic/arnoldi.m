function [V,Rsylv,W,Lsylv] = arnoldi(E,A,B,varargin)
% ARNOLDI - Arnoldi algorithm for Krylov subspaces with multiple shifts
% 
% Syntax:
%       V					= ARNOLDI(E,A,B,s0)
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
%                       IP=@(x,y) (x'*I*y)
%
%       which requires E to be a positive definite matrix.
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
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  26 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Define execution parameters
if isstruct(varargin{end});
    %Options defined
    Opts = varargin{end};
    varargin = varargin(1:end-1);
end

Def.makeOrth = 1; %make orthogonal?
Def.makeReal = 1; %keep the projection matrices real?
Def.reorth = 'gs'; %gram schmidt reorthogonalization {0, 'gs','qr'}
        
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
        if length(varargin) == 5
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

m = size(B,2); if hermite, p = size(C,1); end
if exist('Rt','var') && ~isempty(Rt)
    if length(s0) ~= size(Rt,2),
        error('Rt must have the same columns as s0')
    end
    %   The reduced order is equivalent to the number of shifts
    q = length(s0);
else
    %   Block Krylov subspaces will be performed
    q = length(s0)*m;
end

if exist('Lt','var') && ~isempty(Lt)
    if length(s0) ~= size(Lt,2),
        error('Lt must have the same columns as s0')
    end
end

if Opts.makeReal
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
        nS0c = length(s0c); %number of complex shifts
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

nS0 = length(s0); %number of shifts for the computations

%   Tangential directions
if ~exist('Rt', 'var') || isempty(Rt)%   Compute block Krylov subspaces
    if m == 1; %SISO -> tangential directions are scalars
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
        
        s0old = s0; s0 = [];
        for iShift = 1:nS0
            s0 = [s0, s0old(iShift)*ones(1,m)];
        end
        Rt = repmat(speye(m,m),1,nS0);
        
        %update the number of shifts
        nS0     = length(s0);
        if exist('nS0c','var'), nS0c = m*nS0c; end
        
    end
    if hermite
        if m ~=p 
            error('Block Krylov for m~=p is not supported in arnoldi');
        else
            Lt = Rt;
        end
    end
end

%%  Define variables that might have not been passed to the function
%   IP
if ~exist('IP', 'var') 
    IP=@(x,y) (x'*y); 
end
%%  Compute the Krylov subspaces
% preallocate memory
V=zeros(length(B),q);
Rsylv=zeros(m,q);
if hermite, W = zeros(length(B),q); Lsylv = zeros(p,q);end
for jCol=1:nS0
    % new basis vector
    tempV=B*Rt(:,jCol); newlu=1; 
    Rsylv(:,jCol) = Rt(:,jCol);
    if hermite, tempW = C'*Lt(:,jCol); Lsylv(:,jCol) = Lt(:,jCol); end
    if jCol>1
        if s0(jCol)==s0(jCol-1)
            newlu=0;
            if Rt(:,jCol) == Rt(:,jCol-1)
                % Higher order moments, for the SISO and MIMO case
                newtan = 0;
                tempV = V(:,jCol-1); %overwrite
                Rsylv(:,jCol)=zeros(m,1);
                if hermite
                    tempW = W(:,jCol-1); 
                    Lsylv(:,jCol)=zeros(p,1); 
                end
            else
                newtan = 1;
            end
        end
    end
    
    if isinf(s0(jCol)) %Realization problem (match Markov parameters)
        if newlu==0
            tempV=A*tempV;
        end
        if newlu==1
            try
                % compute Cholesky factors of E
                clear L U a o S
                [R,~,S] = chol(sparse(E));
%                 R = chol(sparse(E));
            catch err
                if (strcmp(err.identifier,'MATLAB:posdef'))
                    % E is not pos. def -> use LU instead
                    [L,U,a,o,S]=lu(sparse(E),'vector');
                else
                    rethrow(err);
                end
            end
        end
        if exist('U', 'var')
            tempV(o,:) = U\(L\(S(:,a)\tempV)); %LU x(o,:) = S(:,a)\b 
        else
            tempV = S*(R\(R'\(S'*tempV)));
        end
    else %Rational Krylov
        if newlu==0
            if m==1 %SISO
                tempV=E*tempV;
                if hermite, tempW = E'*tempW; end
            elseif newtan==0
                % Tangential matching of higher order moments
                tempV=E*tempV;
                if hermite, tempW = E'*tempW; end
            end
        end
        if newlu==1
            % vector LU for sparse matrices
            [L,U,a,o,S]=lu(sparse(A-s0(jCol)*E),'vector');
        end
        % Solve the linear system of equations
        tempV(o,:) = U\(L\(S(:,a)\tempV)); %LU x(o,:) = S(:,a)\b 
        if hermite, tempW = (S(:,a)).'\(L.'\(U.'\(tempW(o,:)))); end %U'L'S(:,a) x = c'(o,:) 
    end 

    % split complex conjugate columns into real (->j) and imag (->j+length(s0c)/2
    if Opts.makeReal
        if ~isreal(s0(jCol))
            V(:,jCol+nS0c/2)=imag(tempV); 
            tempV=real(tempV);
            Rsylv(:,jCol+nS0c/2) = imag(Rsylv(:,jCol));
            Rsylv(:,jCol) = real(Rsylv(:,jCol));
            if hermite, 
                W(:,jCol+nS0c/2)=imag(tempW);tempW=real(tempW); 
                Lsylv(:,jCol+nS0c/2) = imag(Lsylv(:,jCol));
                Lsylv(:,jCol) = real(Lsylv(:,jCol));
            end
        end
    end

    if Opts.makeOrth
        % orthogonalize vectors
        for iCol=1:jCol-1
          h=IP(tempV,V(:,iCol));
          tempV=tempV-V(:,iCol)*h;
          Rsylv(:,jCol)=Rsylv(:,jCol)-h*Rsylv(:,iCol);
          if hermite
            h=IP(tempW,W(:,iCol));
            tempW=tempW-W(:,iCol)*h;
            Lsylv(:,jCol)=Lsylv(:,jCol)-h*Lsylv(:,iCol);
          end 
        end

        % normalize new basis vector
        h = sqrt(IP(tempV,tempV));
        V(:,jCol)=tempV/h;
        Rsylv(:,jCol) = Rsylv(:,jCol)/h;
        if hermite
            h = sqrt(IP(tempW,tempW));
            W(:,jCol)=tempW/h;
            Lsylv(:,jCol) = Lsylv(:,jCol)/h;
        end
    else
        V(:,jCol) = tempV; W(:,jCol) = tempW;
    end
end

%orthogonalize columns from imaginary components
if Opts.makeOrth
    for jCol=length(s0)+1:q
        tempV=V(:,jCol);
        if hermite, tempW=W(:,jCol);end
        for iCol=1:jCol-1
          h=IP(tempV, V(:,iCol));
          tempV=tempV-h*V(:,iCol);
          Rsylv(:,jCol) = Rsylv(:,jCol)-h*Rsylv(:,iCol);
          if hermite        
            h=IP(tempW, W(:,iCol));
            tempW=tempW-h*W(:,iCol);
            Lsylv(:,jCol) = Lsylv(:,jCol)-h*Lsylv(:,iCol);
          end
        end
        h = sqrt(IP(tempV,tempV));
        V(:,jCol)=tempV/h;
        Rsylv(:,jCol) = Rsylv(:,jCol)/h;
        if hermite
            h = sqrt(IP(tempW,tempW));
            W(:,jCol)=tempW/h;
            Lsylv(:,jCol) = Lsylv(:,jCol)/h;
        end
    end
end

%% Reorthogonalization
%{   
   Even modified Gram-Schmidt is not able to yield an orthonormal basis
   if the dimensions are high. Therefore, a reorthogonalization might be
   needed. On can choose to run modified GS again. From a theoretical 
   standpoint, this does not change the basis. However,
   numerically it is necessary to keep the numerics well behaved if the 
   reduced order is large
   The QR algorithm is much faster, however it does change the basis
%}
if Opts.reorth
   switch Opts.reorth
       case 'gs' %reorthogonalized GS
            for jCol = 2:q
                tempV = V(:,jCol);
                if hermite, tempW = W(:,jCol);end
                for iCol = 1:jCol-1
                     h=IP(tempV, V(:,iCol));
                     tempV=tempV-h*V(:,iCol);
                     Rsylv(:,jCol)=Rsylv(:,jCol)-h*Rsylv(:,iCol);
                     if hermite
                        h=IP(tempW, W(:,iCol));
                        tempW=tempW-h*W(:,iCol);
                        Lsylv(:,jCol)=Lsylv(:,jCol)-h*Lsylv(:,iCol);
                     end
                end
                h = sqrt(IP(tempV,tempV));
                V(:,jCol)=tempV/h;
                Rsylv(:,jCol) = Rsylv(:,jCol)/h;
                if hermite
                    h = sqrt(IP(tempW,tempW));
                    W(:,jCol)=tempW/h;
                    Lsylv(:,jCol) = Lsylv(:,jCol)/h;
                end
            end
       case 'qr' 
           [V,~] = qr(V,0); if hermite, [W,~] = qr(W,0); end
       otherwise
           error('The orthogonalization chosen is incorrect or not implemented')
   end
end

