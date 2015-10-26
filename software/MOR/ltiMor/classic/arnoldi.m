function [V,Rsylv,W,Lsylv] = arnoldi(E,A,B,varargin)
% ARNOLDI - Arnoldi algorithm using multiple expansion points
% -------------------------------------------------------------------------
% [V,Ct]        = ARNOLDI(E,A,B,s0)
% [V,Ct]        = ARNOLDI(E,A,B,s0,IP)
% [V,Ct]        = ARNOLDI(E,A,B,s0,Rt)
% [V,Ct]        = ARNOLDI(E,A,B,s0,Rt,IP)
% [V,Ct,W,Bt]   = ARNOLDI(E,A,B,C,s0)
% [V,Ct,W,Bt]   = ARNOLDI(E,A,B,C,s0,IP)
% [V,Ct,W,Bt]   = ARNOLDI(E,A,B,C,s0,Rt,Lt)
% [V,Ct,W,Bt]   = ARNOLDI(E,A,B,C,s0,Rt,Lt,IP)
%
% Inputs:       * E,A,B,C: System matrices
%               * s0:    Vector of expansion points
%               * Rt,Lt:   (opt.) Matrix of right/left tangential directions
%               * IP:    (opt.) function handle for inner product
% Outputs:      * V:    Orthonormal basis spanning the input Krylov subsp.
%               * TODO Change the name of tangential directions!
%               * Rsylv:   Right tangential directions of Sylvester Eq.
%               * W:    Orthonormal basis spanning the output Krylov subsp.
%               * Lsylv:   Left tangential directions of Sylvester Eq.
% -------------------------------------------------------------------------
% USAGE:  This function is used to compute the matrix V spanning the 
% input Krylov subspace corresponding to E, A, b and s0 [1,2].
%
% s0 must be a vector of complex frequencies closed under conjugation. In
% case of MIMO systems, if matrices of tangential directions Rt (and Lt) are
% defined, they must have the same number of columns as the shifts, so that
% for each tangential direction it is clear to which shift it belongs. If
% not tangential directions are specified, then block Krylov subspaces are
% computed.
%
% NOTE that for MIMO systems, not all possible Krylov subspaces can be
% build with this function. For example, block Krylov subpspaces with
% multiplicities in the shifts are not supported so far.
%
% The columns of V build an orthonormal basis of the input Krylov 
% subspace. The orthogonalization is conducted using a reorthogonalized
% modified Gram-Schmidt procedure [3] with respect to the inner product
% defined in IP (optional). If no inner product is specified, then the
% elliptic product corresponding to E is chosen by default:
%                       IP=@(x,y) (x'*E*y)
% which requires E to be a positive definite matrix.
%
% See also RK.
%
% ------------------------------------------------------------------
% REFERENCES:
% [1] Grimme (1997), Krylov projection methods for model reduction
% [2] Antoulas (2005), Approximation of large-scale dynamical systems
%TODO: Reference for the duality between Krylov and Sylvester 
% [3] Giraud (2005), The loss of orthogonality in the Gram-Schmidt...
% ------------------------------------------------------------------
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                   -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Heiko Panzer, Alessandro Castagnotto 
%               (a.castagnotto@tum.de)
% Last Change:  26 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%%  Parse input

if nargin == 4
    % usage: ARNOLDI(E,A,B,s0)
    s0 = varargin{1};
    hermite = 0; % same shifts for input and output Krylov?
elseif nargin > 4
    %   Do the classification depending on the properties of the objects
    %   ARNOLDI(E,A,B,s0,...) or ARNOLDI(E,A,B,C,...)
    if size(varargin{1},2) == size(A,1)
        % usage: ARNOLDI(E,A,B,C,s0,...)
        hermite = 1;
        C = varargin{1};
        s0 = varargin{2};
        if nargin == 6
            % usage: ARNOLDI(E,A,B,C,s0,IP)
            IP = varargin{3};
        elseif nargin == 7
            % usage: ARNOLDI(E,A,B,C,s0,Rt,Lt)
            Rt = varargin{3};
            Lt = varargin{4};
        elseif nargin == 8
            % usage: ARNOLDI(E,A,B,C,s0,Rt,Lt,IP)
            Rt = varargin{3};
            Lt = varargin{4};
            IP = varargin{5};
        end
    else
        % usage: ARNOLDI(E,A,B,s0,...)
        hermite = 0;
        s0 = varargin{1};
        if nargin == 5
            if size(varargin{2},2) == size(s0,2)
                % usage: ARNOLDI(E,A,B,s0,Rt)
                Rt = varargin{2};
            else   
                % usage: ARNOLDI(E,A,B,s0,IP)
                IP = varargin{2};
            end
        else
            % usage: ARNOLDI(E,A,b,s0,Rt,IP)
            Rt = varargin{3};
            IP = varargin{4};
        end
    end
end

if size(s0,1)>1
    error('s0 must be a vector containing the expansion points.')
end

m = size(B,2);

if exist('Rt','var')
    if length(s0) ~= size(Rt,2),
        error('Rt must have the same columns as s0')
    end
    %   The reduced order is equivalent to the number of shifts
    q = length(s0);
else
    %   Block Krylov subspaces will be performed
    q = length(s0)*m;
end

if exist('Lt','var')
    if length(s0) ~= size(Lt,2),
        error('Lt must have the same columns as s0')
    end
end

%%  Define variables that might have not been passed to the function

%   IP
if ~exist('IP', 'var') 
    if abs(condest(E))<Inf % 
        IP=@(x,y) (x'*E*y); 
    else
        IP=@(x,y) (x'*y); 
    end
end


% remove one element of complex pairs (must be closed under conjugation)
k=find(imag(s0));
if ~isempty(k)
    s0c = cplxpair(s0(k));
    s0(k) = [];
    s0 = [s0 s0c(1:2:end)];
end

nS0 = length(s0); %number of shifts for the computations
nS0c = length(s0c); %number of complex shifts

%   Tangential directions
if ~exist('Rt', 'var') %   Compute block Krylov subspaces
    if m == 1; %SISO -> tangential directions are scalars
        Rt = ones(1,nS0);
    else %MIMO -> fill up s0 and define tangential blocks
        
        % tangential matching of higher order moments not implemented so
        % far! Therefore, if two shifts are the same, give an error
        if any(diff(s0)==0)
            error(['tangential matching of higher order moments with block'...
                'Krylov not implemented yet']);
        end
        
        s0old = s0; s0 = [];
        for iShift = 1:nS0
            s0 = [s0, s0old(iShift)*ones(1,m)];
        end
        Rt = repmat(speye(m,m),1,nS0);
        
        %update the number of shifts
        nS0     = length(s0);
        nS0c    = m*nS0c; 
        
    end
    if hermite
        p = size(C,1); 
        if m ~=p 
            error('Block Krylov for m~=p is not supported in arnoldi');
        else
            Lt = Rt;
        end
    end
else % Rt (and Lt) were specified
    % make sure they have the right length
end

reorth = 'gs'; %0, 'gs','qr'
% lseSol = 'lu'; %'lu', '\'


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
                if m == 1 %SISO
                    tempV = V(:,jCol-1); %overwrite
                    Rsylv(:,jCol)=zeros(m,1);
                else
                    error('Higher order moments for MIMO Systems not implemented yet');
                end
            end
            if hermite, tempW = W(:,jCol-1); Lsylv(:,jCol)=zeros(p,1); end
        end
    end
    
    if isinf(s0(jCol)) %Realization problem (match Markov parameters)
        if newlu==0
            tempV=A*tempV;
        end
        if newlu==1
            try
                % compute Cholesky factors of E
                clear L U o p S
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
            else
                % MIMO with different tangential direction: do nothing
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
    if ~isreal(s0(jCol))
        V(:,jCol+nS0c/2)=imag(tempV); 
        tempV=real(tempV);
        if hermite, W(:,jCol+nS0c/2)=imag(tempW);tempW=real(tempW); end
    end

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
   
end

%orthogonalize columns from imaginary components
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
if reorth
   switch reorth
       case 'gs' %reorthogonalized GS
            for jCol = 2:q
                tempV = V(:,jCol);
                if hermite, tempW = W(:,jCol);end
                for iCol = 1:jCol-1
                     h=IP(tempV, V(:,iCol));
                     tempV=tempV-h*V(:,iCol);
                     if hermite
                        h=IP(tempW, W(:,iCol));
                        tempW=tempW-h*W(:,iCol);
                     end
                end
                h = sqrt(IP(tempV,tempV));
                V(:,jCol)=tempV/h;
                if hermite
                    h = sqrt(IP(tempW,tempW));
                    W(:,jCol)=tempW/h;
                end
            end
       case 'qr' 
           V = qr(V,0); if hermite, W = qr(W,0); end
       otherwise
           error('The orthogonalization chosen is incorrect or not implemented')
   end
end

