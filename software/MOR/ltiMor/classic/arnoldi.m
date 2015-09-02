function [V,Ct,W,Bt] = arnoldi(E,A,b,varargin)
% ARNOLDI - Arnoldi algorithm using multiple expansion points
% ------------------------------------------------------------------
% [V,Ct]        = ARNOLDI(E,A,b,s0,IP)
% [V,Ct,W,Bt]   = ARNOLDI(E,A,b,c,s0,IP)
% Inputs:       * E,A,b,c: System matrices
%               * s0:    Vector of expansion points
%               * IP:    (opt.) function handle for inner product
% Outputs:      * V:    Orthonormal basis spanning the input Krylov subsp.
%               * Ct:   Right tangential directions of Sylvester Eq.
%               * W:    Orthonormal basis spanning the output Krylov subsp.
%               * Bt:   Left tangential directions of Sylvester Eq.
% ------------------------------------------------------------------
% USAGE:  This function is used to compute the matrix V spanning the 
% input Krylov subspace corresponding to E, A, b and s0 [1,2].
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
% This file is part of MORLab, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                  -> MORLab@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Heiko Panzer, Alessandro Castagnotto 
%               (a.castagnotto@tum.de)
% Last Change:  22 Jul 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%%  Parse input

if nargin == 4
    % usage: ARNOLDI(E,A,b,s0)
    s0 = varargin{1};
    hermite = 0; % same shifts for input and output Krylov?
elseif nargin == 5
    if size(varargin{1},2) == size(A,1)
        % usage: ARNOLDI(E,A,b,c,s0)
        c = varargin{1};
        s0 = varargin{2};
        hermite = 1;
    else
        % usage: ARNOLDI(E,A,b,s0,IP)
        s0 = varargin{1};
        IP = varargin{2};
        hermite = 0; % same shifts for input and output Krylov?
    end
elseif nargin == 6
    % usage: ARNOLDI(E,A,b,c,s0,IP)
    c = varargin{1};
    s0 = varargin{2};
    IP = varargin{3};
    hermite = 1;
else
    error('Wrong number of inputs')
end
    
if ~exist('IP', 'var') 
    if abs(condest(E))<Inf % 
        IP=@(x,y) (x'*E*y); 
    else
        IP=@(x,y) (x'*y); 
    end
end

if size(s0,1)>1
    error('s0 must be a vector containing the expansion points.')
end

q=length(s0); % order of the reduced model
reorth = 'gs'; %0, 'gs','qr'
% lseSol = 'lu'; %'lu', '\'


% % default execution options
% Def.reorth = 'gs'; %Orthogonalization algorithm (0, 'gs','qr')
% Def.lseSol = 'lu'; %Solution of LSE ('lu', '\')
% 
% % create the options structure
% if ~exist('Opts','var') || isempty(Opts)
%     Opts = Def;
% else
%     Opts = parseOpts(Opts,Def);
% end


%%  Compute the Krylov subspaces
% remove one element of complex pairs
k=find(imag(s0));
if ~isempty(k)
    s0c = cplxpair(s0(k));
    s0(k) = [];
    s0 = [s0 s0c(1:2:end)];
end

% preallocate memory
V=zeros(length(b),q);
Ct=eye(1,q); %**
if hermite, W = zeros(length(b),q); Bt = eye(1,q)';end

for jCol=1:length(s0)
    % new basis vector
    tempV=b; newlu=1; 
    Ct(jCol)=1; %**
    if hermite, tempW = c'; Bt(jCol) = 1; end ;
    if jCol>1
        if s0(jCol)==s0(jCol-1)
            tempV=V(:,jCol-1);
            newlu=0;
            Ct(jCol)=0; %**
            if hermite, tempW = W(:,jCol-1); Bt(jCol)=0; end
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
                    [L,U,p,o,S]=lu(sparse(E),'vector');
                else
                    rethrow(err);
                end
            end
        end
        if exist('U', 'var')
            tempV = invsolve(tempV,L,U,p,o,S);
        else
            tempV = S*(R\(R'\(S'*tempV)));
        end
    else %Rational Krylov
        if newlu==0
            tempV=E*tempV;
            if hermite, tempW = E'*tempW; end
        end
        if newlu==1
            % vector LU for sparse matrices
            [L,U,p,o,S]=lu(sparse(A-s0(jCol)*E),'vector');
        end
        % Solve the linear system of equations
        tempV(o,:) = U\(L\(S(:,p)\tempV)); %LU x(o,:) = S(:,p)\b 
        if hermite, tempW = (S(:,p)).'\(L.'\(U.'\(tempW(o,:)))); end %U'L'S(:,p) x = c'(o,:) 
    end 

    % split complex conjugate columns into real (->j) and imag (->j+length(s0c)/2
    if ~isreal(s0(jCol))
        V(:,jCol+length(s0c)/2)=imag(tempV); 
        tempV=real(tempV);
        if hermite, W(:,jCol+length(s0c)/2)=imag(tempW);tempW=real(tempW); end
    end

    % orthogonalize vectors
    for iCol=1:jCol-1
      h=IP(tempV,V(:,iCol));
      tempV=tempV-V(:,iCol)*h;
      Ct(jCol)=Ct(jCol)-h*Ct(iCol);
      if hermite
        h=IP(tempW,W(:,iCol));
        tempW=tempW-W(:,iCol)*h;
        Bt(jCol)=Bt(jCol)-h*Bt(iCol);
      end 
    end

    % normalize new basis vector
    h = sqrt(IP(tempV,tempV));
    V(:,jCol)=tempV/h;
    Ct(jCol) = Ct(jCol)/h;
    if hermite
        h = sqrt(IP(tempW,tempW));
        W(:,jCol)=tempW/h;
        Bt(jCol) = Bt(jCol)/h;
    end
   
end

%orthogonalize columns from imaginary components
for jCol=length(s0)+1:q
    tempV=V(:,jCol);
    if hermite, tempW=W(:,jCol);end
    for iCol=1:jCol-1
      h=IP(tempV, V(:,iCol));
      tempV=tempV-h*V(:,iCol);
      Ct(jCol) = Ct(jCol)-h*Ct(iCol);
      if hermite        
        h=IP(tempW, W(:,iCol));
        tempW=tempW-h*W(:,iCol);
        Bt(jCol) = Bt(jCol)-h*Bt(iCol);
      end
    end
    h = sqrt(IP(tempV,tempV));
    V(:,jCol)=tempV/h;
    Ct(jCol) = Ct(jCol)/h;
    if hermite
        h = sqrt(IP(tempW,tempW));
        W(:,jCol)=tempW/h;
        Bt(jCol) = Bt(jCol)/h;
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

