function [V, Sv, Rv, W, Sw, Lw] = arnoldi(E,A,B,varargin)
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
%                           [{'sparse'} / 'full' / 'hess']
%           -.dgksTol:      tolerance for dgks orthogonalization
%                           [{1e-12} / positive float]
%           -.krylov:       standard or cascaded krylov basis
%                           [{'standardKrylov'} / 'cascadedKrylov']
%
% Output Arguments:
%       -V:        Orthonormal basis spanning the input Krylov subsp. 
%       -Sv:       Matrix of input Sylvester Eq. (1)
%       -Rv:       Right tangential directions of Sylvester Eq. (1), (mxq) matrix
%       -W:        Orthonormal basis spanning the output Krylov subsp.
%       -Sw:       Matrix of output Sylvester Eq. (2)
%       -Lw:       Left tangential directions of Sylvester Eq. (2), (pxq) matrix
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
% Last Change:  13 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
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
end

% Tangential directions
if ~exist('Rt', 'var') || isempty(Rt)%   Compute block Krylov subspaces
    s0 = tangentialDirection(s0);
end

% Compute the Krylov subspaces
if hermite
%     [V, Sv, Rv, W, Sw, Lw] = krylovSubspace(s0, q);
    [V, W, Sv, Rv, Sw, Lw] = solveLse(A,B,C,E,s0,Rt,Lt,IP,Opts);
    Sw=Sw.';
else
    [V, Sv, Rv] = solveLse(A,B,E,s0,Rt,IP,Opts);
end

%% ------------------ AUXILIARY FUNCTIONS --------------------------

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

end
