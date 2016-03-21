function [sysr, V, W, Bb, SRsylv, Rsylv, Cb, SLsylv, Lsylv] = rk(sys, s0_inp, varargin)
% RK - Model Order Reduction by Rational Krylov
%
% Syntax:
%       sysr = RK(sys, s0_inp)
%       sysr = RK(sys, s0_inp, Rt)
%       sysr = RK(sys, [], s0_out)
%       sysr = RK(sys, [], s0_out, [], Lt)
%       sysr = RK(sys, s0_inp, s0_out)
%       sysr = RK(sys, s0_inp, s0_out ,IP)
%       sysr = RK(sys, s0_inp, s0_out, Rt, Lt)
%       sysr = RK(sys, s0_inp, s0_out, Rt, Lt, IP)
%
%       [sysr, V, W]                                        = RK(sys,s0_inp,...)
%		[sysr, V, W, Bb, SRsylv, Rsylv]                     = RK(sys,s0_inp,...)
%       [sysr, V, W, Bb, SRsylv, Rsylv, Cb, SLsyslv, Lsylv]	= RK(sys,s0_inp, s0_out, ...)
%		[sysr,...]                                          = RK(sys, s0_inp, ..., Opts)
%
% Description:
%       Reduction by Rational Krylov subspace methods. 
%
%       This function computes input (and output) Krylov subspaces
%       corresponding to the shifts s0_inp and s0_out. For MIMO systems,
%       tangential directions Rt and Lt can be defined. Otherwise, block
%       Krylov subspaces will be computed.
%
%       s0 may either be horizontal vectors containing the desired
%       expansion points, e.g. [1 2 3] matches one moment about 1, 2 and 3,
%       respectively. [1+1j 1-1j 5 5 5 5 inf inf] matches one moment about 1+1j,
%       1-1j, 4 moments about 5 and 2 Markov parameters.
% 
%       An alternative notation for s0 is a two-row matrix, containing the
%       expansion points in the first row and their multiplicity in the second,
%       e.g. [4 pi inf; 1 20 10] matches one moment about 4, 20 moments about pi
%       and 10 Markov parameters.
% 
%       To perform one-sided RK, set s0_inp or s0_out to [], respectively.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys:       an sss-object containing the LTI system
%       -s0_inp:    expansion points for input Krylov subspace
%
%       *Optional Input Arguments:*
%       -s0_out:            expansion points for output Krylov subspace
%       -Rt/Lt:             right/left tangential directions
%       -IP:                inner product (optional)
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
%       -sysr:              reduced system
%       -V,W:               projection matrices spanning Krylov subspaces
%       -Bb,SRsylv,Rsylv:   resulting matrices of the input Sylvester equation
%       -Cb,SLsylv,Lsylv:   resulting matrices of the output Sylvester equation
%
% Examples:
%       This code reduces the benchmark model build by orthogonal
%       projection using the input Krylov subspace corresponding to the
%       shifts s0 = [1 1 2 2 2] (i.e. matching two moments about 1 and
%       three moments about 2)
%
%> sys = loadSss('building');
%> s0  = [1 1 2 2 2];
%> sysr = rk(sys,s0);
%
%       The same result can be achieved by specyfing the shifts as a matrix
%       with two rows:
%
%> s0 = [1 2; 2 3];
%> sysr = rk(sys,s0);
%
%       For two-sided reduction, specify shifts for the output Krylov
%       subspace
%
%> sysr = rk(sys,s0,s0);
%
%       To compute an orthogonal projection by output Krylov subspace, set
%       the input shifts to []
%
%> sysr = rk(sys,[],s0);
%
%       For MIMO systems, specify tangential directions
%
%> sys = loadSss('CDplayer'); n = 5;
%> s0 = rand(1,n); Rt = rand(sys.m,n); Lt = rand(sys.p,n);
%> sysr = rk(sys, s0, s0, Rt, Lt);
%
%       If no tangential directions are specified, block Krylov reduction
%       is conducted.
%
%> sysr = rk(sys, s0, s0); 
%
%//Note: In the block Krylov case, the reduced order is in general higher
%       than the lenght of s0.
%
% See Also: 
%       arnoldi, irka
%
% References:
%       * *[1] Grimme (1997)*, Krylov projection methods for model reduction
%       * *[2] Antoulas (2010)*, Interpolatory model reduction of
%              large-scale dynamical Systems
%       * *[3] Beattie et al. (2014)*, Model reduction by rational interpolation 
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
% Last Change:  09 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Parsing
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
else
    Opts = struct();
end

if ~isempty(varargin)
    if isempty(s0_inp) || all(size(varargin{1}) == size(s0_inp));
        %usage: RK(sys, s0_inp, s0_out)
        s0_out = varargin{1};
        if length(varargin) == 2
            %usage: RK(sys, s0_inp, s0_out ,IP)
            IP = varargin{2};
        elseif length(varargin) > 2
            %usage: RK(sys, s0_inp, s0_out, Rt, Lt)
            Rt = varargin{2};
            Lt = varargin{3};
            if length(varargin) == 4
                %usage: RK(sys, s0_inp, s0_out, Rt, Lt, IP)
                IP = varargin{4};
            end
        end
    elseif size(varargin{1},1) == size(sys.B,2);
        %usage: RK(sys, s0_inp, Rt)
        Rt = varargin{1};
    else
        error('Input not compatible with current rk implementation');
    end
end

%%  Check the inputs
if  (~exist('s0_inp', 'var') || isempty(s0_inp)) && ...
    (~exist('s0_out', 'var') || isempty(s0_out))
    error('No expansion points assigned.');
end

if exist('s0_inp', 'var')
    s0_inp = s0_vect(s0_inp);
    % sort expansion points & tangential directions
    s0old = s0_inp;
    s0_inp = cplxpair(s0_inp);
    if exist('Rt','var') && ~isempty(Rt)
        [~,cplxSorting] = ismember(s0_inp,s0old); 
        Rt = Rt(:,cplxSorting);
    end
clear s0old
else
    s0_inp = [];
end
if exist('s0_out', 'var')
    s0_out = s0_vect(s0_out);
        % sort expansion points & tangential directions
    s0old = s0_out;
    s0_out = cplxpair(s0_out);
    if exist('Lt','var') && ~isempty(Lt)
        [~,cplxSorting] = ismember(s0_out,s0old); 
        Lt = Lt(:,cplxSorting);
    end
else
    s0_out = [];
end

if exist('Rt', 'var')
    if size(Rt,2) ~= length(s0_inp),error('Inconsistent size of Rt');end
else
    Rt = [];
end
if exist('Lt', 'var')
    if size(Lt,2) ~= length(s0_out),error('Inconsistent size of Lt');end
else
    Lt = [];
end

if ~isempty(s0_inp) && ~isempty(s0_out)
    % check if number of input/output expansion points matches
    if length(s0_inp) ~= length(s0_inp)
        error('Inconsistent length of expansion point vectors.');
    end
end
%%  Define execution variables
if ~exist('IP', 'var'), 
    IP=@(x,y) (x'*y); %seems to be better conditioned
end
%%  Computation
if isempty(s0_out)
    % input Krylov subspace
    
    % SISO Arnoldi
    [V, SRsylv, Rsylv] = arnoldi(sys.E, sys.A, sys.B, s0_inp, Rt, IP, Opts);
    W = V;
    sysr = sss(V'*sys.A*V, V'*sys.B, sys.C*V, sys.D, V'*sys.E*V);
    sysr.Name = sprintf('%s_rk_inp',sys.Name);
    if nargout>3
        Bb = sys.B - sys.E*V*(sysr.E\sysr.B);
        Cb = []; Lsylv = [];  SLsylv=[];
    end
elseif isempty(s0_inp)
    % output Krylov subspace
    
    % SISO Arnoldi
    [W, SLsylv, Lsylv] = arnoldi(sys.E', sys.A', sys.C', s0_out, Lt, IP, Opts);
    V = W;
    sysr = sss(W'*sys.A*W, W'*sys.B, sys.C*W, sys.D, W'*sys.E*W);
    sysr.Name = sprintf('%s_rk_out',sys.Name);
    if nargout>3
        Cb = sys.C - sysr.C/sysr.E*W'*sys.E;
        Bb = []; Rsylv = []; SRsylv=[];
    end

else
    if all(s0_inp == s0_out) %use only 1 LU decomposition for V and W
        [V, SRsylv, Rsylv, W, SLsylv, Lsylv] = arnoldi(sys.E, sys.A, sys.B, sys.C,...
                            s0_inp,Rt, Lt, IP, Opts);
                        
        sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);
        sysr.Name = sprintf('%s_rk_herm',sys.Name);

    else
        [V, SRsylv, Rsylv] = arnoldi(sys.E, sys.A, sys.B, s0_inp, Rt, IP, Opts);
        [W, SLsylv, Lsylv] = arnoldi(sys.E', sys.A', sys.C', s0_out, Lt, IP, Opts);
        sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);
        sysr.Name = sprintf('%s_rk_2sided',sys.Name);
    end


    if nargout > 3
        Bb = sys.B - sys.E*V*(sysr.E\sysr.B);
        Cb = sys.C - sysr.C/sysr.E*W'*sys.E;
    end
end


%% ----------- AUXILIARY --------------
function s0=s0_vect(s0)
    % change two-row notation to vector notation
    if size(s0,1)==2
        temp=zeros(1,sum(s0(2,:)));
        for j=1:size(s0,2)
            k=sum(s0(2,1:(j-1))); k=(k+1):(k+s0(2,j));
            temp(k)=s0(1,j)*ones(1,s0(2,j));
        end
        s0=temp;
    end

    if size(s0,1)>size(s0,2)
        s0=transpose(s0);
    end
