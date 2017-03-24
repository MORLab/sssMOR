function [sysr, V, W, B_, Sv, Rv, C_, Sw, Lw] = rk(sys, s0_inp, varargin)
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
%       [sysr, V, W]                         = RK(sys,s0_inp,...)
%		[sysr, V, W, B_, Sv, Rv]             = RK(sys,s0_inp,...)
%       [sysr, V, W, B_, Sv, Rv, C_, Sw, Lw] = RK(sys,s0_inp, s0_out, ...)
%		[sysr,...]                           = RK(sys, s0_inp, ..., Opts)
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
%                           [{true} / false]
%           -.(refer to help arnoldi for other options)
%           
%
% Output Arguments:
%       -sysr:              reduced system
%       -V,W:               projection matrices spanning Krylov subspaces
%       -B_,Sv,Rv:          resulting matrices of the input Sylvester
%                           equation (Rv is a (mxq) matrix, where q is the reduced order)
%       -C_,Sw,Lw:          resulting matrices of the output Sylvester
%                           equation (Lw is a (pxq) matrix, where q is the reduced order)
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
% Last Change:  06 Apr 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Parsing

% create the options structure
Def.real = true; %keep the projection matrices real?       

% use hess if sys is ssRed object
if isa(sys,'ssRed'), 
    Def.lse='hess'; 
    if isempty(sys.E), sys.E = eye(size(sys.A)); end %ssRed robust compatibility
else
    Def.lse = 'sparse'; 
end

if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
    
    Opts = parseOpts(Opts,Def);
else
    Opts = Def;
end       

% check usage and inputs
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
    error('sssMOR:rk:NoExpansionPoints','No expansion points assigned.');
end

if exist('s0_inp', 'var')
    s0_inp = shiftVec(s0_inp);
    % sort expansion points & tangential directions
    s0old = s0_inp;
    if Opts.real, 
        s0_inp = cplxpair(s0_inp); %make sure shifts can be paired 
    else
        s0_inp = sort(s0_inp);
    end
    if exist('Rt','var') && ~isempty(Rt)
        if size(Rt,2) ~= length(s0_inp),error('Inconsistent size of Rt');end
        
        [~,cplxSorting] = ismember(s0_inp,s0old); 
        Rt = Rt(:,cplxSorting);
    else
        Rt = [];
    end
clear s0old
else
    s0_inp = [];
end
if exist('s0_out', 'var')
    s0_out = shiftVec(s0_out);
        % sort expansion points & tangential directions
    s0old = s0_out;
    if Opts.real, 
        s0_out = cplxpair(s0_out); %make sure shifts can be paired 
    else
        s0_out = sort(s0_out);
    end
    if exist('Lt','var') && ~isempty(Lt)
        if size(Lt,2) ~= length(s0_out),error('Inconsistent size of Lt');end
        
        [~,cplxSorting] = ismember(s0_out,s0old); 
        Lt = Lt(:,cplxSorting);
    else
        Lt = [];
    end
else
    s0_out = [];
end
if length(s0_inp)> sys.n || length(s0_out)>sys.n
    error('sssMOR:arnoldi:reducedOrderExceedsOriginal','The desired reduced order exceeds the original order');
end

if ~isempty(s0_inp) && ~isempty(s0_out)
    % check if number of input/output expansion points matches
    if length(s0_inp) ~= length(s0_inp)
        error('Inconsistent length of expansion point vectors.');
    end
end
%%  Define execution variables
if ~exist('IP', 'var'), 
    IP=@(x,y) (x.'*y);
end
%%  Computation
if isempty(s0_out)
    % input Krylov subspace
    
    % SISO Arnoldi
    [V, Sv, Rv] = arnoldi(sys.E, sys.A, sys.B, s0_inp, Rt, IP, Opts);
    W = V;
    sysr = projectiveMor(sys,V,W);
    sysr.Name = sprintf('%s_%i_rk_inp',sys.Name,sysr.n);
    if nargout>3
        B_ = sys.B - sys.E*V*(sysr.E\sysr.B);
        C_ = []; Lw = [];  Sw=[];
    end
elseif isempty(s0_inp)
    % output Krylov subspace
    
    % SISO Arnoldi
    [W, Sw, Lw] = arnoldi(sys.E.', sys.A.', sys.C.', s0_out, Lt, IP, Opts);
    V = W;
    sysr = projectiveMor(sys,V,W);
    sysr.Name = sprintf('%s_%i_rk_out',sys.Name,sysr.n);
    if nargout>3
        C_ = sys.C - sysr.C/sysr.E*W.'*sys.E;
        B_ = []; Rv = []; Sv=[];
    end

else
    if all(s0_inp == s0_out) % use only 1 LU decomposition for V and W
        [V, Sv, Rv, W, Sw, Lw] = arnoldi(sys.E, sys.A, sys.B, sys.C,...
                            s0_inp,Rt, Lt, IP, Opts);
                        
        sysr = projectiveMor(sys,V,W);
        sysr.Name = sprintf('%s_%i_rk_herm',sys.Name,sysr.n);

    else
        [V, Sv, Rv] = arnoldi(sys.E, sys.A, sys.B, s0_inp, Rt, IP, Opts);
        [W, Sw, Lw] = arnoldi(sys.E.', sys.A.', sys.C.', s0_out, Lt, IP, Opts);
        if size(V,2)<size(W,2)
            V=[V,W(:,size(V,2)+1:size(W,2))];
        elseif size(V,2)>size(W,2)
            W=[W,V(:,size(W,2)+1:size(V,2))];
        end
        sysr = projectiveMor(sys,V,W);
        sysr.Name = sprintf('%s_%i_rk_2sided',sys.Name,sysr.n);
    end

    if nargout > 3
        B_ = sys.B - sys.E*V*(sysr.E\sysr.B);
        C_ = sys.C - sysr.C/sysr.E*W'*sys.E;
    end
end

%%  Storing additional parameters
%Stroring additional information about thr reduction in the object 
%containing the reduced model:
%   1. Define a new field for the Opts struct and write the information
%      that should be stored to this field
%   2. Adapt the method "parseParamsStruct" of the class "ssRed" in such a
%      way that the new defined field passes the check
Opts.originalOrder = sys.n;
if ~isfield(Opts,'orth') Opts.orth = '2mgs'; end
if ~isfield(Opts,'reorth') Opts.reorth = 0; end
if ~isfield(Opts,'lse') Opts.lse = 'sparse'; end
if ~isfield(Opts,'dgksTol') Opts.dgksTol = 1e-12; end
if ~isfield(Opts,'krylov') Opts.krylov = 0; end
Opts.IP = IP;
if ~exist('Rt','var') Opts.Rt = []; else Opts.Rt = Rt; end
if ~exist('Lt','var') Opts.Lt = []; else Opts.Lt = Lt; end
if ~exist('s0_inp','var') Opts.s0_inp = []; else Opts.s0_inp = s0_inp; end
if ~exist('s0_out','var') Opts.s0_out = []; else Opts.s0_out = s0_out; end

% Convert to ssRed-object
sysr = ssRed('rk',Opts,sysr);