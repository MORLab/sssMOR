function [sysr, V, W, Bb, Rsylv, Cb, Lsylv] = rk(sys, s0_inp, varargin)
% RK - Model Order Reduction by Rational Krylov
% ------------------------------------------------------------------
% Syntax:
%   RK(A)
%   [sysr, V, W] = RK(sys, s0_inp)
%   [sysr, V, W] = RK(sys, s0_inp, Rt)
% 
%   [sysr, V, W] = RK(sys, s0_inp, s0_out)
%   [sysr, V, W] = RK(sys, s0_inp, s0_out ,IP)
%   [sysr, V, W] = RK(sys, s0_inp, s0_out, Rt, Lt)
%   [sysr, V, W] = RK(sys, s0_inp, s0_out, Rt, Lt, IP)
%
% Inputs:       * sys: an sss-object containing the LTI system
%               * s0_inp: Expansion points for Input Krylov Subspace
%               * s0_out: Expansion points for Output Krylov Subspace
%               * IP: Inner product (optional)
% Outputs:      * sysr: reduced system
%               * V, W: Projection matrices spanning Krylov subspaces
% ------------------------------------------------------------------
% USAGE:  s0 may either be horizontal vectors containing the desired
% expansion points, e.g. [1 2 3] matches one moment about 1, 2 and 3,
% respectively. [1+1j 1-1j 5 5 5 5 inf inf] matches one moment about 1+1j,
% 1-1j, 4 moments about 5 and 2 Markov parameters.
%
% An alternative notation for s0 is a two-row matrix, containing the
% expansion points in the first row and their multiplicity in the second,
% e.g. [4 pi inf; 1 20 10] matches one moment about 4, 20 moments about pi
% and 10 Markov parameters.
%
% To perform one-sided RK, set s0_inp or s0_out to [], respectively.
%
% See also ARNOLDI, RK.
%
% ------------------------------------------------------------------
% REFERENCES:
% [1] Grimme (1997), Krylov Projection Methods for Model Reduction
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Heiko Panzer, Alessandro Castagnotto 
% Last Change:  26 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%%  Parsing
if nargin > 2
    if size(varargin{1}) == size(s0_inp);
        %usage: RK(sys, s0_inp, s0_out)
        s0_out = varargin{1};
        if nargin == 4
            %usage: RK(sys, s0_inp, s0_out ,IP)
            IP = varargin{2};
        elseif nargin > 4
            %usage: RK(sys, s0_inp, s0_out, Rt, Lt)
            Rt = varargin{2};
            Lt = varargin{3};
            if nargin == 6
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

if ~exist('IP', 'var'), 
    if ispd(sys.E) %assign IP to speed-up computations
        IP=@(x,y) (x'*sys.E*y); 
    else
        IP=@(x,y) (x'*y);
    end
end

if  (~exist('s0_inp', 'var') || isempty(s0_inp)) && ...
    (~exist('s0_out', 'var') || isempty(s0_out))
    error('No expansion points assigned.');
end

if exist('s0_inp', 'var')
    s0_inp = s0_vect(s0_inp);
else
    s0_inp = [];
end
if exist('s0_out', 'var')
    s0_out = s0_vect(s0_out);
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

%%  Computation
if isempty(s0_out)
    % input Krylov subspace

    % SISO Arnoldi
    [V,Rsylv] = arnoldi(sys.E, sys.A, sys.B, s0_inp, Rt, IP);
    W = V;
    sysr = sss(V'*sys.A*V, V'*sys.B, sys.C*V, sys.D, V'*sys.E*V);
    Bb = sys.B - sys.E*V*(sysr.E\sysr.B);
    Cb = []; Lsylv = [];  
    
elseif isempty(s0_inp)
    % output Krylov subspace
    
    % SISO Arnoldi
    [W,Lsylv] = arnoldi(sys.E', sys.A', sys.C', s0_out, Lt, IP);
    V = W;
    sysr = sss(W'*sys.A*W, W'*sys.B, sys.C*W, sys.D, W'*sys.E*W);
    Cb = sys.C - sysr.C/sysr.E*W'*sys.E;
    Bb = []; Rsylv = [];
else
    if all(s0_inp == s0_out) %use only 1 LU decomposition for V and W
        [V,Rsylv,W,Lsylv] = arnoldi(sys.E, sys.A, sys.B, sys.C,...
                            s0_inp,Rt, Lt, IP);
    else
        [V,Rsylv] = arnoldi(sys.E, sys.A, sys.B, s0_inp, Rt, IP);
        [W,Lsylv] = arnoldi(sys.E', sys.A', sys.C', s0_out, Lt, IP);
    end
    sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);

    if nargout > 3
        Bb = sys.B - sys.E*V*(sysr.E\sysr.B);
        Cb = sys.C - sysr.C/sysr.E*W'*sys.E;
    end
end

%----------- AUXILIARY --------------
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

%     s0 = cplxpair(s0);
    if size(s0,1)>size(s0,2)
        s0=transpose(s0);
    end



